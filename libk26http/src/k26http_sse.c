#include "k26http.h"

#include <curl/curl.h>
#include <stdlib.h>
#include <string.h>

/* Minimal SSE client: buffers WRITEFUNCTION chunks, splits on blank
 * line (\n\n), extracts `data:` payloads, and delivers each payload to
 * the user callback. Anything else in the frame (event:, id:, retry:,
 * comments starting with ':') is ignored — callers that need those
 * metadata fields can layer their own parser on top of the raw bytes
 * they see in their own CURLOPT if they ever need to.
 *
 * Framing is lossless across chunk boundaries: we only pop complete
 * frames (ones with a trailing \n\n) out of the raw buffer per write. */

typedef struct {
    char  *data;
    size_t len;
    size_t cap;
} SB;

static bool sb_append(SB *b, const char *s, size_t n)
{
    if (b->len + n + 1 > b->cap) {
        size_t ncap = b->cap ? b->cap : 2048;
        while (ncap < b->len + n + 1) ncap *= 2;
        char *p = realloc(b->data, ncap);
        if (!p) return false;
        b->data = p; b->cap = ncap;
    }
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
    return true;
}

typedef struct {
    K26HTTPSSECb cb;
    void        *user;
    SB           raw;
    bool         aborted;
    bool       (*should_cancel)(void *);
    void        *cancel_user;
} SSECtx;

/* Deliver one logical frame (multi-line key:value text) to the callback,
 * stripping the SSE framing and presenting just the `data:` payload. If
 * the frame contains multiple `data:` lines, they're concatenated with
 * '\n' — matching the SSE spec. */
static void deliver_frame(SSECtx *ctx, char *frame)
{
    /* Collect data payloads. */
    SB payload = {0};
    char *line = frame;
    while (line && *line) {
        char *nl = strchr(line, '\n');
        if (nl) *nl = '\0';
        if (line[0] == ':' || line[0] == '\0') {
            /* comment / blank — skip */
        } else if (strncmp(line, "data:", 5) == 0) {
            const char *p = line + 5;
            if (*p == ' ') p++;
            if (payload.len) sb_append(&payload, "\n", 1);
            sb_append(&payload, p, strlen(p));
        }
        if (!nl) break;
        line = nl + 1;
    }
    if (payload.data && payload.len > 0) {
        if (!ctx->cb(payload.data, ctx->user)) ctx->aborted = true;
    }
    free(payload.data);
}

static size_t sse_write(char *ptr, size_t size, size_t nmemb, void *user)
{
    SSECtx *ctx = user;
    if (ctx->aborted) return 0;
    if (ctx->should_cancel && ctx->should_cancel(ctx->cancel_user)) {
        ctx->aborted = true;
        return 0;
    }

    size_t n = size * nmemb;
    if (!sb_append(&ctx->raw, ptr, n)) return 0;

    /* Pop every complete frame (\n\n separated). */
    for (;;) {
        char *split = strstr(ctx->raw.data, "\n\n");
        if (!split) break;
        size_t flen = (size_t)(split - ctx->raw.data);
        char *frame = malloc(flen + 1);
        if (!frame) return 0;
        memcpy(frame, ctx->raw.data, flen);
        frame[flen] = '\0';

        size_t remain = ctx->raw.len - flen - 2;
        memmove(ctx->raw.data, split + 2, remain + 1);
        ctx->raw.len = remain;

        deliver_frame(ctx, frame);
        free(frame);
        if (ctx->aborted) return 0;
    }
    return n;
}

K26HTTPStatus k26http_sse(const K26HTTPRequest *req, K26HTTPSSECb cb, void *user)
{
    if (!req || !req->url || !cb) return K26HTTP_ERR_INTERNAL;

    CURL *h = curl_easy_init();
    if (!h) return K26HTTP_ERR_INTERNAL;

    struct curl_slist *hdr = NULL;
    /* SSE servers expect `Accept: text/event-stream`. We unconditionally
     * prepend it so callers don't have to. */
    hdr = curl_slist_append(hdr, "Accept: text/event-stream");
    if (req->headers) {
        for (size_t i = 0; req->headers[i]; i++) {
            struct curl_slist *next = curl_slist_append(hdr, req->headers[i]);
            if (!next) { curl_slist_free_all(hdr); curl_easy_cleanup(h); return K26HTTP_ERR_INTERNAL; }
            hdr = next;
        }
    }

    SSECtx ctx = {
        .cb = cb, .user = user,
        .should_cancel = req->should_cancel, .cancel_user = req->user,
    };

    curl_easy_setopt(h, CURLOPT_URL, req->url);
    curl_easy_setopt(h, CURLOPT_USERAGENT,
                     req->user_agent ? req->user_agent : "libk26http/1.0");
    curl_easy_setopt(h, CURLOPT_HTTPHEADER, hdr);
    curl_easy_setopt(h, CURLOPT_WRITEFUNCTION, sse_write);
    curl_easy_setopt(h, CURLOPT_WRITEDATA, &ctx);
    curl_easy_setopt(h, CURLOPT_TIMEOUT_MS,
                     (long)(req->timeout_ms > 0 ? req->timeout_ms : 90000));
    curl_easy_setopt(h, CURLOPT_CONNECTTIMEOUT_MS,
                     (long)(req->connect_timeout_ms > 0 ? req->connect_timeout_ms : 5000));
    curl_easy_setopt(h, CURLOPT_NOSIGNAL, 1L);

    if (req->method && strcmp(req->method, "POST") == 0) {
        curl_easy_setopt(h, CURLOPT_POST, 1L);
        curl_easy_setopt(h, CURLOPT_POSTFIELDS, req->body ? req->body : "");
        curl_easy_setopt(h, CURLOPT_POSTFIELDSIZE,
                         (long)(req->body ? strlen(req->body) : 0));
    } else {
        curl_easy_setopt(h, CURLOPT_HTTPGET, 1L);
    }

    CURLcode cc = curl_easy_perform(h);
    long code = 0;
    curl_easy_getinfo(h, CURLINFO_RESPONSE_CODE, &code);
    curl_slist_free_all(hdr);
    curl_easy_cleanup(h);
    free(ctx.raw.data);

    if (ctx.aborted) return K26HTTP_ERR_CANCELLED;
    if (cc == CURLE_OPERATION_TIMEDOUT) return K26HTTP_ERR_TIMEOUT;
    if (cc == CURLE_COULDNT_RESOLVE_HOST
     || cc == CURLE_COULDNT_RESOLVE_PROXY
     || cc == CURLE_COULDNT_CONNECT
     || cc == CURLE_INTERFACE_FAILED) return K26HTTP_ERR_OFFLINE;
    if (cc != CURLE_OK) return K26HTTP_ERR_INTERNAL;

    if (code == 429)                return K26HTTP_ERR_429;
    if (code >= 200 && code < 300)  return K26HTTP_OK;
    if (code >= 400 && code < 500)  return K26HTTP_ERR_4XX;
    if (code >= 500 && code < 600)  return K26HTTP_ERR_5XX;
    return K26HTTP_ERR_INTERNAL;
}
