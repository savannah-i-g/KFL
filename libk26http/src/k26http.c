#include "k26http.h"

#include <curl/curl.h>
#include <stdlib.h>
#include <string.h>

#define K26HTTP_BUF_INITIAL_CAP 4096
#define K26HTTP_BUF_HARD_CAP    (16u * 1024u * 1024u)
#define K26HTTP_UA_DEFAULT      "libk26http/1.0"

/* -------- buffer ----------------------------------------------- */

void k26http_buf_init(K26HTTPBuf *b)
{
    if (!b) return;
    b->data = NULL; b->len = 0; b->cap = 0;
}

void k26http_buf_free(K26HTTPBuf *b)
{
    if (!b) return;
    free(b->data);
    b->data = NULL; b->len = 0; b->cap = 0;
}

bool k26http_buf_append(K26HTTPBuf *b, const void *src, size_t n)
{
    if (!b || (!src && n)) return false;
    size_t need = b->len + n + 1;
    if (need > K26HTTP_BUF_HARD_CAP) return false;
    if (need > b->cap) {
        size_t cap = b->cap ? b->cap : K26HTTP_BUF_INITIAL_CAP;
        while (cap < need) cap *= 2;
        if (cap > K26HTTP_BUF_HARD_CAP) cap = K26HTTP_BUF_HARD_CAP;
        char *p = realloc(b->data, cap);
        if (!p) return false;
        b->data = p; b->cap = cap;
    }
    memcpy(b->data + b->len, src, n);
    b->len += n;
    b->data[b->len] = '\0';
    return true;
}

/* -------- libcurl global refcount ------------------------------ */

static int g_curl_refs = 0;

void k26http_lib_ref(void)
{
    if (g_curl_refs == 0) curl_global_init(CURL_GLOBAL_DEFAULT);
    g_curl_refs++;
}
void k26http_lib_unref(void)
{
    if (--g_curl_refs == 0) curl_global_cleanup();
    if (g_curl_refs < 0)    g_curl_refs = 0;
}

/* -------- one-shot --------------------------------------------- */

typedef struct {
    K26HTTPBuf *out;
    size_t      max_bytes;
    bool        truncated;
    bool      (*should_cancel)(void *);
    void       *user;
} WriteCtx;

static size_t write_cb(char *ptr, size_t size, size_t nmemb, void *userp)
{
    WriteCtx *wc = userp;
    if (wc->should_cancel && wc->should_cancel(wc->user)) return 0;
    size_t n = size * nmemb;
    if (wc->max_bytes && wc->out->len + n > wc->max_bytes) {
        size_t room = (wc->out->len < wc->max_bytes)
                      ? (wc->max_bytes - wc->out->len) : 0;
        if (room) k26http_buf_append(wc->out, ptr, room);
        wc->truncated = true;
        return 0;
    }
    if (!k26http_buf_append(wc->out, ptr, n)) return 0;
    return n;
}

static K26HTTPStatus map_http_status(long code)
{
    if (code == 0)                          return K26HTTP_ERR_INTERNAL;
    if (code >= 200 && code < 300)          return K26HTTP_OK;
    if (code == 429)                        return K26HTTP_ERR_429;
    if (code >= 400 && code < 500)          return K26HTTP_ERR_4XX;
    if (code >= 500 && code < 600)          return K26HTTP_ERR_5XX;
    return K26HTTP_ERR_INTERNAL;
}

static K26HTTPStatus map_curl_err(CURLcode cc)
{
    switch (cc) {
    case CURLE_OPERATION_TIMEDOUT: return K26HTTP_ERR_TIMEOUT;
    case CURLE_COULDNT_RESOLVE_HOST:
    case CURLE_COULDNT_RESOLVE_PROXY:
    case CURLE_COULDNT_CONNECT:
    case CURLE_INTERFACE_FAILED:   return K26HTTP_ERR_OFFLINE;
    case CURLE_ABORTED_BY_CALLBACK:
    case CURLE_WRITE_ERROR:        return K26HTTP_ERR_CANCELLED;
    default:                       return K26HTTP_ERR_INTERNAL;
    }
}

K26HTTPStatus k26http_do(const K26HTTPRequest *req, K26HTTPBuf *out, long *http_status)
{
    if (!req || !req->url || !out) return K26HTTP_ERR_INTERNAL;
    if (http_status) *http_status = 0;

    CURL *c = curl_easy_init();
    if (!c) return K26HTTP_ERR_INTERNAL;

    struct curl_slist *hdr = NULL;
    if (req->headers) {
        for (size_t i = 0; req->headers[i]; i++) {
            struct curl_slist *next = curl_slist_append(hdr, req->headers[i]);
            if (!next) {
                curl_slist_free_all(hdr);
                curl_easy_cleanup(c);
                return K26HTTP_ERR_INTERNAL;
            }
            hdr = next;
        }
    }

    WriteCtx wc = {
        .out = out,
        .max_bytes = req->max_bytes,
        .should_cancel = req->should_cancel,
        .user = req->user,
    };

    curl_easy_setopt(c, CURLOPT_URL, req->url);
    curl_easy_setopt(c, CURLOPT_USERAGENT,
                     req->user_agent ? req->user_agent : K26HTTP_UA_DEFAULT);
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_cb);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &wc);
    curl_easy_setopt(c, CURLOPT_TIMEOUT_MS,
                     (long)(req->timeout_ms > 0 ? req->timeout_ms : 30000));
    curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT_MS,
                     (long)(req->connect_timeout_ms > 0 ? req->connect_timeout_ms : 5000));
    curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c, CURLOPT_NOSIGNAL, 1L);
    if (hdr) curl_easy_setopt(c, CURLOPT_HTTPHEADER, hdr);

    if (req->method && strcmp(req->method, "POST") == 0) {
        curl_easy_setopt(c, CURLOPT_POST, 1L);
        curl_easy_setopt(c, CURLOPT_POSTFIELDS, req->body ? req->body : "");
        curl_easy_setopt(c, CURLOPT_POSTFIELDSIZE,
                         (long)(req->body ? strlen(req->body) : 0));
    } else {
        curl_easy_setopt(c, CURLOPT_HTTPGET, 1L);
    }

    CURLcode cc = curl_easy_perform(c);
    long code = 0;
    curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    if (http_status) *http_status = code;

    curl_slist_free_all(hdr);
    curl_easy_cleanup(c);

    if (wc.truncated) return K26HTTP_ERR_TRUNCATED;

    if (cc != CURLE_OK) {
        /* If the callback aborted by returning 0, curl reports WRITE_ERROR
         * or ABORTED_BY_CALLBACK — map to CANCELLED vs INTERNAL by checking
         * whether the user asked for cancellation. */
        if ((cc == CURLE_WRITE_ERROR || cc == CURLE_ABORTED_BY_CALLBACK)
            && req->should_cancel && req->should_cancel(req->user)) {
            return K26HTTP_ERR_CANCELLED;
        }
        return map_curl_err(cc);
    }
    return map_http_status(code);
}

/* -------- thin helpers ----------------------------------------- */

K26HTTPStatus k26http_get(const char *url, const char *const *headers,
                          int timeout_ms, K26HTTPBuf *out, long *http_status)
{
    K26HTTPRequest r = { .method = "GET", .url = url, .headers = headers,
                         .timeout_ms = timeout_ms };
    return k26http_do(&r, out, http_status);
}

K26HTTPStatus k26http_post_json(const char *url, const char *json_body,
                                const char *const *extra_headers,
                                int timeout_ms, K26HTTPBuf *out, long *http_status)
{
    /* Build a headers[] with Content-Type prepended. Cap at 16 user headers. */
    const char *ct = "Content-Type: application/json";
    const char *arr[18] = { ct };
    size_t n = 1;
    if (extra_headers) {
        for (size_t i = 0; extra_headers[i] && n < 17; i++) arr[n++] = extra_headers[i];
    }
    arr[n] = NULL;
    K26HTTPRequest r = { .method = "POST", .url = url, .body = json_body,
                         .headers = arr, .timeout_ms = timeout_ms };
    return k26http_do(&r, out, http_status);
}

char *k26http_header_bearer(const char *token)
{
    if (!token || !token[0]) return NULL;
    size_t n = strlen(token);
    size_t cap = n + sizeof("Authorization: Bearer ");
    char *s = malloc(cap);
    if (!s) return NULL;
    snprintf(s, cap, "Authorization: Bearer %s", token);
    return s;
}

const char *k26http_status_str(K26HTTPStatus s)
{
    switch (s) {
    case K26HTTP_OK:             return "OK";
    case K26HTTP_ERR_OFFLINE:    return "OFFLINE";
    case K26HTTP_ERR_TIMEOUT:    return "TIMEOUT";
    case K26HTTP_ERR_4XX:        return "4XX";
    case K26HTTP_ERR_429:        return "RATE LIMITED";
    case K26HTTP_ERR_5XX:        return "5XX";
    case K26HTTP_ERR_TRUNCATED:  return "TRUNCATED";
    case K26HTTP_ERR_CANCELLED:  return "CANCELLED";
    case K26HTTP_ERR_INTERNAL:   return "INTERNAL";
    }
    return "?";
}
