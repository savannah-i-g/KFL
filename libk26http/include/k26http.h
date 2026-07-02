#ifndef K26HTTP_H
#define K26HTTP_H

#include <stdbool.h>
#include <stddef.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * libk26http — provider-neutral HTTP client for K26 programs.
 *
 * Wraps libcurl with:
 *   - A one-request-one-response primitive (k26http_do).
 *   - Thin GET/POST-JSON convenience wrappers.
 *   - An SSE streamer that delivers each `data:` payload to a callback.
 *   - A layered network-status probe that says WHICH layer failed
 *     (iface/DNS/route) so TUIs can give an actionable hint.
 *
 * Explicitly NOT in scope: cJSON wrappers (callers parse), cookie
 * jars, multipart upload, OAuth, WebSockets. Keep the surface small.
 *
 * Thread model: each k26http_do() creates and destroys its own libcurl
 * easy handle — no shared state across calls. The global libcurl init
 * is refcounted via k26http_lib_ref/unref; every process that uses the
 * library should call ref() once at startup and unref() at shutdown.
 * ================================================================== */

typedef enum {
    K26HTTP_OK = 0,
    K26HTTP_ERR_OFFLINE,     /* no iface / DNS / connect — caller should fall
                                back to cache and surface an OFFLINE state */
    K26HTTP_ERR_TIMEOUT,
    K26HTTP_ERR_4XX,         /* generic 4xx — caller inspects http_status   */
    K26HTTP_ERR_429,         /* explicit — rate-limit governor must back off */
    K26HTTP_ERR_5XX,
    K26HTTP_ERR_TRUNCATED,   /* response body exceeded max_bytes            */
    K26HTTP_ERR_CANCELLED,   /* should_cancel returned true                 */
    K26HTTP_ERR_INTERNAL,
} K26HTTPStatus;

/* Growable buffer used as the response sink. Caller initialises with
 * k26http_buf_init, reads via data/len, frees with k26http_buf_free. */
typedef struct {
    char  *data;
    size_t len;
    size_t cap;
} K26HTTPBuf;

void k26http_buf_init(K26HTTPBuf *);
void k26http_buf_free(K26HTTPBuf *);
bool k26http_buf_append(K26HTTPBuf *, const void *src, size_t n);

typedef struct {
    const char *method;               /* "GET"|"POST" — NULL defaults to GET */
    const char *url;                  /* required                           */
    const char *body;                 /* request body (NULL for GET)        */
    const char *const *headers;       /* NULL-terminated full-line headers,
                                         e.g. "Authorization: Bearer <tok>" */
    int         timeout_ms;           /* 0 => 30000                         */
    int         connect_timeout_ms;   /* 0 => 5000                          */
    size_t      max_bytes;            /* 0 => 16 MiB; body > cap => TRUNCATED */
    const char *user_agent;           /* NULL => "libk26http/1.0" */

    /* Optional progress cancel-check. Returning true aborts mid-transfer
     * and k26http_do returns K26HTTP_ERR_CANCELLED. */
    bool (*should_cancel)(void *user);
    void *user;
} K26HTTPRequest;

K26HTTPStatus k26http_do      (const K26HTTPRequest *, K26HTTPBuf *out, long *http_status);
K26HTTPStatus k26http_get     (const char *url, const char *const *headers,
                               int timeout_ms, K26HTTPBuf *out, long *http_status);
K26HTTPStatus k26http_post_json(const char *url, const char *json_body,
                                const char *const *extra_headers,
                                int timeout_ms, K26HTTPBuf *out, long *http_status);

/* Build an "Authorization: Bearer <token>" header line. Returns malloc'd
 * string; caller frees. NULL on OOM or empty token. */
char *k26http_header_bearer(const char *token);

/* ---- Net probe -------------------------------------------------- */

typedef enum {
    K26HTTP_NET_ONLINE = 0,
    K26HTTP_NET_NO_IFACE,
    K26HTTP_NET_NO_DNS,
    K26HTTP_NET_NO_ROUTE,
} K26HTTPNetState;

/* Non-blocking probe: does this host+port accept a TCP connect within
 * timeout_ms? host/port may be NULL to use the library defaults
 * ("1.1.1.1" / "443"). timeout_ms==0 => 2000 ms default. */
K26HTTPNetState k26http_net_probe(const char *host, const char *port, int timeout_ms);
const char     *k26http_net_state_str(K26HTTPNetState);

/* ---- SSE streaming ---------------------------------------------- */

/* Called once per `data:` payload received from the server. The string
 * is the payload with SSE framing stripped (no "data:" prefix, no
 * trailing newline). Return false to abort the stream; the request
 * then resolves with K26HTTP_ERR_CANCELLED. */
typedef bool (*K26HTTPSSECb)(const char *event_data, void *user);

/* Perform a POST/GET whose response is SSE framed (server:
 *   `Content-Type: text/event-stream`). Emits event_data callbacks
 * frame-by-frame; the response body is NOT buffered into out. */
K26HTTPStatus k26http_sse(const K26HTTPRequest *, K26HTTPSSECb, void *user);

/* ---- Lifecycle + diagnostics ------------------------------------ */

/* Refcounted wrappers around curl_global_init / _cleanup. Safe to call
 * once per-process; nested calls are cheap. */
void        k26http_lib_ref  (void);
void        k26http_lib_unref(void);

const char *k26http_status_str(K26HTTPStatus);

#ifdef __cplusplus
}
#endif

#endif /* K26HTTP_H */
