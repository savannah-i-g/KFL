#include "k26http.h"

#include <errno.h>
#include <fcntl.h>
#include <ifaddrs.h>
#include <net/if.h>
#include <netdb.h>
#include <netinet/in.h>
#include <poll.h>
#include <stdbool.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#define K26HTTP_NET_DEFAULT_TIMEOUT_MS 2000
#define K26HTTP_NET_DEFAULT_HOST       "1.1.1.1"
#define K26HTTP_NET_DEFAULT_PORT       "443"

static bool has_iface_addr(void)
{
    struct ifaddrs *list = NULL;
    if (getifaddrs(&list) != 0) return false;
    bool found = false;
    for (struct ifaddrs *ia = list; ia; ia = ia->ifa_next) {
        if (!ia->ifa_addr) continue;
        if (!(ia->ifa_flags & IFF_UP)) continue;
        if (ia->ifa_flags & IFF_LOOPBACK) continue;
        int fam = ia->ifa_addr->sa_family;
        if (fam == AF_INET || fam == AF_INET6) { found = true; break; }
    }
    freeifaddrs(list);
    return found;
}

static long now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long)ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
}

static int remaining_ms(long deadline, long now)
{
    long r = deadline - now;
    if (r < 0) return 0;
    if (r > 60000) return 60000;
    return (int)r;
}

K26HTTPNetState k26http_net_probe(const char *host, const char *port, int timeout_ms)
{
    if (!host || !*host) host = K26HTTP_NET_DEFAULT_HOST;
    if (!port || !*port) port = K26HTTP_NET_DEFAULT_PORT;
    if (timeout_ms <= 0) timeout_ms = K26HTTP_NET_DEFAULT_TIMEOUT_MS;

    long start    = now_ms();
    long deadline = start + timeout_ms;

    if (!has_iface_addr()) return K26HTTP_NET_NO_IFACE;

    struct addrinfo hints = { .ai_socktype = SOCK_STREAM };
    struct addrinfo *res = NULL;
    int gai = getaddrinfo(host, port, &hints, &res);
    if (gai != 0 || !res) {
        if (res) freeaddrinfo(res);
        return K26HTTP_NET_NO_DNS;
    }

    K26HTTPNetState out = K26HTTP_NET_NO_ROUTE;
    for (struct addrinfo *ai = res; ai; ai = ai->ai_next) {
        int fd = socket(ai->ai_family, ai->ai_socktype | SOCK_NONBLOCK, ai->ai_protocol);
        if (fd < 0) continue;

        int rc = connect(fd, ai->ai_addr, ai->ai_addrlen);
        if (rc == 0) { close(fd); out = K26HTTP_NET_ONLINE; break; }
        if (errno != EINPROGRESS) { close(fd); continue; }

        struct pollfd pfd = { .fd = fd, .events = POLLOUT };
        int wait = remaining_ms(deadline, now_ms());
        int pr   = poll(&pfd, 1, wait);
        if (pr > 0) {
            int so_err = 0;
            socklen_t slen = sizeof so_err;
            if (getsockopt(fd, SOL_SOCKET, SO_ERROR, &so_err, &slen) == 0 && so_err == 0) {
                close(fd);
                out = K26HTTP_NET_ONLINE;
                break;
            }
        }
        close(fd);
    }
    freeaddrinfo(res);
    return out;
}

const char *k26http_net_state_str(K26HTTPNetState s)
{
    switch (s) {
    case K26HTTP_NET_ONLINE:   return "ONLINE";
    case K26HTTP_NET_NO_IFACE: return "NO IFACE";
    case K26HTTP_NET_NO_DNS:   return "NO DNS";
    case K26HTTP_NET_NO_ROUTE: return "NO ROUTE";
    }
    return "?";
}
