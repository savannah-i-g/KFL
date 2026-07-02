#include "k26geo.h"

const char *k26geo_status_str(K26GeoStatus s)
{
    switch (s) {
    case K26GEO_OK:              return "ok";
    case K26GEO_ERR_INVAL:       return "invalid argument";
    case K26GEO_ERR_OOM:         return "out of memory";
    case K26GEO_ERR_IO:          return "i/o error";
    case K26GEO_ERR_PARSE:       return "parse error";
    case K26GEO_ERR_RANGE:       return "value out of range";
    case K26GEO_ERR_UNSUPPORTED: return "unsupported";
    case K26GEO_ERR_INTERNAL:    return "internal error";
    }
    return "unknown";
}
