#include "k26compute.h"

const char *k26compute_status_str(K26CStatus s)
{
    switch (s) {
    case K26C_OK:           return "ok";
    case K26C_ERR_OOM:      return "out of memory";
    case K26C_ERR_INVAL:    return "invalid argument or dimension mismatch";
    case K26C_ERR_SINGULAR: return "matrix is singular";
    case K26C_ERR_CONVERGE: return "iterative method failed to converge";
    case K26C_ERR_RANGE:    return "result out of representable range";
    case K26C_ERR_INTERNAL: return "internal error";
    }
    return "unknown status";
}
