/* v4 — vec4 operations. v2/v3 hot paths live inline in the public
 * header; v4 ops are out-of-line because they appear less frequently
 * in the renderer's inner loops and because keeping the header tight
 * helps compile time for kflc-emitted TUs. */
#include "k26m3d.h"

K26V4 k26m3d_v4_add(K26V4 a, K26V4 b)
{
    return k26m3d_v4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

K26V4 k26m3d_v4_sub(K26V4 a, K26V4 b)
{
    return k26m3d_v4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

K26V4 k26m3d_v4_scale(K26V4 a, double s)
{
    return k26m3d_v4(a.x * s, a.y * s, a.z * s, a.w * s);
}

double k26m3d_v4_dot(K26V4 a, K26V4 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
