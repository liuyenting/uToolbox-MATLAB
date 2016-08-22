#include <cuda_runtime_api.h>   // cuda*

#include "coeff.cuh"            // w0, w1, w2, w3

#include "interp.cuh"

__device__
float cubicInterp(float x,
                  float c0, float c1, float c2, float c3) {
    float r;
    r = c0 * w0(x);
    r += c1 * w1(x);
    r += c2 * w2(x);
    r += c3 * w3(x);
    return r;
}

__device__
float pixelCubicLookup(cudaTextureObject_t texObj,
                       float x, float y, float layer) {
    #define texlookup(x, y) tex2DLayered<float>(texObj, x, y, layer)

    // Integer pixel location.
    float ix = floor(x);
    float iy = floor(y);
    // Fraction for interpolation.
    float fx = x - ix;
    float fy = y - iy;

    // Perform two cubic interpolations.
    return cubicInterp(fy,
                       cubicInterp(fx, texlookup(ix-1, iy-1), texlookup(ix, iy-1), texlookup(ix+1, iy-1), texlookup(ix+2, iy-1)),
                       cubicInterp(fx, texlookup(ix  , iy  ), texlookup(ix, iy  ), texlookup(ix+1, iy  ), texlookup(ix+2, iy  )),
                       cubicInterp(fx, texlookup(ix-1, iy+1), texlookup(ix, iy+1), texlookup(ix+1, iy+1), texlookup(ix+2, iy+1)),
                       cubicInterp(fx, texlookup(ix-1, iy+2), texlookup(ix, iy+2), texlookup(ix+1, iy+2), texlookup(ix+2, iy+2)));
}
