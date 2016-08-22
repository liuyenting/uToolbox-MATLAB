#ifndef KERNEL_SHEAR_CUH
#define KERNEL_SHEAR_CUH

#include <cuda_runtime_api.h>       // cuda*

__constant__ size_t devNewSize[2];
__constant__ float devShFact;

__global__
void shearLayer(float *out,
                cudaTextureObject_t texObj, const size_t layer,
                const bool isReversed) {
    // Calculate the worker location.
    const size_t u = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t v = blockIdx.y*blockDim.y + threadIdx.y;
    // Output size.
    const size_t nu = devNewSize[0];
    const size_t nv = devNewSize[1];
    // Make sure we are inside the bounding box of new image.
    if ((u >= nu) || (v >= nv)) {
        return;
    }
    // Calculate the linear memory location.
    const size_t i = (layer * (nu*nv)) + (v * nu) + u;

    // Calculate the shifted location.
    const float offset = (isReversed) ? ((nv-1)-v) : v;
    const float x = u - devShFact*offset;
    const float y = v;

    out[i] = tex2DLayered<float>(texObj, x, y, layer);
}

#endif
