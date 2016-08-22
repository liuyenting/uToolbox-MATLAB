#ifndef CATMULLROM_INTERP_CUH
#define CATMULLROM_INTERP_CUH

#include <cuda_runtime_api.h>   // cuda*

__device__
float pixelCubicLookup(cudaTextureObject_t, float, float, float);

#endif
