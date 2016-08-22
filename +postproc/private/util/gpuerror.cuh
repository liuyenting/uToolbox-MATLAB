#ifndef GPUERROR_CUH
#define GPUERROR_CUH

#include <cuda_runtime_api.h>

#define CUDA_MSGID "postproc:cuda"

inline void gpuCheckError(const char *errMsg) {
    cudaError_t err = cudaGetLastError() ;
    if (err != cudaSuccess) {
        mexErrMsgIdAndTxt(CUDA_MSGID, errMsg);
        cudaDeviceReset();
        exit(EXIT_FAILURE);
    }
}

#endif
