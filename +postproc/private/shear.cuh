#ifndef SHEAR_CUH
#define SHEAR_CUH

#include "mex.h"

#include <cmath>                // std::cosf
#include <cstdint>              // uint16_t
#include <cuda_runtime_api.h>   // cuda*

// Directly use the constant without defining _USE_MATH_DEFINES first.
#define M_PI 3.14159265358979323846f

#define CORE_MSGID "postproc:shear:core"

#define deg2rad(d) (d * M_PI/180.0f)

struct AcqParams {
    float objAngle, zStepWidth, pixelWidth;

    AcqParams() {
        objAngle = zStepWidth = pixelWidth = 0.0f;
    }

    float shearFactor() {
        float cosVal = std::cosf( deg2rad(objAngle) );
        return zStepWidth * cosVal / pixelWidth;
    }

    float isReversed() {
        return (objAngle < 0);
    }
};

class Shear {
private:
    AcqParams acqParams;

    cudaMemcpy3DParms cpParam;

    size_t oldSize[2];
    size_t newSize[2];
    size_t nLayers;

    cudaArray *devOldImg;
    float *devNewImg;

    cudaTextureObject_t texObj;

    const size_t blkSize = 16;

public:
    Shear();

    void setAcqParam(const mxArray *);
    void loadStack(const mxArray *);
    void execute();
    void retrieveResult(mxArray **);

    ~Shear();

private:
    void preallocateWorkspace(const size_t *);
    void generateWorkspace();
    void releaseWorkspace();

    // Pre-calculate the workspace size.
    void saveOldSize(const size_t *);
    void estimateNewSize();

    // Device access.
    void uploadWorkspaceSizeToDevice();
    void setupImageCopyParameter();
    void copyImageToDevice(const uint16_t *);
    void bindImageToTexture();
};

#endif
