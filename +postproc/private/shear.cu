#include "mex.h"
#include "gpu/mxGPUArray.h"

#include <cstdint>              // uint16_t, uint8_t
#include <cassert>              // assert
#include <algorithm>            // std::fill
#include <cuda_runtime_api.h>   // cuda*

#include "util/structparser.hpp"    // parseFloatField
#include "util/gpuerror.cuh"        // gpuCheckError
#include "kernels/shear.cuh"        // shearLayer

#include "shear.cuh"

Shear::Shear() {
    // Reset the acquisition paramters.
    acqParams = {};

    // Reset copy related parameters.
    cpParam = {};

    // Reset the layer counter only, since this will trigger the copy.
    nLayers = 0;

    // Initialize MathWorks GPU API.
    mxInitGPU();
}

void Shear::setAcqParam(const mxArray *in) {
    acqParams.objAngle = parseFloatField(in, "ObjectiveAngle");
    acqParams.zStepWidth = parseFloatField(in, "ZStepWidth");
    acqParams.pixelWidth = parseFloatField(in, "PixelWidth");

    mexPrintf(" ** Acquisition paramter parsed...\n    - Angle: %f[deg]\n    - Z Width: %f[rad]\n    - Px Size: %f[rad]\n",
              acqParams.objAngle, acqParams.zStepWidth, acqParams.pixelWidth);

    // Force reset in order to trigger reallocate the workspace.
    nLayers = 0;
    // TODO: Free the WS if cudaArray is not NULL.
}

void Shear::loadStack(const mxArray *in) {
    // Turn mxArray to mxGPUArray.
    const mxGPUArray *inArr = mxGPUCreateFromMxArray(in);

    if (nLayers == 0) {
        const size_t *inSize = (const size_t *)mxGPUGetDimensions(inArr);
        preallocateWorkspace(inSize);
        setupImageCopyParameter();
    }

    const uint16_t *inImg = (const uint16_t *)mxGPUGetDataReadOnly(inArr);
    copyImageToDevice(inImg);

    mxGPUDestroyGPUArray(inArr);
}

void Shear::preallocateWorkspace(const size_t *_oriSize) {
    // Pre-calculate the image sizes along the pipeline.
    saveOldSize(_oriSize);
    estimateNewSize();

    mexPrintf(" ** Input size [%d, %d, %d] ... Output size [%d, %d, %d]\n",
              oldSize[0], oldSize[1], nLayers,
              newSize[0], newSize[1], nLayers);

    generateWorkspace();
}

void Shear::saveOldSize(const size_t *size) {
    std::copy(size, size+2, oldSize);
    nLayers = size[2];
}

void Shear::estimateNewSize() {
    // Alias for the original dimension.
    size_t nx = oldSize[0];
    size_t ny = oldSize[1];

    // Calculate amount of shifted pixels.
    size_t diff = static_cast<size_t>(acqParams.shearFactor()*(ny-1));

    // Fill the values.
    newSize[0] = nx + diff;
    newSize[1] = ny;
}

/*
void Shear::estimateRotatedSize() {
    // Check whether we have to bypass rotating.
    if (acqParams.noRotate) {
        std::copy(shSize, shSize+2, rotSize);
        return;
    }

    // Alias for the sheared size.
    size_t nx = shSize[0];
    size_t ny = shSize[1];

    // Unit length.
    float unitSize[4][2] = { {1, 1}, {-1, 1}, {-1, -1}, {1, -1} };
    // Generate reverse matrix.
    fillRotationMatrix(-acqParams.stackRotatedAngle());
    // Apply rotation to all the end points.
    for (uint8_t i = 0; i < 4; i++) {
        float x = unitSize[i][0] - 1.0f;
        float y = unitSize[i][1] - 1.0f;

        float u = rotMat[0]*x + rotMat[1]*y;
        float v = rotMat[2]*x + rotMat[3]*y;

        unitSize[i][0] = u;
        unitSize[i][1] = v;
    }
    // Find the minimum and maximum value of each axis.
    float uMin, uMax, vMin, vMax;
    uMin = uMax = unitSize[0][0];
    vMin = vMax = unitSize[0][1];
    for (uint8_t i = 1; i < 4; i++) {
        if (unitSize[i][0] < uMin) {
            uMin = unitSize[i][0];
        } else if (unitSize[i][0] > uMax) {
            uMax = unitSize[i][0];
        }

        if (unitSize[i][1] < vMin) {
            vMin = unitSize[i][1];
        } else if (unitSize[i][1] > vMax) {
            vMax = unitSize[i][1];
        }
    }
    // Calculate the length.
    float nu = ceil(uMax - uMin);
    assert(nu > 0 );
    float nv = ceil(vMax - vMin);
    assert(nv > 0);

    // Fill back to the array.
    rotSize[0] = (size_t)(nu * nx/2.0f);
    rotSize[1] = (size_t)(nv * ny/2.0f);

    // TODO: Calcualte the cropping factor.

    mexWarnMsgIdAndTxt(CORE_MSGID,
                       " >> Rotated size: [%d, %d, %d]",
                       rotSize[0], rotSize[1], nLayers);
}
*/

void Shear::generateWorkspace() {
    cudaMemcpyToSymbol(devNewSize, newSize, sizeof(size_t)*2,
                       0, cudaMemcpyHostToDevice);
    gpuCheckError("Failed to copy workspace sizes to device constant memory.");

    const float shFact = acqParams.shearFactor();
    cudaMemcpyToSymbol(devShFact, &shFact, sizeof(float),
                       0, cudaMemcpyHostToDevice);

    // NOTE: cudaChannelFormatDesc cannot reuse, otherwise it will cause memory
    // deallocation error.
    cudaChannelFormatDesc oldChDesc =
        cudaCreateChannelDesc(16, 0, 0, 0, cudaChannelFormatKindUnsigned);
    cudaMalloc3DArray(&devOldImg,
                      &oldChDesc,
                      make_cudaExtent(oldSize[0], oldSize[1], nLayers),
                      cudaArrayLayered);
    gpuCheckError("Failed to allocate old image buffere on device.");

    const size_t nElem = newSize[0]*newSize[1]*nLayers;
    cudaMalloc((void **)&devNewImg, sizeof(float)*nElem);
    gpuCheckError("Failed to allocate new image buffer on device.");
}

/*
void Shear::fillRotationMatrix(const float angle) {
    float sinVal = std::sinf(angle);
    float cosVal = std::cosf(angle);

    rotMat[0] = cosVal;
    rotMat[1] = -sinVal;
    rotMat[2] = sinVal;
    rotMat[3] = cosVal;
}
*/

void Shear::setupImageCopyParameter() {
    // Setup non-changing parameters.
    cpParam.srcPos = make_cudaPos(0, 0, 0);
    cpParam.dstPos = make_cudaPos(0, 0, 0);
    cpParam.dstArray = devOldImg;
    cpParam.extent = make_cudaExtent(oldSize[0], oldSize[1], nLayers);
    cpParam.kind = cudaMemcpyDeviceToDevice;
}

void Shear::copyImageToDevice(const uint16_t *in) {
    // Alias for the original dimension.
    const size_t nx = oldSize[0];
    const size_t ny = oldSize[1];
    // Complete the copy parameters.
    cpParam.srcPtr =
        make_cudaPitchedPtr((void *)in, sizeof(uint16_t)*nx, nx, ny);
    // Start the copy.
    cudaMemcpy3D(&cpParam);
    gpuCheckError("Failed to copy original image.");
}

void Shear::execute() {
    bindImageToTexture();

    const size_t nu = newSize[0];
    const size_t nv = newSize[1];
    const size_t nw = nLayers;
    // Setup the threads, grid size is ceiled.
    const dim3 blockSize(blkSize, blkSize, 1);
    const dim3 gridSize((nu+blockSize.x-1)/blockSize.x,
                        (nv+blockSize.y-1)/blockSize.y);

    // Iterate through all the layers.
    bool reversed = acqParams.isReversed();
    for (size_t iw = 0; iw < nw; iw++) {
        shearLayer<<<gridSize, blockSize>>>(devNewImg, texObj, iw, reversed);
    }
}

void Shear::bindImageToTexture() {
    // Specify texture resource.
    cudaResourceDesc resDesc = {};
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = devOldImg;

    // Specify texture paramters, voxels are accessed through raw coordinates.
    cudaTextureDesc texDesc = {};
    texDesc.addressMode[0] = cudaAddressModeBorder;
    texDesc.addressMode[1] = cudaAddressModeBorder;
    // Nearest neighbor: Point, Interpolation: Linear.
    texDesc.filterMode = cudaFilterModePoint;
    // Output as normalized float [0, 1] per pixel.
    texDesc.readMode = cudaReadModeNormalizedFloat;
    // Access with original coordinates.
    texDesc.normalizedCoords = false;

    // Bind the array to texture.
    cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);
    gpuCheckError("Failed to bind device data to texture.");
}

void Shear::retrieveResult(mxArray **out) {
    const size_t outSize[3] = { newSize[0], newSize[1], nLayers };
    mxGPUArray *outArr =
        mxGPUCreateGPUArray(3, outSize, mxSINGLE_CLASS, mxREAL,
                            MX_GPU_DO_NOT_INITIALIZE);
    float *outImg = (float *)mxGPUGetData(outArr);
    // NOTE: result is assumed to match the size of outImg.
    const size_t nElem = outSize[0]*outSize[1]*outSize[2];
    cudaMemcpy(outImg, devNewImg, sizeof(float)*nElem,
               cudaMemcpyDeviceToDevice);
    gpuCheckError("Failed to move result from device to host.");

    // Output the result to MATLAB.
    *out = mxGPUCreateMxArrayOnGPU(outArr);
    mxGPUDestroyGPUArray(outArr);
}

Shear::~Shear() {
    releaseWorkspace();
    mexPrintf(" ** GPU resources released\n");
}

void Shear::releaseWorkspace() {
    cudaFreeArray(devOldImg);
    cudaFree(devNewImg);
    gpuCheckError("Failed to free the resources on device.");
}
