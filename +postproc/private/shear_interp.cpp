#include "mex.h"

#include <cstring>              // strcmp

#include "gateway_interp.hpp"   // interpCmd
#include "handle_wrapper.hpp"   // convertPtr2Mat, convertMat2Ptr
#include "shear.cuh"            // Shear

#define INTERP_MSGID "postproc:shear:interp"

#define CMDSTR prhs[0]
#define HANDLE prhs[1]
#define IMGOUT plhs[0]

#define isCommand(str) !strcmp(str, cmd)

int interpCmd(const char *cmd,
              int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]) {
    if (isCommand("new")) {
        if (nlhs != 1) {
            mexErrMsgIdAndTxt(INTERP_MSGID,
                              "One output expected.");
        }
        // Retrieve the class instance pointer.
        plhs[0] = convertPtr2Mat<Shear>(new Shear);
        return 0;
    }

    // If we are not creating an object, then the 2nd agrument should be the
    // object handle. Therefore, at least 2 arguments should at RHS.
    if (nrhs < 2) {
        mexErrMsgIdAndTxt(INTERP_MSGID,
                          "At least two inputs should be provided.");
    }

    if (isCommand("delete")) {
        // Destroy the object.
        destroyObject<Shear>(HANDLE);
        // Ignore rest of the inputs.
        if ((nlhs > 0) || (nrhs > 2)) {
            mexWarnMsgIdAndTxt(INTERP_MSGID,
                               "Unexpected arguments ignored.");
        }
        return 0;
    }

    // Retrieve the class instance since we are neither creating nor deleting
    // it.
    Shear *objInstance = convertMat2Ptr<Shear>(HANDLE);

    /*
    * Call member functions.
    */
    if (isCommand("setacqparam")) {
        if (nrhs < 3) {
            mexWarnMsgIdAndTxt(INTERP_MSGID,
                               "Parameter not provided.");
        }
        objInstance->setAcqParam(prhs[2]);
    } else if (isCommand("loadstack")) {
        if (nrhs < 3) {
            mexWarnMsgIdAndTxt(INTERP_MSGID,
                               "Image stack not provided.");
        }
        objInstance->loadStack(prhs[2]);
    } else if (isCommand("execute")) {
        objInstance->execute();
    } else if (isCommand("retrieveresult")) {
        if (nlhs != 1) {
            mexErrMsgIdAndTxt(INTERP_MSGID,
                              "One output expected.");
        }
        objInstance->retrieveResult(&IMGOUT);
    } else {
        return -1;
    }
    return 0;
}
