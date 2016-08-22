#include "mex.h"

#include <cstring>          // strcmp

#include "shear.cuh"        // Shear

#include "class_handle.hpp"

#define INTERFACE_MSGID "postproc:shear:interface"

#define CMDSTR prhs[0]
#define HANDLE prhs[1]
#define IMGOUT plhs[0]

#define isCommand(str) !strcmp(str, cmd)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]) {
    // Validate the input.
    if (nrhs < 1) {
        mexErrMsgIdAndTxt(INTERFACE_MSGID,
                          "Not enough arguments.");
    } else if (mxGetM(CMDSTR) != 1) {
        mexErrMsgIdAndTxt(INTERFACE_MSGID,
                          "One command string is expected.");
    }

    char *cmd;
    int status;
    // Determine the length of input command.
    size_t cmdLen = mxGetN(CMDSTR) + 1;
    // Allocate the memory for instruction.
    cmd = (char *)mxCalloc(cmdLen, sizeof(char));
    if (cmd == NULL) {
        mexErrMsgIdAndTxt(INTERFACE_MSGID,
                          "Not enough heap space to hold the command.");
    }
    // Get the string.
    status = mxGetString(CMDSTR, cmd, cmdLen);
    if (status != 0) {
        mexErrMsgIdAndTxt(INTERFACE_MSGID,
                          "Could not convert the command string.");
    }

    /*
     * New
     */
    if (isCommand("new")) {
        if (nlhs != 1) {
            mexErrMsgIdAndTxt(INTERFACE_MSGID,
                              "One output expected.");
        }
        // Retrieve the class instance pointer.
        plhs[0] = convertPtr2Mat<Shear>(new Shear);
        return;
    }

    // If we are not creating an object, then the 2nd agrument should be the
    // object handle. Therefore, at least 2 arguments should at RHS.
    if (nrhs < 2) {
        mexErrMsgIdAndTxt(INTERFACE_MSGID,
                          "At least two inputs should be provided.");
    }

    /*
     * Delete
     */
    if (isCommand("delete")) {
        // Destroy the object.
        destroyObject<Shear>(HANDLE);
        // Ignore rest of the inputs.
        if ((nlhs > 0) || (nrhs > 2)) {
            mexWarnMsgIdAndTxt(INTERFACE_MSGID,
                               "Unexpected arguments ignored.");
        }
        return;
    }

    // Retrieve the class instance since we are neither creating nor deleting
    // it.
    Shear *objInstance = convertMat2Ptr<Shear>(HANDLE);

    /*
     * Call member functions.
     */
    if (isCommand("setacqparam")) {
        if (nrhs < 3) {
            mexWarnMsgIdAndTxt(INTERFACE_MSGID,
                               "Parameter not provided.");
        }
        objInstance->setAcqParam(prhs[2]);
    } else if (isCommand("loadstack")) {
        if (nrhs < 3) {
            mexWarnMsgIdAndTxt(INTERFACE_MSGID,
                               "Image stack not provided.");
        }
        objInstance->loadStack(prhs[2]);
    } else if (isCommand("execute")) {
        objInstance->execute();
    } else if (isCommand("retrieveresult")) {
        if (nlhs != 1) {
            mexErrMsgIdAndTxt(INTERFACE_MSGID,
                              "One output expected.");
        }
        objInstance->retrieveResult(&IMGOUT);
    } else {
        // Unrecognized command.
        mexErrMsgIdAndTxt(INTERFACE_MSGID,
                          "Command not recognized.");
    }
}
