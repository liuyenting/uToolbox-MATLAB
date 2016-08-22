#include "mex.h"

#include "gateway_interp.hpp"

#define GATEWAY_MSGID "util:cppbridge:gateway"

#define CMDSTR prhs[0]
#define HANDLE prhs[1]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]) {
    // Validate the input.
    if (nrhs < 1) {
        mexErrMsgIdAndTxt(GATEWAY_MSGID,
                          "Not enough arguments.");
    } else if (mxGetM(CMDSTR) != 1) {
        mexErrMsgIdAndTxt(GATEWAY_MSGID,
                          "One command string is expected.");
    }

    char *cmd;
    int status;
    // Determine the length of input command.
    size_t cmdLen = mxGetNumberOfElements(CMDSTR) + 1;
    // Allocate the memory for the command.
    cmd = (char *)mxCalloc(cmdLen, sizeof(char));
    if (cmd == NULL) {
        mexErrMsgIdAndTxt(GATEWAY_MSGID,
                          "Not enough heap space to hold the command.");
    }
    // Get the string.
    status = mxGetString(CMDSTR, cmd, cmdLen);
    if (status != 0) {
        mexErrMsgIdAndTxt(GATEWAY_MSGID,
                          "Could not extract the command string.");
    }

    status = interpCmd(cmd, nlhs, plhs, nrhs, prhs);
    if (status != 0) {
        // Unrecognized command.
        mexErrMsgIdAndTxt(GATEWAY_MSGID,
                          "Command not recognized.");
    }
}
