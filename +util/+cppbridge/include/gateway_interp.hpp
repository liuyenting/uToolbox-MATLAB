#ifndef GATEWAY_INTERP_HPP
#define GATEWAY_INTERP_HPP

// Interpret the command string. Return 0 if success.
int interpCmd(const char *cmd,
              int nlhs, mxArray *plhs[],
              int nrhs, mxArray *prhs[]);

#endif
