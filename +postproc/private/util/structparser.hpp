#ifndef STRUCTPARSER_HPP
#define STRUCTPARSER_HPP

#include "mex.h"

#include <cassert>

float parseFloatField(const mxArray *in, const char *fname) {
    const mxArray *fptr;
    float *dptr;

    fptr = mxGetField(in, 0, fname);
    assert(fptr != NULL);
    dptr = (float *)mxGetData(fptr);

    // Dereferenced.
    return *dptr;
}

#endif
