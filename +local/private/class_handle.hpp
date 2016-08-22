#ifndef CLASS_HANDLE_HPP
#define CLASS_HANDLE_HPP

#include "mex.h"

#include <cstdint>      // uint64_t
#include <string>       // std::string
#include <cstring>      // strcmp
#include <typeinfo>     // typeid

#define BASECLASS_MSGID "postproc:class_handle"

#define CLASS_HANDLE_SIGNATURE 0xFF00F0A5

template <typename Base>
class ClassHandle {
private:
    uint32_t signature;
    std::string name;
    Base *classPtr;

public:
    ClassHandle(Base *_ptr)
        : classPtr(_ptr), name(typeid(Base).name()) {
        signature = CLASS_HANDLE_SIGNATURE;
    }

    bool isValid() {
        return ((signature == CLASS_HANDLE_SIGNATURE) &&
                !strcmp(name.c_str(), typeid(Base).name()));
    }

    Base * ptr() {
        return classPtr;
    }

    ~ClassHandle() {
        signature = 0;
        delete classPtr;
    }
};

template <typename Base>
inline mxArray * convertPtr2Mat(Base *ptr) {
    // To ensure the resource is remained in memory.
    mexLock();
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(out)) =
        reinterpret_cast<uint64_t>(new ClassHandle<Base>(ptr));
    return out;
}

template <typename Base>
inline ClassHandle<Base> * convertMat2HandlePtr(const mxArray *in)
{
    // Verify input type is a 64-bit pointer adderss.
    if (mxGetNumberOfElements(in) != 1 ||
        mxGetClassID(in) != mxUINT64_CLASS ||
        mxIsComplex(in)) {
        mexErrMsgIdAndTxt(BASECLASS_MSGID,
                          "Input must be a real uint64 scalar.");
    }
    ClassHandle<Base> *ptr =
        reinterpret_cast<ClassHandle<Base> *>(*((uint64_t *)mxGetData(in)));
    if (!ptr->isValid()) {
        mexErrMsgIdAndTxt(BASECLASS_MSGID,
                          "Handle not valid.");
    }
    return ptr;
}

template <typename Base>
inline Base * convertMat2Ptr(const mxArray *in) {
    return convertMat2HandlePtr<Base>(in)->ptr();
}

template <typename Base>
inline void destroyObject(const mxArray *in) {
    delete convertMat2HandlePtr<Base>(in);
    mexUnlock();
}

#endif
