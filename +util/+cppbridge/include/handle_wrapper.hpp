#ifndef HANDLE_WRAPPER_HPP
#define HANDLE_WRAPPER_HPP

#include "mex.h"

#include <cstdint>      // uint64_t
#include <string>       // std::string
#include <cstring>      // strcmp
#include <typeinfo>     // typeid

#define HANDLEWRAPPER_MSGID "util:cppbridge:handlewrapper"

#define CLASS_SIGNATURE 0xFF00F0A5

template <typename Base>
class HandleWrapper {
private:
    uint32_t signature;
    std::string name;
    Base *classPtr;

public:
    HandleWrapper(Base *_ptr)
        : classPtr(_ptr), name(typeid(Base).name()) {
        signature = CLASS_SIGNATURE;
    }

    bool isValid() {
        return ((signature == CLASS_SIGNATURE) &&
                !strcmp(name.c_str(), typeid(Base).name()));
    }

    Base * ptr() {
        return classPtr;
    }

    ~HandleWrapper() {
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
        reinterpret_cast<uint64_t>(new HandleWrapper<Base>(ptr));
    return out;
}

template <typename Base>
inline HandleWrapper<Base> * convertMat2HandlePtr(const mxArray *in)
{
    // Verify input type is a 64-bit pointer adderss.
    if (mxGetNumberOfElements(in) != 1 ||
        mxGetClassID(in) != mxUINT64_CLASS ||
        mxIsComplex(in)) {
        mexErrMsgIdAndTxt(HANDLEWRAPPER_MSGID,
                          "Input must be a real uint64 scalar.");
    }
    HandleWrapper<Base> *ptr =
        reinterpret_cast<HandleWrapper<Base> *>(*((uint64_t *)mxGetData(in)));
    if (!ptr->isValid()) {
        mexErrMsgIdAndTxt(HANDLEWRAPPER_MSGID,
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
