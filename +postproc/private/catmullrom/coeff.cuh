#ifndef CATMULLROM_COEFF_CUH
#define CATMULLROM_COEFF_CUH

__host__ __device__
float w0(float a) {
    // -0.5f*a + a*a - 0.5f*a*a*a
    return a*(-0.5f + a*(1.0f - 0.5f*a));
}

__host__ __device__
float w1(float a) {
    // 1.0f - 2.5f*a*a + 1.5f*a*a*a
    return 1.0f + a*a*(-2.5f + 1.5f*a);
}

__host__ __device__
float w2(float a) {
    // 0.5f*a + 2.0f*a*a - 1.5f*a*a*a
    return a*(0.5f + a*(2.0f - 1.5f*a));
}

__host__ __device__
float w3(float a) {
    // -0.5f*a*a + 0.5f*a*a*a
    return a*a*(-0.5f + 0.5f*a);
}

#endif
