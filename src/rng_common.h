#ifndef _MMC_RNG_COMMON_H
#define _MMC_RNG_COMMON_H

#define TWO_PI     (M_PI*2.0)

// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen(RandType t[RAND_BUF_LEN]){
    float ran=rand_uniform01(t);
    return ((ran==0.f)?LOG_RNG_MAX:(-logf(ran)));
}
// generate [0,1] random number for the next arimuthal angle
__device__ float rand_next_aangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_next_zangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for reflection test
__device__ float rand_next_reflect(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_do_roulette(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
#endif
