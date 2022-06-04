#ifndef _LATTICE_GENERATOR_H_
#define _LATTICE_GENERATOR_H_

#include <gsl/gsl_rng.h>

/********************* lattice_generator() *********************/
/** generates a L-size list of spins                          **/
/** with values -1 and 1 and probability 0.5 for each         **/
/**                                                           **/
/** PARAMETERS                                                **/
/**        L: number of spins                                 **/
/**     Spin: list of spin values                             **/
/**   gslran: list of random numbers                          **/
/**                                                           **/
/** RETURNS:                                                  **/
/**        0 -> OK                                            **/
/***************************************************************/


void lattice_generator(int L, int Spin[L], gsl_rng *gslran){
  for (int i = 0; i < L; ++i) {
    double s = gsl_rng_uniform(gslran);
    if (s < 0.5) {
      Spin[i] = -1;
    } else {
      Spin[i] = 1;
    }
  }
}

#endif /* _LATTICE_GENERATOR_H_ */