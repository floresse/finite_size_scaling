#ifndef _ENERGY_1D_LRI_H_EXT_H_
#define _ENERGY_1D_LRI_H_EXT_H_

/********************* energy_1d_lri_hext() ********************/
/** computes the energy for 1D long-range interaction systems **/
/** with external magnetic field                              **/
/**                                                           **/
/** PARAMETERS                                                **/
/**        L: ith spin                                        **/
/**     Spin: relations of bonds                              **/
/**        M: system magnetization                            **/
/**        h: external magnetic field                         **/ 
/**     dist: distance between spins ith and jth              **/  
/**                                                           **/
/** RETURNS:                                                  **/
/**        System Energy                                      **/
/***************************************************************/


double energy_1d_lri_h_ext(int L, int Spin[L],int M,double h,double dist[L])
{
  int i,j;
  double E = 0;
  for (i = 0; i < L; ++i) for (j = i+1; j < L; ++j ) {
    E += -Spin[i]*Spin[j]*dist[j-i-1];
  } 
  E += -h*M;
  return E;
}

#endif /* _ENERGY_1D_LRI_H_EXT_H_ */