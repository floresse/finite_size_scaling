#ifndef _FIND_ROOT_H_
#define _FIND_ROOT_H_

/********************* find_root() ******************************/
/** finds the cluster root for a selected spin                **/
/** see: M. E. J. Newman and R. M. Ziff. Fast Monte Carlo     **/
/** algorithm for site or bond percolation.                   **/
/** Phys. Rev. E, 64, 016706, (2001).                         **/
/**                                                           **/
/** PARAMETERS                                                **/
/**        i: ith spin                                        **/
/**      clu: list of spin roots                              **/
/**                                                           **/
/** RETURNS:                                                  **/
/**        r: spin root                                       **/
/***************************************************************/


int find_root(int i, int *clu){
  int r,s;
  r = s = i;
  while (clu[r] >= 0) {
    clu[s] = clu[r];
    s = r;
    r = clu[r];
  }
  return r;
}

#endif /* _FIND_ROOT_H_ */