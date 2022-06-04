#ifndef _ALIAS_TABLE_H_
#define _ALIAS_TABLE_H_

/********************* alias_table() ***************************/
/** finds the cluster root for a selected spin                **/
/** see: D. E. Knuth. The Art of Computer Programming,        **/
/** Volume 2 (3rd Ed.): Seminumerical Algorithms.             **/
/** Addison-Wesley Longman Publishing Co., Inc.,              **/
/** Boston, MA, USA, (1997).                                  **/
/**                                                           **/
/** PARAMETERS                                                **/
/**        L: ith spin                                        **/
/**        a: relations of bonds                              **/
/**       bp: bond selection probability                      **/
/**     Jtot: sum of all interaction couplings                **/     
/**     dist: distance between spins ith and jth              **/  
/**                                                           **/
/** RETURNS:                                                  **/
/**        update the a and bp                                **/
/***************************************************************/


void alias_table(int L, int *a, double *bp, double Jtot, double *dist)
{
  int cal,coun,al,i;
  int rn[L];
  double p[L],re[L]; 
  int num = 0;
  for (i = 0; i < L-1; ++i) {
   p[i] = dist[i]/Jtot;
   bp[i] = p[i]*L*(L-1)/2;
   num += 1;
   rn[i] = 0;
   re[i] = 0;
  }
  cal = -1;
  coun = num-1;
  for (i = 0; i < L-1; ++i) {
    if ( bp[i] >= 1) {
      cal += 1;
      re[cal] = bp[i];
      rn[cal] = i;
    }
  }
  for (i = L-2; i > (-1); --i) {
    if ( bp[i] < 1 ) {
      cal += 1;
      coun -= 1;
      re[cal] = bp[i];
      rn[cal] = i;
    }
  }
  al = num-1;
  while ((al > -1 && coun > -1) && (coun != al)) {
    a[rn[al]] = rn[coun];
    re[coun] = re[coun]-1+re[al];
    al -= 1;
    if (re[coun] < 1) {
      coun -= 1; 
    }
  }
  for(i=0;i<num;++i) bp[rn[i]]=re[i];
}

#endif /* _ALIAS_TABLE_H_ */