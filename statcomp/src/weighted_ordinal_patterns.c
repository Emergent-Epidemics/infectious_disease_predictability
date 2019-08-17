#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rdefines.h>
void weighted_ordinal_pattern_loop(double *x, int *nx, int *ndemb, double *ifrec, int *nifrec)
{
    int ipa[*ndemb][*ndemb];
    int nd;
    int has_na;
    double xvec[*ndemb];
    double epsilon=1e-10;
    double average, variance, sum = 0, sum1 = 0;
    //if (ISNA(x[0])) x[0]=0;
    //x[2]=10;
    //ifrec[0]=1;
    for (int nv=0;nv<*nx-*ndemb+1;nv++) {
      has_na=0;
      for (int ixvec=0;ixvec<*ndemb;ixvec++) {
        if (ISNA(x[nv+ixvec])) {has_na=1; break;}
        xvec[ixvec]=x[nv+ixvec]; /* Fill xvec */
      }
      if (has_na) continue;
      for (int i=0;i<*ndemb;i++) for (int j=0;j<*ndemb;j++) ipa[i][j]=0; /* Reset ipa matrix */

      for (int ilag=1;ilag<*ndemb;ilag++) {
        for (int itime=ilag;itime<*ndemb;itime++) {
          ipa[itime][ilag] = ipa[itime-1][ilag-1];
          if ((xvec[itime] <= xvec[itime - ilag] ) || ( fabs( xvec[itime - ilag] - xvec[itime]) < epsilon))
            ipa[itime][ilag] = ipa[itime][ilag] + 1;
        }
      }
      nd = ipa[*ndemb-1][1];
      for (int ilag=2;ilag<*ndemb;ilag++) {
        nd =(ilag+1) * nd + ipa[*ndemb-1][ilag];
      }
        /* compute variance of the pattern */
        sum=0;
        for (int i = 0; i < *ndemb; i++) sum = sum + xvec[i];
        average = sum / *ndemb;
        /*  Compute  variance:  */
        sum1=0;
        for (int i = 0; i < *ndemb; i++) {
            sum1 = sum1 + pow((xvec[i] - average), 2);
        }
        variance = sum1 / *ndemb;
        
       ifrec[nd] = ifrec[nd] + variance;
        /* ifrec[nd] = ifrec[nd] + 1; */
        
        /* *Nweighted = *Nweighted + variance; N_weighted: Calculates denominator of Eq. 4 in Fadlallah et al (2013) */ 
    }
}


// static R_NativePrimitiveArgType ordinalArgs_t[] = {
//    REALSXP, INTSXP, INTSXP, INTSXP,INTSXP
// };

// static R_CMethodDef cMethods[] = {
//   {"weighted_ordinal_pattern_loop", (DL_FUNC) &weighted_ordinal_pattern_loop, 5, ordinalArgs_t},
//   {NULL, NULL, 0}
// };

// void R_init_mylib(DllInfo *info)
// {
//  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
// }

// void R_unload_mylib(DllInfo *info)
// {
      /* Release resources. */
// }
