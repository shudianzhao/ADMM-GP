/*
function [T, b, hash, gamma] = tri_sep_kc( y, n, maxtri, h_tmp);
return most violated triangular ineq. for k-cut problem
input: matrix y (n,n) spaltenweise als vektor /columnwise for vectors
maxtri .. gewuenschte anzahl an dreiecken / count for tri.
(maxtri=0) -> nur check ob Y dreiecke verletzt /only checK if there is violated tri. in Y
h_tmp...sortiertes altes hash-array /sorted old hash array
output: [T,b,hash] data structure of new triangles
gamma: violations
*/

/*
- This function only works when n<=1000 due to the structure of hash values
- y in [0, 1] y = vec(X)

- type 0: x_ik + x_jk <= x_ij +1
- type 1: x_ij + x_jk <= x_ik +1
- type 2: x_ik + x_ij <= x_jk +1

*/



#include <math.h>
#include <stdio.h>
#include "mex.h"
/*#include "heapsort.h"*/

#define MAX_INEQ 20000

#define COMPARE(x,y) (((x)<(y)) ? -1: ((x)==(y))? 0: 1)


/*
* function prototypes
*/

int binsearch(unsigned long list[], unsigned long num, int left, int right);


/*
* subroutines
*/

void heapify_n(i,n,ind,d) // generate a min heap, item at root has the minimum value
int i,n,*ind;
double *d;
{
	int k,h,j;

	ind--; i++;
	k=i;
	while(2*k<=n){
		if (2*k==n){
			if (d[ind[k]]>d[ind[2*k]]){
				h=ind[k];
				ind[k]=ind[2*k];
				ind[2*k]=h;
			}
			break;
		}
		if ((d[ind[k]]>d[ind[2*k]])||
		(d[ind[k]]>d[ind[2*k+1]])){
			if (d[ind[2*k]]<d[ind[2*k+1]]) j=2*k;
			else j=2*k+1;
			h=ind[k];
			ind[k]=ind[j];
			ind[j]=h;
			k=j;
		}
		else break;
	}
}


int binsearch(unsigned long list[],unsigned long num, int left, int right)
{

	/* search num through list[0]<=list[1]<=...<=list[n-1]
	return 1 if its found ; return 0 otherwise
	*/
    
	int middle;
	while (left <= right) {
		middle=floor((left+right)/2);
		switch(COMPARE(list[middle],num))
		{ case -1: left = middle+1;break;
			case 0 : return 1;right=0;
			case 1 : right=middle-1;
		}
	}
	return 0;
}



void mexFunction( int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
	{
		double *y, *n, *maxtri, *T, *b, *gamma;
        double *h_tmp,*hash; 
		int m, l, vmax, ftest=0;

		/* right input has 3 rhs arguments and 4 lhs arguments */

		if (nlhs != 4 || nrhs != 4)
		{
			mexErrMsgTxt("input-output arguments inconsistent.");
		}

		y = mxGetPr(prhs[0]);
		n = mxGetPr(prhs[1]);
		maxtri = mxGetPr(prhs[2]);
		h_tmp = mxGetPr(prhs[3]);
		m=mxGetM(prhs[3]);//  the number of rows in h_tmp

    

		vmax = maxtri[0]; l = n[0];
		//   vmax: number of triangle inequalites to be returned
		//  l: the dimension of the matrix
		if (vmax<1) {ftest = 1; vmax = 1;}
		//   only test for violated inequalities
		if (vmax>MAX_INEQ) {vmax = MAX_INEQ;}
		if (mxGetM(prhs[0]) != l*l) {
			mexErrMsgTxt("y and n inconsistent.");
		}

		{
			double threshold,x, yij, yik, yjk, yi,yj,yk,yy,violl[MAX_INEQ];
			int ind[MAX_INEQ], i, j, k, ii, cntr;
			int T_tmp[MAX_INEQ*4], b_tmp[MAX_INEQ];

			unsigned long  h_new[MAX_INEQ];

			int flag,i1,typ, check;
			unsigned long  hash_curr_old;
            unsigned long hash_curr;
			unsigned long *h_tmpC;

            
           
			/* dynamic memory allocation */
			h_tmpC=(unsigned long*)mxCalloc(m,sizeof(unsigned long));
            
     

			for (j=0;j<m;j++)                  /* move h_tmp to h_tmpC */
			h_tmpC[j]=(unsigned long)h_tmp[j]; /* (unsigned-long-array) */


			for (i=0; i<vmax; i++) {
				violl[i] = -1.0;
				ind[i] = i;
			}


            

            cntr = 0;
			threshold = .001;
			if (l>300)	//  large thershold for large matrix ?
			threshold = .001;
			for (i=1; i<=l-2; i++)
			{ for (j=i+1; j<=l-1; j++)
				{ for (k=j+1; k<=l; k++)
					{ yij = y[(i-1)*l+j-1];
						yik = y[(i-1)*l+k-1];
						yjk = y[(j-1)*l+k-1];

						/* Typ 0 */
						x =  -(yij - yik - yjk) -1 ; // if x > 0, constr. violated
						typ=0;
						hash_curr = (unsigned long)((typ*1000+i)*1000+j)*1000+k;
						//       check data type
						if (m<1) flag=0;
						else if (hash_curr<h_tmpC[0]) flag=0;
						else if (hash_curr>h_tmpC[m-1]) flag=0;
						else if (x>threshold)
						flag=binsearch( h_tmpC, hash_curr, 0, m-1);
						else
						flag=0;

						if ((x>threshold) && (flag<1))
						{ if (x>violl[ind[0]])
							{ ii = ind[0];
								violl[ii] = x;
								T_tmp[ii] = i; T_tmp[vmax+ii] = j;
								T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = typ;
								h_new[ii] = hash_curr;
								b_tmp[ii] = 0;
								//            maybe is the rhs? need?
								heapify_n(0,vmax,ind,violl);
								cntr++;
							}
						}
						/* Typ 1 */
						x =  -(- yij + yik - yjk)-1;
						typ=1;
						hash_curr = (unsigned long) ((typ*1000+i)*1000+j)*1000+k;
                        
						if (m<1) flag=0;
						else if (hash_curr<h_tmpC[0]) flag=0;
						else if (hash_curr>h_tmpC[m-1]) flag=0;
						else if (x>threshold)
						flag=binsearch( h_tmpC, hash_curr, 0, m-1);
						else
						flag=0;

						if ((x>threshold) && (flag<1))
						{ if (x>violl[ind[0]])
							{ ii = ind[0];
								violl[ii] = x;
								T_tmp[ii] = i; T_tmp[vmax+ii] = j;
								T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = typ;
								h_new[ii] = hash_curr;
								b_tmp[ii] = 0;
								heapify_n(0,vmax,ind,violl);
								cntr++;
							}
						}
						/* Typ 2 */
						x =  -(- yij - yik + yjk)-1;
						typ=2;
						hash_curr = (unsigned long)((typ*1000+i)*1000+j)*1000+k;
						if (m<1) flag=0;
						else if (hash_curr<h_tmpC[0]) flag=0;
						else if (hash_curr>h_tmpC[m-1]) flag=0;
						else if (x>threshold)
                            flag=binsearch( h_tmpC, hash_curr, 0, m-1);
						else
						flag=0;
                
						if ((x>threshold) && (flag<1))
						{ if (x>violl[ind[0]])
							{
                              if (hash_curr > 4294967295) {
                                  printf("hash value %lu \n", hash_curr);
                                  printf("hash value is supposed to be %lu \n",((typ*1000+i)*1000+j)*1000+k );
                                  hash_curr_old = 4294967295;
                                  printf("hash value beyond unsigned long limit %lu \n",hash_curr_old); 
                                  goto ende;
                              }
                              ii = ind[0];
								violl[ii] = x;
								T_tmp[ii] = i; T_tmp[vmax+ii] = j;
								T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = typ;
//                                 printf("hash value %lu \n", hash_curr);
//                                 hash_curr_old = ((typ*1000+i)*1000+j)*1000+k;
//                                 printf("((typ*1000+i)*1000+j)*1000+k: %lu \n", hash_curr_old);
								h_new[ii] = hash_curr;
								b_tmp[ii] = 0;
								heapify_n(0,vmax,ind,violl);
								cntr++;
							}
						}

						if (ftest>0 && cntr >0) goto ende;
						// ftest: if there is anything violated.
					}
				}
			}

			ende:
  
            
                
			if (cntr > vmax) cntr= vmax;  /* and 4 lhs arguments            */
			plhs[0] = mxCreateDoubleMatrix(cntr*4, 1, mxREAL);
			plhs[1] = mxCreateDoubleMatrix(cntr, 1, mxREAL);
			plhs[2] = mxCreateDoubleMatrix(cntr, 1, mxREAL);
			plhs[3] = mxCreateDoubleMatrix(cntr, 1, mxREAL);

			T = mxGetPr(plhs[0]);
			b = mxGetPr(plhs[1]);
			hash = mxGetPr(plhs[2]);
			gamma = mxGetPr(plhs[3]);

			i = 0;
			for (k=0; k<vmax; k++) {
				j = ind[k];
				if (violl[j]>-1.) {
					gamma[i] = violl[j];
					hash[i] = h_new[j];
//                     printf("hash value %lf \n", hash[i]);
					T[i] = (double) T_tmp[j];
					T[cntr+i] = (double) T_tmp[vmax+j];
					T[2*cntr+i] = (double) T_tmp[vmax*2+j];
					T[3*cntr+i] = (int) T_tmp[vmax*3+j];
					b[i] = b_tmp[j];
					i++;}
				}
			}
		}
