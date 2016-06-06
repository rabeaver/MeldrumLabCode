#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"

/*
    calling syntax: [spectrum, compte] = TwoDLaplaceMainLoop(E, indice, spectrum, data, W)
*/

/* Input Arguments */

#define E_IN prhs[0]
#define indice_IN prhs[1]
#define spectrum_IN prhs[2]
#define data_IN prhs[3]
#define W_IN prhs[4]

/* Output Arguments */

#define spectrum_OUT plhs[0]
#define compte_OUT plhs[1]

int do_calculations(int, double *, int *, int, double *, double *, double *);
void qrdecmp(double *, int *, int *, int *, double *);
void lineareqsolve(double *, int *, int *, double *);
void deletefromset(int *, int *, int, int);
double house(int, int, double *, int *, double *);
void orthomult(double *, int *, double *, double);
int findlambda(double *, double *, double *, int *, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *spectrum, *spectrum2, *E, *data, *W;
	int E_dim[2];
	int compte, indice, Xdim;
    
/* Check for proper number of arguments */
  
	if (nrhs != 5) 
	{
		mexErrMsgTxt("This program requires FIVE input arguments.");
	} 
	else if (nlhs > 2) 
	{
		mexErrMsgTxt("This program has only TWO output arguments.");
  	}

/* Assign the various parameters */
  
  	E = mxGetPr(E_IN);
  	E_dim[0] = mxGetM(E_IN);
  	E_dim[1] = mxGetN(E_IN);
	Xdim=E_dim[1]-1;
  	indice = ((int) *mxGetPr(indice_IN))-1;
   spectrum = mxGetPr(spectrum_IN);
   data = mxGetPr(data_IN);
   W = mxGetPr(W_IN);
  	
   compte = do_calculations(Xdim, E, E_dim, indice, spectrum, data, W);
   spectrum_OUT = mxCreateDoubleMatrix(1,Xdim,mxREAL);
	spectrum2 = mxGetPr(spectrum_OUT);
	memcpy(spectrum2,spectrum,Xdim*sizeof(double));
	compte_OUT = mxCreateDoubleScalar(compte);
    /* Replaced mxCreateScalarDouble with mxCreateDoubleScalar on line above, TKM 31 Mar 2016. Compiles. */
}

int do_calculations(int Xdim, double *E, int *E_dim, int indice, double *spectrum, double *data, double *W)
{
	double *soln, *Es, *answ, lambda[1], *u;
	int *setp, *oldsetp, *setz;
	double biggest;
	int compte = 0, i, j, newindice, fini=0, repeat=1, ind, indextozero;
	
	if((setp=mxCalloc(Xdim,sizeof(int)))==NULL) mexErrMsgTxt("Not enough memory for setp.\n");
	if((setz=mxCalloc(Xdim,sizeof(int)))==NULL) mexErrMsgTxt("Not enough memory for setz.\n");
	if((oldsetp=mxCalloc(Xdim,sizeof(int)))==NULL) mexErrMsgTxt("Not enough memory for oldsetp.\n");
	if((soln=mxCalloc(Xdim,sizeof(double)))==NULL) mexErrMsgTxt("Not enough memoryfor soln.\n");
	if((Es=mxCalloc(E_dim[0]*(E_dim[1]-1),sizeof(double)))==NULL) mexErrMsgTxt("Not enough memory for Es.\n");
	if((answ=mxCalloc(E_dim[0],sizeof(double)))==NULL) mexErrMsgTxt("Not enough memoryfor answ.\n");
	if((u=mxCalloc(E_dim[0],sizeof(double)))==NULL) mexErrMsgTxt("Not enough memory for u.\n");

	memcpy(Es,E,E_dim[0]*(E_dim[1]-1)*sizeof(double));
	memset(soln,0,Xdim*sizeof(double));
	memset(setp,0,Xdim*sizeof(int));
	for(i=0;i<Xdim;i++)
		setz[i]=1;
	
	while(fini == 0)
	{
		compte++;
/* step #5 of the L&H algorithm */
/* move the index "indice" from the set setz to the set setp */
		memcpy(oldsetp, setp, Xdim*sizeof(int));
		ind=0;
		while(ind<Xdim)
			if(setp[ind]!=0)
				ind++;
			else
				break;
		setp[ind]=indice;
		setz[indice]=-1;
		
		qrdecmp(E,E_dim,setp,oldsetp,u);
		lineareqsolve(E,E_dim,setp,soln);
		if(soln[indice] < 0) mexWarnMsgTxt("correction made for roundoff error.\n");
/* only if soln[indice] < 0, find new indice */
		while(soln[indice] < 0)
	   {
	   	W[indice] = 0;
	   	biggest = 0;
	   	newindice = 0;
	   	for(i = 0;i < Xdim;i++)
	   		if((W[i] > biggest)&&(setz[i] == 1))
	   		{
	   			biggest = W[i];
	   			newindice = i;
	   		}
	   	if(newindice == 0)
	   	{
				fini = 1;
				break;
			}
			setz[indice]=1;
			setz[newindice]=-1;
			indice = newindice;
			ind=0;
			while(ind<Xdim)
				if(setp[ind]!=0)
					ind++;
				else
					break;
			setp[ind-1]=indice;
			
			qrdecmp(E,E_dim,setp,oldsetp,u);
			lineareqsolve(E,E_dim,setp,soln);
		}

		if(fini) break;
		repeat = 1;
	    
/* step #6 L&H */
/* step #7: if Zj >0 for all j in setp, set x=z (spectrum = soln) and go to step #2, else ... */
		while(repeat)
		{
			indextozero = findlambda(lambda, spectrum, soln, setp, Xdim);
/* step #8 #9. find indextozero=q in setp, such that xq/(xq-zq)=min. */
/* If there is a negative value in soln, lambda won't be zero */
			if(indextozero != 0)
			{
/* if there is a negative value in soln, we have to keep removingindices from setp until all elt in soln are positive */
				if(lambda[0] == 0)     /* not necessary */
				{
					mexWarnMsgTxt("algorithm fault - infinite loop.\n");
					repeat = 0;
					fini = 1;
					break;
				}
/* step #10 set x = x + lambda(z-x) */
				for(i=0;i<Xdim;i++)
					spectrum[i] += lambda[0]*(soln[i] - spectrum[i]);
/* steps #11. Move from setp to setz all indices j in setp for which xj=0. go to step #6 */
				memcpy(oldsetp, setp, Xdim*sizeof(int));
				deletefromset(setp,setz,indextozero,Xdim);
				spectrum[indextozero] = 0;
				soln[indextozero] = 0;
/* step # 6 */
				qrdecmp(E,E_dim,setp,oldsetp,u);     							/*  Get to step 6 in the algorithm */
				lineareqsolve(E,E_dim,setp,soln);                      /* until no more soln is negative */
			}
			else 
/* Step #7 : x=z, go to step #2 */
			{
				repeat = 0;
				memcpy(spectrum, soln, Xdim*sizeof(double));
			}
		}
		if(fini == 1) break;
		
/* step #2 of the L&H algorithm. first step: spectrum is all zero */
		memcpy(answ,data,Xdim*sizeof(double));
		memset(answ+Xdim,0,(E_dim[0]-Xdim)*sizeof(double));
		for(i=0;i<E_dim[0];i++)
			for(j=0;j<Xdim;j++)
				answ[i] -= Es[j*E_dim[0]+i]*spectrum[j];
/* w = E'(f-Ex) */
		memset(W,0,Xdim*sizeof(double));
		for(i = 0; i < Xdim; i++)
			for(j = 0; j < E_dim[0]; j++)
				W[i] += Es[j+i*E_dim[0]]*answ[j];
		
/* step #4 of the L&H algorithm */
/* find index of bigger number in W, whose index is still in setz */
	    indice = -1;
	    biggest = 0;
	    for(i = 0; i < Xdim; i++)
	        if((W[i] > biggest)&&(setz[i] == 1))
	        {
	            biggest = W[i];
	            indice = i;
				}
	    
	    if(indice == -1)
	        fini = 1;
	    if((fini != 1)&&(setp[Xdim-1] != 0))
	        fini = 1;
	}
	mxFree(setp);
	mxFree(setz);
	mxFree(oldsetp);
	mxFree(soln);
	mxFree(Es);
	mxFree(answ);
	mxFree(u);
	return compte;
}

void qrdecmp(double *E, int *E_dim, int *setp, int *oldsetp, double *u)
{
/* series of Householder tranformation designed specifically for the LS problem */
	int setlength = 0, added = 1, removed;
	double b;
	
	while((setlength<E_dim[1]-1)&&(setp[setlength]!=0)&&(oldsetp[setlength]!=0))
			setlength++;
/* find a zero to either setp or oldsetp */
	
	if(setp[setlength]==0)
		added = 0;
	
	if(added == 1)
/* householder transformation to most recently chosen column, setp[setlength] in the matrix */
	{
		b = house(setp[setlength],setlength,E,E_dim,u);
		orthomult(E,E_dim,u,b);
	}
	else
/* householder transformation to all columns listed in setp */
	{
		removed = 0;
		while((removed<E_dim[1]-1)&&(setp[removed] == oldsetp[removed]))
			removed++;
		while((removed<E_dim[1]-1)&&(setp[removed] != 0))
		{
				b = house(setp[removed],removed,E, E_dim,u);
				orthomult(E,E_dim,u,b);
				removed++;
		}
	}
}

double house(int column, int pivot, double *E, int *E_dim, double *u)
{
/* rotates a vector from "mat" (indexed by column) */
/* the outputs "u" (vector) and "b" (float) allow us to apply orthogonal tranform elsewhere */
	double som, s, b;
	int i;

	memset(u,0,pivot*sizeof(double));
	memcpy(u+pivot+1,E+pivot+1+column*E_dim[0], (E_dim[0]-pivot-1)*sizeof(double));
	som = 0;
	for(i=pivot+1; i<E_dim[0]; i++)
		som += u[i]*u[i];
	if(E[pivot+column*E_dim[0]]<0)
		s = sqrt(som + E[pivot+column*E_dim[0]]*E[pivot+column*E_dim[0]]);
	else
		s = - sqrt(som + E[pivot+column*E_dim[0]]*E[pivot+column*E_dim[0]]);
	u[pivot] = E[pivot+column*E_dim[0]]-s;
	b = s*u[pivot];
	return b;
}

void orthomult(double *E, int *E_dim, double *u, double b)
{
/* function used to apply the householder transformation to each column in some "mat" */
	int i, j;
	double c;
	
	for(i = 0; i < E_dim[1]; i++)
	{
		c = 0;
		for(j=0; j<E_dim[0]; j++)
			c += u[j]*E[i*E_dim[0]+j]/b;
		for(j=0; j<E_dim[0]; j++)
			E[i*E_dim[0]+j] += c*u[j];
	}
}

void lineareqsolve(double *E, int *E_dim, int *setp, double *soln)
{
	double som;
	int i, j, setlength=0;

	memset(soln,0,(E_dim[1]-1)*sizeof(double));
	while((setlength<E_dim[1]-1)&&(setp[setlength]!=0))
			setlength++;
	setlength--;

	soln[setp[setlength]] = E[setlength+(E_dim[1]-1)*E_dim[0]]/E[setlength+setp[setlength]*E_dim[0]];
	for(i = setlength-1; i >= 0; i--)
	{
		som = 0;
		for(j=setlength; j > i; j--)
			som += E[i+setp[j]*E_dim[0]]*soln[setp[j]];
		soln[setp[i]] = (E[i+(E_dim[1]-1)*E_dim[0]]-som)/E[i+setp[i]*E_dim[0]];
	}
}

void deletefromset(int *setp, int *setz, int index, int Xdim)
{
	int i;
	
	setz[index] = 1;
	i = 0;
	while((i<Xdim)&&(setp[i] != index))
			i++;
	memmove(setp+i,setp+(i+1),(Xdim-(i+1))*sizeof(int));
	setp[Xdim-1]=0;
}

int findlambda(double *lambda, double *spectrum, double *soln, int *setp, int Xdim)
{
	int i, zeroindice;
	
/* step #8 */
/* Find an indice q (zeroindice) in setp such that xq/(xq-zq)=min */
/* x=spectrum and z = soln (solution of step 6) */
	zeroindice=0;
	lambda[0]=0.0;
	for(i = 0; i < Xdim; i++)
		if(setp[i] != 0)
			if(soln[setp[i]]<0.0)
			{
/* run through all indices in setp,  */
/* calculate the lambda value for indices that correspond  */
/* to negative coefficients in soln, see which lambda is smallest */
				if(lambda[0] == 0.0)
				{
/* for the first time in the loop */
					lambda[0] = spectrum[setp[i]]/(spectrum[setp[i]]-soln[setp[i]]);
					zeroindice = setp[i];
				}
				if(spectrum[setp[i]]/(spectrum[setp[i]]-soln[setp[i]]) < lambda[0])
				{
/* successive loops compared to current lambda */
					lambda[0] = spectrum[setp[i]]/(spectrum[setp[i]]-soln[setp[i]]);
					zeroindice = setp[i];
				}
			}
	return zeroindice;
}
