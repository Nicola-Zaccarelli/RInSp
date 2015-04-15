# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <R.h>
# include <Rinternals.h>
# include <Rdefines.h>
# include <Rmath.h>
// The procedure replicates the two strategies of Monte Carlo resampling available
// in Emc.c, PSicalc.c, WTcMC.c, and WTdMC.c.
//
// Author: Dr. Nicola ZACCARELLI (nicola.zaccarelli@gmail.com)
//
// Version 1.3
// Date: 2015-02-03

SEXP MCprocedure(SEXP resources, SEXP MCType, SEXP DietPop, SEXP nreplicates)
{

// Declare all variables
int   NRows, NCols, i, j, sim, nreps, typeMC;
double cumulativep, lowerbound, item;
double temp, newi, newj, x;
double *rdata, *risult, *popdiet, *totaldieti, *Nitems;
double **mcdata, **realdata;
SEXP Rris, Rdim; 

Rdim = getAttrib(resources, R_DimSymbol);
NRows = INTEGER(Rdim)[0];
NCols = INTEGER(Rdim)[1];

// Coerce data matrix into vector for data access form C
PROTECT(resources = coerceVector(resources, REALSXP));
rdata = REAL(resources);

// Create relevant objects
MCType = coerceVector(MCType, INTSXP);
typeMC = INTEGER(MCType)[0];
DietPop = coerceVector(DietPop, REALSXP);
popdiet = REAL(DietPop);
Nitems = calloc(NRows, sizeof(double *));

// Create the matrix for Monte Carlo resempling
realdata = calloc(NRows, sizeof(double *));
for (i=0; i < NRows; i++) realdata[i] = calloc(NCols, sizeof(double));

// Create the matrix for Monte Carlo resempling
mcdata = calloc(NRows, sizeof(double *));
for (i=0; i < NRows; i++) mcdata[i] = calloc(NCols, sizeof(double));
totaldieti = calloc(NRows, sizeof(double));

// Attention to the right address for the elements
for (i=0; i<NCols; i++) {
  for (j=0; j<NRows; j++){
  realdata[j][i] = rdata[j + NRows*i];}
  }
// Create the return matrix
nreplicates = coerceVector(nreplicates, INTSXP);
nreps = INTEGER(nreplicates)[0];
PROTECT(Rris = allocMatrix(REALSXP, (NRows * NCols), nreps));
risult = REAL(Rris);

// read in (or create) .Random.seed for R random generation fuctions
GetRNGstate();

// Calculate individual diet or number of items from real data
for (i=0; i < NRows; i++)
     {
	 totaldieti[i] = 0;
	 Nitems[i] = 0;
	 for (j=0; j < NCols; j++) {
	     totaldieti[i] = totaldieti[i] + realdata[i][j];
		 if (realdata[i][j] > 0) Nitems[i] = Nitems[i] + 1;
		 }
	}

// Monte Carlo simulation and writing of results
for (sim = 0; sim < nreps; sim++)
{
// initialize MC matrix
for (i=0; i<NRows; i++) {
  for (j=0; j<NCols; j++){
   mcdata[i][j] = 0;
}}
// integer data type
if (typeMC == 1)
{
for (i=0; i<NRows; i++)
     {for (x=0; x< totaldieti[i]; x++)
         {item =  unif_rand();    // Using R random function from Rmath.h
          cumulativep = 0;
          for (j=0; j<NCols; j++)
             {lowerbound = cumulativep;
              cumulativep = cumulativep + popdiet[j];
              if (item>=lowerbound && item<cumulativep)
                 {mcdata[i][j] = mcdata[i][j] + 1;}}
          } 
      };
} else 
{
for (i=0; i<NRows; i++)
     {
      for (j=0; j<Nitems[i]; j++)
         {
         temp = 0;
         while (temp == 0)
             {
              newi = (double)NRows * unif_rand();    // Using R random function from Rmath.h
              newj = (double)NCols * unif_rand();    // Using R random function from Rmath.h
              temp = realdata[(int)newi][(int)newj];
              }
          mcdata[i][j] = temp;
         }
     }
}

for (i=0; i<NCols; i++)
     {
      for (j=0; j<NRows; j++)
       { risult[(NCols*NRows*sim)+ (i*NRows) + j] = mcdata[j][i];}
     }

} // for of sim

// write .Random.seed out after use
PutRNGstate();

UNPROTECT(2);
free(totaldieti);
return Rris;
}
