# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <R.h>
# include <Rinternals.h>
# include <Rdefines.h>
# include <Rmath.h>
// Calculate the proportional similarity index Psi
// Bolnick et al. (2003) describe another measure of individual specialization, based on
// he average pairwise overlap of the niche distribution of individuals and the population.
//
// Author: Dr. Nicola ZACCARELLI (nicola.zaccarelli@gmail.com)
//
// Version 1.3
// Date: 2015-02-03

SEXP PSicalc(SEXP comcalc, SEXP popdiet, SEXP nreplicates)
{
// Declare all variables
int   NRows, NCols, i, j, nreps, popD, k, l, Rep;
double IS, bspoptotal;
double item, cumulativep, lowerbound, bsresourcejtotal;
double ps, t1, t2, t3, t4;
double *dati, *results, *totaldieti, *populationdiet;
double *bspopulationdiet, *resfin;
double **bsdata, **bsproportions;
double *PSresults, *VARps, *plessthanq, *qlessthanp;

// Rsin is a matrix containing results with the first replicate as the actual
// result from the data matrix
SEXP Rris, Rdim;
Rdim = getAttrib(comcalc, R_DimSymbol);
NRows = INTEGER(Rdim)[0];
NCols = INTEGER(Rdim)[1];

// Coerce data matrix into vector for data access form C
PROTECT(comcalc = coerceVector(comcalc, REALSXP));
dati = REAL(comcalc);

// Create objects and return matrix from R
popdiet = coerceVector(popdiet, INTSXP);
popD = INTEGER(popdiet)[0];
nreplicates = coerceVector(nreplicates, INTSXP);
nreps = INTEGER(nreplicates)[0];
PROTECT(Rris = allocMatrix(REALSXP, (2 * NRows + 1), (nreps + 1)));
resfin = REAL(Rris);

// Let's initialize the final vector
for (i = 0; i < ((2 * NRows + 1)*(nreps + 1)); i++) { resfin[i] = 0; }

// Create the matrix for bootstrapping
bsdata = calloc(NRows, sizeof(double *));
for (i=0; i < NRows; i++) bsdata[i] = calloc(NCols, sizeof(double));

// Creation of necessary matrices in C
bsproportions = calloc(NRows, sizeof(double *));
for (i=0; i < NRows; i++) bsproportions[i] = calloc(NCols, sizeof(double));
bspopulationdiet = calloc(NCols, sizeof(double));
results = calloc(NRows, sizeof(double));
totaldieti = calloc(NRows, sizeof(double));
populationdiet = calloc(NCols, sizeof(double));
plessthanq = calloc(NCols, sizeof(double));
qlessthanp = calloc(NCols, sizeof(double));

for (j=0; j < NCols; j++)
  { populationdiet[j] = 0;
    bspopulationdiet[j] = 0;
    plessthanq[j] = 0;
    qlessthanp[j] = 0;
  }

VARps = calloc(NRows, sizeof(double));
PSresults = calloc(NRows, sizeof(double));

// read in (or create) .Random.seed for R random generation fuctions
GetRNGstate();

// Attention to the right address for the elements
for (j=0; j<NCols; j++)
  { for (i=0; i<NRows; i++)
      { bsdata[i][j] = (double)dati[i + NRows*j];
        bsproportions[i][j] = 0; }
   }

// Calculate individual diet and population diet

bspoptotal = 0;

for (i=0; i < NRows; i++)
  { totaldieti[i] = 0;
    // bstotaldieti[i] = 0;
    for (j=0; j < NCols; j++)
      { totaldieti[i] = totaldieti[i] + bsdata[i][j]; }
      bspoptotal = bspoptotal + totaldieti[i];
   }

// Start with calculate real data and than to simulate
for (Rep = 0; Rep < (nreps + 1); Rep++)
  {
    if (Rep > 0)
    {
      for (i=0; i<NRows; i++) { for (j=0; j<NCols; j++) { bsdata[i][j] = 0; } }
      for (i=0; i<NRows; i++)
        { for (k=0; k< totaldieti[i]; k++)
            { item = unif_rand();    // Using R random function from Rmath.h
              cumulativep = 0;
              for (j=0; j<NCols; j++)
                { lowerbound = cumulativep;
                  cumulativep = cumulativep + populationdiet[j];
                  if (item>=lowerbound && item<cumulativep) { bsdata[i][j]=bsdata[i][j] + 1; }
                } // end for j
            } // for x
        }; // for i
    } // end if (R>0)
// As the marginal total for rows are fixed the proportions calculation is easier
   for (i=0; i < NRows; i++)
     { for (j=0; j < NCols; j++)
         { bsproportions[i][j] = (bsdata[i][j] / totaldieti[i]); }
     }
// CALCULATE POPULATION DIET
   switch (popD)
        { case 0 :
               {
                for (j = 0; j < NCols; j++)
                  {
                    bsresourcejtotal = 0;
                    for (i = 0; i<NRows; i++) { bsresourcejtotal = bsresourcejtotal + bsdata[i][j];}
                    bspopulationdiet[j] = ( bsresourcejtotal/bspoptotal );
                    if (Rep == 0) { populationdiet[j] = bspopulationdiet[j];}
                  }
           break;}
           case 1:
                {
                 for (j = 0; j<NCols; j++)
                   { bsresourcejtotal = 0;
                     for (i=0; i<NRows; i++) { bsresourcejtotal = bsresourcejtotal + bsproportions[i][j]; }
                     bspopulationdiet[j] = (bsresourcejtotal / NRows);
                     if (Rep == 0) { populationdiet[j] = bspopulationdiet[j]; }
                    }
           break;}
        } // end switch

    for (i=0; i< NRows; i++) { results[i] = 0; }

    for (i = 0; i<NRows; i++)
      { 
        PSresults[i] = 0;
        VARps[i] = 0;
        ps = 0;
// CALCULATE PSi
        for (j = 0; j<NCols; j++)
          { if (bsproportions[i][j] <= bspopulationdiet[j])
             { plessthanq[j] = bsproportions[i][j];
               qlessthanp[j] = 0; }
            else
             { plessthanq[j] = 0;
               qlessthanp[j] = bspopulationdiet[j];}
          }
        for (j = 0; j<NCols; j++) { ps = ps + plessthanq[j] + qlessthanp[j]; }
        PSresults[i] = ps;
// CALCULATE VAR(PSi)
        t1 = 0;
        t2 = 0;
        t3 = 0;
        t4 = 0;
        for (k = 0; k<NCols; k++)
         {
          t1 = t1 + (plessthanq[k]*(1 - plessthanq[k]));
          t3 = t3 + (qlessthanp[k]*(1 - qlessthanp[k]));
         }
        for (k = 0; k<NCols; k++)
        { for (l = (k+1); l<NCols; l++)
            { t2 = t2 + (2 * (plessthanq[k] * plessthanq[l]));
              t4 = t4 + (2 * (qlessthanp[k] * qlessthanp[l]));
            }
        }
        VARps[i] = ( (1 / totaldieti[i]) * (t1 - t2) + ( 1 / bspoptotal)*( t3 - t4));
    } // end Claculation of Var(PSi)
// CALCULATE IS
    IS = 0;
    for (i=0; i<NRows; i++) { IS = IS + PSresults[i]; }
    IS = IS / NRows;
// Write results
     for (k = 0; k < NRows; k++)
      { 
        resfin[Rep * (2*NRows + 1) + k] = PSresults[k];
        resfin[Rep * (2*NRows + 1) + NRows + k] = VARps[k]; }
    resfin[Rep * (2*NRows + 1) + 2 * NRows] = IS;
	} // End For Rep
UNPROTECT(2);
free(results);
free(totaldieti);
free(populationdiet);
free(bspopulationdiet);
free(bsproportions);
free(bsdata);
free(PSresults);
free(VARps);
free(plessthanq);
free(qlessthanp);

// write .Random.seed out after use
PutRNGstate();

return Rris;
}
