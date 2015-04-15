# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <R.h>
# include <Rinternals.h>
# include <Rdefines.h>
# include <Rmath.h>
// the E measure of interindividual variation and the value of The program calculates
// the Total Niche Width (TNW), and breaks TNW down into its Between Individual Component (BIC)
// and Within Individual Component (WIC). It reports these three statistics and the
// proportion WIC/TNW, or what we call W/T. Niche widths of the population, within individuals,
// and among individuals, are calculated following Roughgarden 1974.
// A Monte Carlo resampling routine is used to calculate a null distribution for each of these
// statistics and generate a p-value testing the null hypothesis that all individuals sample
// equally from the population diet distribution.the \eqn{C_{ws}} measure of modularity
// following Araujo et al. (2008).
//
// Author: Dr. Nicola ZACCARELLI (nicola.zaccarelli@gmail.com)
//
// Version 1.3
// Date: 2015-02-03

SEXP WTcMC(SEXP comcalc, SEXP nreplicates, SEXP weight)
{

/* Declare all variables */
int   NRows, NCols, i, j, t, R, Nitems, Yitems, nreps, weightOpt;
double temp, summation, sumofsquare, wic, bic, meanmeanxi, meanxij, varxij, tnw, meanweight, numTot;
double newi, newj;
double *dati, *risult;
double **bsdata;
double *varxj, *meanxi, *numxi, *numInt;
SEXP Rris, Rdim;

   Rdim = getAttrib(comcalc, R_DimSymbol);
   NRows = INTEGER(Rdim)[0];
   NCols = INTEGER(Rdim)[1];

// Coerce data matrix into vector for data access form C
   PROTECT(comcalc = coerceVector(comcalc, REALSXP));
   dati = REAL(comcalc);

// Create the return matrix
   nreplicates = coerceVector(nreplicates, INTSXP);
   nreps = INTEGER(nreplicates)[0];
   PROTECT(Rris = allocMatrix(REALSXP, (nreps + 1), 4));
   risult = REAL(Rris);

// Weighting option
   weight = coerceVector(weight, INTSXP);
   weightOpt = INTEGER(weight)[0];

// Create the matrix for bootstrapping
   bsdata = calloc(NRows, sizeof(double *));
   for (i=0; i < NRows; i++) bsdata[i] = calloc(NCols, sizeof(double));

// NOW we calculate the WT statistics for real data
varxj  = calloc(NRows, sizeof(double));
meanxi = calloc(NRows, sizeof(double));
numxi  = calloc(NRows, sizeof(double));
numInt  = calloc(NRows, sizeof(double));
numTot = 0;
//for each individual
for (i=0; i<NRows; i++)
     {
      Nitems = 0;
      summation = 0;
//determine its mean niche
      for (j=0; j<NCols; j++)
         {
// A different way of addressing the matrix elements in the vector
         if (dati[i + NRows*j]>0) 
            {
             summation = summation + dati[i + NRows*j];
             Nitems++;
             }
         numxi[i] = Nitems;
         }
      if (Nitems == 0) continue;
      meanxi[i] = summation / (double)Nitems;
//determine its niche variance
      sumofsquare = 0;
      for (j=0; j<NCols; j++)
         {
          if (dati[i + NRows*j]>0)
           {sumofsquare = sumofsquare + ((dati[i + NRows*j] - meanxi[i]) * (dati[i + NRows*j] - meanxi[i]));}
          }
          varxj[i] = sumofsquare/(double)Nitems;
     }
for (i=0; i<NRows; i++) { numTot = numxi[i] + numTot;}
// Initialize variables
wic = 0;
bic = 0;
tnw = 0;

// This change is equivalent to Travis Ingram's solution of weighting each individual
//   regardless of diet item number

if (weightOpt == 1)
     {
// Calculate WIC and BIC when all itmes have an equal weight, TNW is calculated as sum
// Within individual component (WIC) is the average of the niche variances
      for (i=0; i<NRows; i++) { wic = wic + varxj[i];}
      wic = wic / (double)NRows;
// Between individual component (BIC) is the variance of the niche averages
      meanmeanxi = 0;
      for (i=0; i<NRows; i++) { meanmeanxi = meanmeanxi + meanxi[i];}
      meanmeanxi = meanmeanxi / (double)NRows;
      for (i=0; i<NRows; i++) {bic = bic + ((meanxi[i] - meanmeanxi)*(meanxi[i] - meanmeanxi));}
      bic = bic / (double)NRows;
      tnw = bic + wic;
     } else
     {
// Calculate TNW and BIC when itmes have a weight proportional to the number of items, WIC is calculated as difference
// Between individual component (BIC) is the variance of the niche averages
      meanweight = 0;
      for (i=0; i<NRows; i++)
         {
          numInt[i] = numxi[i] / numTot;
          meanweight = meanweight + (numInt[i] * meanxi[i]);
          }
      for (i=0; i<NRows; i++) { bic = bic + ((meanxi[i] - meanweight) *  (meanxi[i] - meanweight) * numInt[i]);}
// Calculate TNW by averaging across all individuals
      meanxij = 0;
      varxij = 0;
      Yitems = 0;
      for (i=0; i<NRows; i++)
         {
          for (j=0; j<NCols; j++)
             {
              if (dati[i + NRows*j] > 0)
                 {
                  meanxij = meanxij + dati[i + NRows*j];
                  Yitems++;
                  }
             }
         }
      meanxij = meanxij/(double)Yitems;
      for (i=0; i<NRows; i++)
         {
          for (j=0; j<NCols; j++)
             {
              if (dati[i + NRows*j] > 0)
                 {
                   varxij = varxij + ((dati[i + NRows*j] - meanxij) * (dati[i + NRows*j] - meanxij));
                  }
             }
         }
// double tnw = varxij / (double)Yitems;
    tnw = varxij / (double)Yitems;
    wic = tnw - bic;
    }

// Assign results to proper place in Rris
R = 0;
risult[R + (nreps + 1) * 0] = wic;
risult[R + (nreps + 1) * 1] = bic;
risult[R + (nreps + 1) * 2] = tnw;
risult[R + (nreps + 1) * 3] = wic/tnw;

// NOW we calculate the WT statistics for real data
// Monte Carlo Resampling

// read in (or create) .Random.seed for R random generation fuctions
GetRNGstate();

for (R = 1; R< (nreps + 1); R++)
    {
//RESAMPLING THE DATA SET
     for (i=0; i<NRows; i++)
         {
// count number of items an individual used
            Nitems = 0;
            for (j=0; j<NCols; j++)
                {if (dati[i + NRows*j]>0) Nitems++;}

// ASSIGNS RANDOMLY CHOSEN DATA POINT FROM ORIGINAL DATA TO CELL IN NEW DATA
            for (j=0; j<Nitems; j++)
                {
                temp = 0;
                while (temp == 0)
                    {
                    newi = (double)NRows * unif_rand();    // Using R random function from Rmath.h
                    newj = (double)NCols * unif_rand();    // Using R random function from Rmath.h
                    temp = dati[(int)newi + NRows*(int)newj];
                    }
                bsdata[i][j] = temp;
                }
            for (j=Nitems; j<NCols; j++)
                {bsdata[i][j] = 0;}
            }
     for (t=0; t<NRows; t++)
         {
         varxj[t]  = 0;
         meanxi[t] = 0;
         }
for (i=0; i<NRows; i++)
     {
      Nitems = 0;
      summation = 0;
//determine its mean niche
      for (j=0; j<NCols; j++)
         {
// A different way of addressing the matrix elements in the vector
         if (bsdata[i][j]>0) 
            {
             summation = summation + bsdata[i][j];
             Nitems++;
             }
// numxi[i] = Nitems;
         }
      if (Nitems == 0) continue;
      meanxi[i] = summation / (double)Nitems;
//determine its niche variance
      sumofsquare = 0;
      for (j=0; j<NCols; j++)
         {
          if (bsdata[i][j]>0)
           {sumofsquare = sumofsquare + ((bsdata[i][j] - meanxi[i]) * (bsdata[i][j] - meanxi[i]));}
          }
          varxj[i] = sumofsquare/(double)Nitems;
     }
wic = 0;
bic = 0;
tnw = 0;

// This change is equivalent to Travis Ingram's solution of weighting each individual
//   regardless of diet item number */

if (weightOpt == 1)
     {
// Calculate WIC and BIC when all itmes have an equal weight, TNW is calculated as sum
// Within individual component (WIC) is the average of the niche variances
      for (i=0; i<NRows; i++) { wic = wic + varxj[i];}
      wic = wic / (double)NRows;
// Between individual component (BIC) is the variance of the niche averages
      meanmeanxi = 0;
      for (i=0; i<NRows; i++) { meanmeanxi = meanmeanxi + meanxi[i];}
      meanmeanxi = meanmeanxi / (double)NRows;
      for (i=0; i<NRows; i++) {bic = bic + ((meanxi[i] - meanmeanxi)*(meanxi[i] - meanmeanxi));}
      bic = bic / (double)NRows;
      tnw = bic + wic;
     } else
     {
// Calculate NW and BIC when itmes have a weight proportional to the number of items, WIC is calculated as difference
// Between individual component (BIC) is the variance of the niche averages
      meanweight = 0;
      for (i=0; i<NRows; i++)
         {
          numInt[i] = numxi[i] / numTot;
          meanweight = meanweight + (numInt[i] * meanxi[i]);
          }
      for (i=0; i<NRows; i++) { bic = bic + ((meanxi[i] - meanweight) *  (meanxi[i] - meanweight) * numInt[i]);}
//calculate TNW by averaging across all individuals
      meanxij = 0;
      varxij = 0;
      Yitems = 0;
      for (i=0; i<NRows; i++)
         {
          for (j=0; j<NCols; j++)
             {
              if (dati[i + NRows*j] > 0)
                 {
                  meanxij = meanxij + dati[i + NRows*j];
                  Yitems++;
                  }
             }
         }
      meanxij = meanxij/(double)Yitems;
      for (i=0; i<NRows; i++)
         {
          for (j=0; j<NCols; j++)
             {
              if (dati[i + NRows*j] > 0)
                 {
                   varxij = varxij + ((dati[i + NRows*j] - meanxij) * (dati[i + NRows*j] - meanxij));
                  }
             }
         }
    tnw = varxij / (double)Yitems;
    wic = tnw - bic;
    }

  risult[R + (nreps + 1) * 0] = wic;
  risult[R + (nreps + 1) * 1] = bic;
  risult[R + (nreps + 1) * 2] = tnw;
  risult[R + (nreps + 1) * 3] = wic / tnw;

} // END Monte Carlo resampling
UNPROTECT(2);
free(bsdata);
free(varxj);
free(meanxi);
free(numxi);
free(numInt);

// write .Random.seed out after use
PutRNGstate();

return Rris;
}

