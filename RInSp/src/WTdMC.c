# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <R.h>
# include <Rinternals.h>
# include <Rdefines.h>
# include <Rmath.h>
// The procedure performs a Monte Carlo resampling under a null hypothesis to calculate
// The program calculates the Total Niche Width (TNW), and breaks TNW down into its
// Between Individual Component (BIC) and Within Individual Component (WIC). It reports
// these three statistics and the proportion WIC/TNW, or what we call W/T. Niche widths
// of the population, within individuals, and among individuals, are calculated using
// the Shannon-Weaver diversity index on count data, following Roughgarden (1979) or
// Bolnick et al (2002).
// A Monte Carlo resampling routine is used to calculate a null distribution for each
// of these statistics and generate a p-value testing the null hypothesis that all individuals
// sample equally from the population diet distribution. Note that the bootstrapping procedure
// is only biologically meaningful when applied to integer data representing counts of
// individual prey within individual predators' diets, and should not be applied to other forms
// of data such as prey mass or proportion of total volume.
//
// Author: Dr. Nicola ZACCARELLI (nicola.zaccarelli@gmail.com)
//
// Version 1.3
// Date: 2015-02-03

SEXP WTdMC(SEXP comcalc, SEXP popdiet, SEXP nreplicates)
{
// Declare all variables
int   NRows, NCols, i, j, nreps, popD, zeroes, x, R;
double popshannonweaver, wpc, indsw, bspoptotal;
double cumulativep, lowerbound, bsresourcejtotal;
double item;
double *dati, *results, *totaldieti, *bstotaldieti, *populationdiet;
double *bspopulationdiet, *resfin;
double **bsdata, **bsproportions;

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
PROTECT(Rris = allocMatrix(REALSXP, (nreps + 1), 5));
resfin = REAL(Rris);
bsdata = calloc(NRows, sizeof(double *));
for (i=0; i < NRows; i++) bsdata[i] = calloc(NCols, sizeof(double));

// Creation of necessary matrices in C
bsproportions = calloc(NRows, sizeof(double *));
for (i=0; i < NRows; i++) bsproportions[i] = calloc(NCols, sizeof(double));
bstotaldieti = calloc(NRows, sizeof(double));
bspopulationdiet = calloc(NCols, sizeof(double));
results = calloc(NRows, sizeof(double));
totaldieti = calloc(NRows, sizeof(double));
populationdiet = calloc(NCols, sizeof(double));
// Maybe this cycle is useless
 for (j=0; j < NCols; j++)
    {populationdiet[j] = 0;
     bspopulationdiet[j] = 0;}

// read in (or create) .Random.seed for R random generation fuctions
GetRNGstate();

// Attention to the right address for the elements
for (j=0; j<NCols; j++) {
  for (i=0; i<NRows; i++){
  bsdata[i][j] = (double)dati[i + NRows*j];
  bsproportions[i][j] = 0;
}}

// Calculate individual diet and population diet
for (i=0; i < NRows; i++)
     {totaldieti[i] = 0;
      bstotaldieti[i] = 0;}
bspoptotal = 0;
for (i=0; i < NRows; i++) {
    for (j=0; j < NCols; j++) {
         totaldieti[i] = totaldieti[i] + bsdata[i][j];}
	bspoptotal = bspoptotal + totaldieti[i];}

// Start with calculate real data and than to simulate
for (R = 0; R < (nreps + 1); R++)
{
if (R > 0) {
	 for (i=0; i<NRows; i++)
         {for (j=0; j<NCols; j++)
               {bsdata[i][j] = 0;}}
     for (i=0; i<NRows; i++)
         {for (x=0; x< totaldieti[i]; x++)
             {item = unif_rand();    // Using R random function from Rmath.h
              cumulativep = 0;
              for (j=0; j<NCols; j++)
                  {lowerbound = cumulativep;
                   cumulativep = cumulativep + populationdiet[j];
                   if (item>=lowerbound && item<cumulativep)
                  {bsdata[i][j]=bsdata[i][j] + 1;}}
          } // for x
      }; // for i
} // end if (R>0)


// As the marginal total for rows are fixed the proportions calculation is easier
for (i=0; i < NRows; i++)
   { for (j=0; j < NCols; j++){
          bsproportions[i][j] = (bsdata[i][j] / totaldieti[i]);}
   }

// CALCULATE POPULATION DIET
switch (popD)
  {case 0 :
    {
     for (j = 0; j<NCols; j++) 
         {bsresourcejtotal = 0;
         for (i=0; i<NRows; i++) {bsresourcejtotal = bsresourcejtotal + bsdata[i][j];}
         bspopulationdiet[j] = ( bsresourcejtotal/bspoptotal );
		 if (R == 0) { populationdiet[j] = bspopulationdiet[j];}}
     break;}
   case 1:
    {for (j=0; j<NCols; j++)
       {bsresourcejtotal = 0;
        for (i=0; i<NRows; i++) {bsresourcejtotal = bsresourcejtotal + bsproportions[i][j];}
        bspopulationdiet[j] = (bsresourcejtotal / NRows);
		if (R == 0) { populationdiet[j] = bspopulationdiet[j];}
		}
     break;}
    }

popshannonweaver = 0;
zeroes = 0;
wpc = 0;

for (i=0; i< NRows; i++) {
    results[i] = 0;
    }

// CALCULATE TNW statistics
for (j=0; j<NCols; j++)
  {if (bspopulationdiet[j] > 0)
      {popshannonweaver = popshannonweaver - (bspopulationdiet[j]*log(bspopulationdiet[j]));}}

// CALCULATE WIC
 for (i=0; i<NRows; i++)
    {indsw = 0;
     for (j=0; j<NCols; j++)
        {if (bsproportions[i][j]>0)
            {indsw = indsw - (bsproportions[i][j]*log(bsproportions[i][j]));}}
     results[i] = indsw;
     if (results[i] == 0) {zeroes++;}
     }

for (i=0; i<NRows; i++)
    {
     wpc = wpc + (totaldieti[i]/ bspoptotal)*results[i];}

resfin[R + (nreps + 1) * 0] = zeroes;
resfin[R + (nreps + 1) * 1] = wpc;
resfin[R + (nreps + 1) * 2] = popshannonweaver - wpc;
resfin[R + (nreps + 1) * 3] = popshannonweaver;
resfin[R + (nreps + 1) * 4] = wpc/popshannonweaver;
}
UNPROTECT(2);
free(bsdata);
free(bsproportions);
free(bspopulationdiet);
free(results);
free(totaldieti);
free(bstotaldieti);
free(populationdiet);

// write .Random.seed out after use
PutRNGstate();

return Rris;
}
