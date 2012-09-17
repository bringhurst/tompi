
/* Information on mappings:
 *    ncol (number of owned columns on local processor)
 *    wcol (given local column number, which global column is it?)
 *    owned (is the given global column owned?  if so, which local # is it?)
 */

/*#define COMM_STATS*/
#ifndef NMAX
#define NMAX 100
#endif
#ifndef ISOUT
#define ISOUT 0
#endif
#define BAND NMAX

#define DEBUG 0
/* Level 0: nothing */
/* Level 1: progress reports (every DEBUG_REPORTLEVEL columns) */
#define DEBUG_REPORTLEVEL 500
/* Level 2: message-passing receive information */
/* Level 3: dependencies between columns */
/* Level 4: counts of columns */
/* (level i ==> level i-1,i-2,...,1 debuggin as well) */

#include <math.h>
#include <stdio.h>
#include "mpi.h"
#ifdef COMM_STATS
#   include "gthread.h"
#endif

#ifdef UseSsend
#   define MPI_Send MPI_Ssend
#endif

/* Inline functions */
#define ra(i,j) (*(((double *)a)+(j)*BAND+(i)))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

#ifdef COMM_STATS
int msgcnt = 0;
long unsigned msgvol = 0;
Mutex mutex;
#endif

void parmp_rsend (int dest, double *data, int size)
{
#ifdef COMM_STATS
   lock_mutex (mutex);
   msgcnt++;
   msgvol += size;
   unlock_mutex (mutex);
#endif
   MPI_Send (data, size, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
}

void parmp_rrecv (int src, double *data, int size)
{
   MPI_Status s;

   if (src < 0)
      MPI_Recv (data, size, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &s);
   else
      MPI_Recv (data, size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, &s);
}

void parmp_isend (int dest, int *data, int size)
{
   MPI_Send (data, size, MPI_INT, dest, 2, MPI_COMM_WORLD);
}

void parmp_irecv (int src, int *data, int size)
{
   MPI_Status s;

   if (src < 0)
      MPI_Recv (data, size, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &s);
   else
      MPI_Recv (data, size, MPI_INT, src, 2, MPI_COMM_WORLD, &s);
}

/* "The simplest method of mapping columns of the matrix to processors is a
 *  _wrap_mapping_" [K. Eswar, P. Sadayappan, C.-H. Huang, V. Visvanathan,
 * Supernodal Sparse Cholesky Factorization on Dist.-Memory Multiprocessors].
 * Assigns the "wrap" map for distributing columns of a matrix.  In it,
 *    map(i) = i mod p
 * where p is the number of processors, and i is from 0 to n-1.
 * p is the number of processors, currproc is the current processor number.
 * If parallelism is not yet available, simply pass 1 for currproc.
 */
void parlinear_map_wrap (int p, int n, int currproc, int *map)
{
   int i;

   for (i = 0; i < n; i++)
      map[i] = i % p;
}

/* Creates necessary information (ncol,wcol,owned) from just map. */
void parlinear_make_map (int p, int n, int currproc, int *map,
   int *ncol, int *wcol, int *owned)
{
   int i;

   *ncol = 0;
   for (i = 0; i < n; i++)
      if (map[i] == currproc)
      {
	 owned[i] = *ncol;
	 wcol[(*ncol)++] = i;
      }
      else
	 owned[i] = -1;
}

/* Alternating map.  Each processor gets a chunk of band/p. */
void parlinear_map_alt (int p, int n, int currproc, int *map, int band)
{
   int i;

   if ((band % p) != 0)
   {
      fprintf (stderr, "parlinear_map_alt: Invalid input.\n");
      return;
   }
   for (i = 0; i < n; i++)
      map[i] = (i % band)/(band/p);
}

/* Parallel fan-out Cholesky factorization technique for banded matrices
 * (general band).  Requires map, ncol, wcol, owned.  The matrix a
 * (n x n) must already be distributed.  band is the upper bandwidth,
 * including the diagonal.  wcol must be sorted.  a(band,*) must be the
 * diagonal.
 * ptemp is a temporary of size max_number_of_processors.
 * coltmp is a temporary of size band+1 (the first one is used for col. #).
 * colcnt is a temporary of size ncol (counts the updates).
 * Doesn't work for diagonal matrix.  Replaces a with L.
 */
void parlinear_llt_fanout (int *map, int ncol, int *wcol, int *owned,
   double **a, int band, int n, int *ptemp, double *coltmp,
   double *coltmp2, int *colcnt)
{
   int p, np, count;
   int i, i1, j, j1, jl, k, keffect, kp1;
   double temp1;

   MPI_Comm_rank (MPI_COMM_WORLD, &p);
   MPI_Comm_size (MPI_COMM_WORLD, &np);
/* Initialization - count the number of columns each owned column will be
 *    updated by.
 */
   for (i = 0; i < ncol; i++)
   {
      j = wcol[i];
      colcnt[i] = j+1 - max (j-band+2, 1);
/*    printf ("local %4d global %4d count %4d\n", i, j, colcnt[i]); */
   }
/* Step 1 */
   count = ncol;
/* Step 2 */
   jl = owned[0];
   if (jl != -1)
   {
      if (jl != 0)
      {
	 fprintf (stderr, "parlinear_llt_fanout: Col 0 must be local col 0\n");
	 return;
      }
/* Step 3 */
      if (ra(band-1,jl) < 0)
      {
	 fprintf (stderr, "parlinear_llt_fanout: Matrix not positive-def!\n");
	 fprintf (stderr, "   (at column 0)\n");
	 return;
      }
      ra(band-1,jl) = sqrt (ra(band-1,jl));
      for (i = 0; i < (band-1); i++)
	 ra(i,jl) = ra(i,jl) / ra(band-1,jl);
/* Step 4 */
      for (i = 0; i < np; i++)
         ptemp[i] = 0;
      for (j = 0; j < band; j++)
	 coltmp[j+1] = ra(j,0);
      coltmp[0] = 0;
      for (i = 1; i < min (band, n); i++)
	 if (ptemp[map[i]] == 0)
	 {
	    parmp_rsend (map[i], coltmp, band+1);
	    ptemp[map[i]] = 1;
	 }
/* Step 5 */
      count--;
   }
/* Step 6 */
   while (count > 0)
   {
/* Step 7 */
#if DEBUG >= 2
      printf ("[p%d] attempting to receive a message (cnt=%d)\n", p, count);
#endif
      parmp_rrecv (-1, coltmp, band+1);
#if DEBUG >= 2
      printf ("[p%d] received column %d\n", p, (int) coltmp[0]);
#endif
      k = coltmp[0];
/* Step 8 */
      keffect = min (k+band-1, n-1);
#if DEBUG >= 3
      printf ("[p%d] received col (%d) effects up to column %d\n", p,k,keffect);
#endif
      kp1 = k+1;
      for (j = k+1; j <= keffect; j++)
      {
	 jl = owned[j];
	 if (jl != -1)
	 {
/* Step 9 */
	    j1 = band-j+k;
/*	    if (j1 >= 0 && j1 < band) { */
	       temp1 = coltmp[j1+1];
	       for (i = j-k; i < band; i++)
		  ra(i,jl) = ra(i,jl) - coltmp[i-j+kp1]*temp1;
/*	    }
	    else
	       printf ("THIS SHOULD *NEVER* HAPPEN j=%d k=%d j1=%d\n", j,k,j1);
 */
	    colcnt[jl]--;
/* Step 10 */
#if DEBUG >= 4
            printf ("[p%d] effected column %d has count %d\n", p,j,colcnt[jl]);
#endif
	    if (colcnt[jl] == 0)
	    {
/* Step 11 */
	       if (ra(band-1,jl) < 0)
	       {
		  fprintf (stderr, "parlinear_llt_fanout: Matrix not PD!\n");
		  fprintf (stderr, "   (at column %d, val %f)\n", j,
		     (float) ra(band-1,jl));
		  return;
	       }
	       ra(band-1,jl) = sqrt(ra(band-1,jl));
	       for (i = 0; i < (band-1); i++)
		  ra(i,jl) = ra(i,jl) / ra(band-1,jl);
/* Step 12 */
               for (i = 0; i < np; i++)
	          ptemp[i] = 0;
#if DEBUG >= 3
	       printf ("[p%d] finished column %d effects up to %d\n", p,j,
	          min(j+band-1,n-1));
#endif
	       for (i = j+1; i <= min(j+band-1,n-1); i++)
		  if (ptemp[map[i]] == 0)
		  {
		     for (i1 = 0; i1 < band; i1++)
		        coltmp2[i1+1] = ra(i1,jl);
		     coltmp2[0] = j;
		     parmp_rsend (map[i], coltmp2, band+1);
		     ptemp[map[i]] = 1;
		  }
/* Step 13 */
	       count--;
#if DEBUG > 0
	       if (j % DEBUG_REPORTLEVEL == 0 || j > (n-2))
	          printf ("[debug] column %d done\n", j);
#endif
	    }
	 }
      }
   }
}

/* Parallel fan-in Cholesky factorization technique for banded matrices
 * (general band).  Requires map, ncol, wcol, owned.  The matrix a
 * (n x n) must already be distributed.  band is the upper bandwidth,
 * including the diagonal.  wcol must be sorted.  a(band,*) must be the
 * diagonal.
 * ptemp is a temporary of size max_number_of_processors.
 * coltmp is a temporary of size band+1 (the first one is used for col. #).
 * colcnt is a temporary of size ncol (counts the updates).
 * Doesn't work for diagonal matrix.  Replaces a with L.
 */
void parlinear_llt_fanin (int *map, int ncol, int *wcol, int *owned,
   double **a, int band, int n, int *ptemp, double *coltmp,
   int *colcnt)
{
   int p, np;
   int i, j, jl, jp, jpl, k, kl, flag;
   int bandmj,bandp1,bandm1,kmj,jm1;
   double temp1;
   MPI_Request req;

   MPI_Comm_rank (MPI_COMM_WORLD, &p);
   MPI_Comm_size (MPI_COMM_WORLD, &np);
   req    = MPI_REQUEST_NULL;
   bandp1 = band + 1;
   bandm1 = band - 1;
/* Step 1 */
   for (i = 0; i < ncol; i++)
   {
      j = wcol[i];
      colcnt[i] = 0;
      for (k = 0; k < np; k++)
         ptemp[k] = 0;
      ptemp[p] = 1;
      for (k = max(j-bandm1, 0); k < j; k++)
	 if (ptemp[map[k]] == 0)
	 {
	    ptemp[map[k]] = 1;
	    colcnt[i]++;
	 }
   }
/* Step 2 */
   for (j = 0; j < n; j++)
   {
/* Step 2.5 */
      bandmj = band-j;
      jm1 = j-1;
      jl = owned[j];
      if (jl != -1)
      {
/* Step 3 */
	 for (k = max(1-bandmj, 0); k < j; k++)
	 {
	    kl = owned[k];
	    if (kl != -1)
	    {
/* Step 4 */
	       temp1 = ra(bandmj+k,kl);
	       kmj = k-j;
	       for (i = 1-kmj; i < band; i++)
		  ra(i,jl) = ra(i,jl) - ra(i+kmj,kl)*temp1;
	    }
	 }
/* Step 5 */
#if DEBUG >= 2
         printf ("[p%d] Collecting data for column %d\n", p, j);
#endif
         while (colcnt[jl] > 0)
	 {
/* Step 6 */
	    parmp_rrecv (-1, coltmp, bandp1);
	    jp  = coltmp[band];
#if DEBUG >= 2
	    printf ("[p%d] Received data for column %d\n", p, jp);
#endif
	    jpl = owned[jp];
	    for (i = 0; i < band; i++)
	       ra(i,jpl) -= coltmp[i];
	    colcnt[jpl]--;
	 }
#if DEBUG >= 2
         printf ("[p%d] Column %d complete\n", p, j);
#endif
/* Step 8 */
	 if (ra(bandm1,jl) < 0)
	 {
	    fprintf (stderr, "parlinear_llt_fanin: Matrix not pos-def!\n");
	    fprintf (stderr, "   (at col %d)\n", j);
	    return;
	 }
	 ra(bandm1,jl) = sqrt (ra(bandm1,jl));
	 for (i = 0; i < bandm1; i++)
	    ra(i,jl) = ra(i,jl) / ra(bandm1,jl);
#if DEBUG > 0
	 if (j % DEBUG_REPORTLEVEL == 0 || j > (n-20))
	    printf ("[debug] column %d done\n", j);
#endif
/* Step 9 */
      } else {
/* Step 10 */
         for (i = 0; i < band; i++)
	    coltmp[i] = 0;
	 flag = 0;
/* Non-blocking */
/*       if (req != MPI_REQEUST_NULL) */
/*          parmp_wait (req); */
/* Step 11 */
         for (k = max (1-bandmj, 0); k <= jm1; k++)
	 {
	    kl = owned[k];
	    if (kl != -1)
	    {
#if DEBUG >= 3
	       printf ("[p%d] Owned column %d effects column %d\n", p, k, j);
#endif
/* Step 12 */
	       temp1 = ra(bandmj+k,kl);
	       kmj = k-j;
	       for (i = 1-kmj; i < band; i++)
		  coltmp[i] += ra(i+kmj,kl)*temp1;
	       flag = 1;
	    }
	 }
/* Step 13 */
	 if (flag != 0)
	 {
	    coltmp[band] = j;
/* Blocking */
	    parmp_rsend (map[j], coltmp, bandp1);
/* Non-blocking */
/*          parmp_rsend_nb (map[j], coltmp, bandp1, req); */
	 }
      }
   }
}

int main (int argc, char **argv)
{
   int ncolmax;
   int ncol;
   int *map, *wcol, *owned, *colcnt, *ptemp;
   double **a, *coltmp, *coltmp2;
   int p, np, i, j;
   double start, end, diff, maxdiff;

#ifdef initial_code
   initial_code
#endif

#ifdef MPII_STACK_SIZE
   MPII_Stack_size = 10000;
#endif
   MPI_Init (&argc, &argv);

   MPI_Comm_rank (MPI_COMM_WORLD, &p);
   MPI_Comm_size (MPI_COMM_WORLD, &np);

   ncolmax = NMAX/np+1;

   a = (double **) malloc (sizeof (double) * ncolmax * BAND);
   map = (int *) malloc (sizeof (int) * NMAX);
   wcol = (int *) malloc (sizeof (int) * (NMAX/np + 1));
   owned = (int *) malloc (sizeof (int) * NMAX);
   colcnt = (int *) malloc (sizeof (int) * (NMAX/np + 1));
   ptemp = (int *) malloc (sizeof (int) * np);
   coltmp = (double *) malloc (sizeof (double) * (BAND+1));
   coltmp2 = (double *) malloc (sizeof (double) * (BAND+1));
   if (a == NULL || map == NULL || wcol == NULL || owned == NULL ||
           colcnt == NULL || ptemp == NULL || coltmp == 0 || coltmp2 == NULL)
   {
      fprintf (stderr, "%s: Out of memory!\n", argv[0]);
      exit (1);
   }

   for (i = 0; i < ncolmax; i++)
   {
      ra(BAND-1,i) = 1;
      for (j = 0; j < (BAND-1); j++)
         ra(j,i) = 0;
   }

   parlinear_map_wrap (np, NMAX, p, map);
   parlinear_make_map (np, NMAX, p, map, &ncol, wcol, owned);

   /* Start timing */
#ifdef PUPPET
   MPI_Pcontrol (3);
#else
   MPI_Barrier (MPI_COMM_WORLD);
   start = MPI_Wtime ();
#endif

#if ISOUT == 0
   parlinear_llt_fanin (map, ncol, wcol, owned, a, BAND, NMAX, ptemp, coltmp,
           colcnt);
#else
   parlinear_llt_fanout (map, ncol, wcol, owned, a, BAND, NMAX, ptemp, coltmp,
           coltmp2, colcnt);
#endif

#ifndef PUPPET
   end = MPI_Wtime ();
   diff = end - start;
   sleep (2);
   if (p == 0)
   {
      int i;
      MPI_Status s;
      maxdiff = diff;
      for (i = 1; i < np; i++)
      {
         MPI_Recv (&diff, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 7, MPI_COMM_WORLD, &s);
         if (diff > maxdiff)
            maxdiff = diff;
      }
      printf ("faninout: Total execution time is %f\n", maxdiff);
   }
   else
      MPI_Send (&diff, 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
#endif

   MPI_Finalize ();

#ifdef COMM_STATS
   if (p == 0)
   {
      printf ("This is the communication-stats version, so ignore the time above\n");
      printf ("Message count: %d      Message volume: %d*%d\n", msgcnt, sizeof (double), msgvol);
   }
#endif

   return 0;
}

