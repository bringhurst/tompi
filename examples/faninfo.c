
/* Information on mappings:
 *    ncol (number of owned columns on local processor)
 *    wcol (given local column number, which global column is it?)
 *    owned (is the given global column owned?  if so, which local # is it?)
 */

#ifndef NMAX
#define NMAX 1000
#endif
#ifndef ISOUT
#define ISOUT 0
#endif
#ifndef NPROC
#define NPROC 1         /* seems broken otherwise */
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
#include "gthread.h"

/* Inline functions */
#define ra(i,j) (*(((double *)a)+(j)*band+(i)))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

int msgcnt = 0, msgvol = 0, flops = 0;
Mutex mutex;

void parmp_rsend (int dest, double *data, int size)
{
   lock_mutex (mutex);
   msgcnt++;
   msgvol += sizeof(double)*size;
   unlock_mutex (mutex);
   MPI_Send (data, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
}

void parmp_rrecv (int src, double *data, int size)
{
   MPI_Status s;

   if (src < 0)
      MPI_Recv (data, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &s);
   else
      MPI_Recv (data, 1, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, &s);
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
      lock_mutex (mutex);
      flops++;
      for (i = 0; i < (band-1); i++)
         flops++;
      unlock_mutex (mutex);
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
               lock_mutex (mutex);
	       for (i = j-k; i < band; i++)
                  flops += 2;
               unlock_mutex (mutex);
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
               lock_mutex (mutex);
               flops++;
	       for (i = 0; i < (band-1); i++)
                  flops++;
               unlock_mutex (mutex);
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
               lock_mutex (mutex);
	       for (i = 1-kmj; i < band; i++)
                  flops += 2;
               unlock_mutex (mutex);
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
            lock_mutex (mutex);
	    for (i = 0; i < band; i++)
               flops++;
            unlock_mutex (mutex);
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
         lock_mutex (mutex);
         flops++;
	 for (i = 0; i < bandm1; i++)
            flops++;
         unlock_mutex (mutex);
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
               lock_mutex (mutex);
	       for (i = 1-kmj; i < band; i++)
                  flops += 2;
               unlock_mutex (mutex);
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
   const int n = NMAX, minp = NPROC, maxp = NPROC, band = BAND;
   const int ncolmax = n/minp+1;
   int map[NMAX], ncol, wcol[NMAX/NPROC+1], owned[NMAX], colcnt[NMAX/NPROC+1];
   int ptemp[NPROC];
   double **a, coltmp[BAND+1], coltmp2[BAND+1];
   int p, np, i, j;
   double start, end, diff, maxdiff;

   MPI_Init (&argc, &argv);

   a = (double **) malloc (sizeof (double) * ncolmax * band);
   if (a == NULL)
   {
      fprintf (stderr, "%s: Out of memory!\n", argv[0]);
      exit (1);
   }

   for (i = 0; i < ncolmax; i++)
   {
      ra(band-1,i) = 1;
      for (j = 0; j < (band-1); j++)
         ra(j,i) = 0;
   }

   MPI_Comm_rank (MPI_COMM_WORLD, &p);
   MPI_Comm_size (MPI_COMM_WORLD, &np);
   if (np != NPROC)
   {
      fprintf (stderr, "Wrong number of processes!\n");
      exit (1);
   }
   parlinear_map_wrap (np, n, p, map);
   parlinear_make_map (np, n, p, map, &ncol, wcol, owned);

   /* Start timing */
#ifdef PUPPET
   MPI_Pcontrol (3);
#else
   MPI_Barrier (MPI_COMM_WORLD);
   start = MPI_Wtime ();
#endif

#if ISOUT == 0
   parlinear_llt_fanin (map, ncol, wcol, owned, (double **) a, band, n,
      ptemp, coltmp, colcnt);
#else
   parlinear_llt_fanout (map, ncol, wcol, owned, (double **) a, band, n,
      ptemp, coltmp, coltmp2, colcnt);
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

   printf ("This is the test version, so times don't mean much.\n");
   printf ("Flops: %d    Msg cnt: %d    Msg vol (bytes): %d\n", flops,
         msgcnt, msgvol);

   return 0;
}

