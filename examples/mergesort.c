/* Parallel merge sort */

#define DEBUG 0
/* Levels:
 *   0: no debugging information
 *   1: small progress report
 *   2: report all message-passing
 */

#include <stdio.h>
#include "mpi.h"

/* Constants */
#ifndef NPROC
#   define QSORT 0                      /* Set to 1 for quicksort  */
#   define NPROC 4                      /* Number of processors    */
#   define TOTALSIZE 5000               /* Size of total array     */
#endif
#define SIZE (TOTALSIZE/NPROC)		/* Size of local array     */
#define APPROXTOTAL (SIZE*NPROC)	/* Accounts for round off  */

static int intcompare (int *i, int *j)  /* For qsort */
{
   return *i - *j;
}

int main (int argc, char *argv[])
{
   int a[APPROXTOTAL], b[APPROXTOTAL], c[APPROXTOTAL], swaptemp;
   int i, j, k, proc, realproc, mult;
   MPI_Status status;
   double start, end, diff, maxdiff;

   /* Initialize the array to sort with random values */
   for (i = 0; i < SIZE; i++)
      a[i] = rand ();

   /* Initialize MPI and start timing */
   MPI_Init (&argc, &argv);

   {
      int np;
      MPI_Comm_size (MPI_COMM_WORLD, &np);
      if (np != NPROC)
      {
         printf ("mergesort: Wrong number of processors (change NPROC)\n");
         exit (1);
      }
   }

#ifndef PUPPET
   MPI_Barrier (MPI_COMM_WORLD);
   start = MPI_Wtime ();
#endif

#if QSORT == 0
   /* Local bubble sort (to use a significant amount of time) */
   for (i = 0; i < SIZE; i++)
      for (j = i + 1; j < SIZE; j++)
         if (a[i] > a[j])
	 {
	    swaptemp = a[i];
	    a[i]     = a[j];
	    a[j]     = swaptemp;
	 }
#else
   /* Local quick sort */
   qsort (a, SIZE, sizeof (int), intcompare);
#endif

   MPI_Comm_rank (MPI_COMM_WORLD, &realproc); /* Get local processor number */
#if DEBUG >= 1
   printf ("[%d] bubble sorted %d elements\n", realproc, SIZE);
#endif

   for (proc = realproc, mult = 1; mult != NPROC; proc /= 2, mult *= 2)
      if (proc % 2 == 0)            /* Collect and merge */
      {
#if DEBUG >= 2
         printf ("[%d] recv.. pend\n", realproc);
#endif
         MPI_Recv (c, SIZE * mult, MPI_INT, MPI_ANY_SOURCE, mult,
	    MPI_COMM_WORLD, &status);
#if DEBUG >= 2
         printf ("[%d] recv.. got it\n", realproc);
#endif
	 for (i = 0; i < (SIZE * mult); i++)
	    b[i] = a[i];            /* Save old array */
	 for (i = 0, j = 0, k = 0; (i < (SIZE * mult)) && (j < (SIZE * mult));
	      k++)
	    if (b[i] < c[j])
	       a[k] = b[i++];
	    else
	       a[k] = c[j++];
         if (i != (SIZE * mult))
	    for (; i < (SIZE * mult); i++, k++)
	       a[k] = b[i];
	 else
	    for (; j < (SIZE * mult); j++, k++)
	       a[k] = c[j];
      } else {                      /* Send */
#if DEBUG >= 2
         printf ("[%d] sending to %d\n", realproc, (realproc/(mult*2))*mult*2);
#endif
         MPI_Send (a, SIZE * mult, MPI_INT, (realproc / (mult*2)) * mult*2,
	    mult, MPI_COMM_WORLD);
         break;
      }

#if DEBUG >= 1
   printf ("[%d] finalizing\n", realproc);
#endif
#ifndef PUPPET
   end = MPI_Wtime ();
   diff = end - start;
   sleep (2);
   if (realproc == 0)
   {
      int i;
      MPI_Status s;
      for (i = 1; i < NPROC; i++)
         MPI_Send (&i, 1, MPI_INT, i, 666, MPI_COMM_WORLD);
      maxdiff = diff;
      for (i = 1; i < NPROC; i++)
      {
         MPI_Recv (&diff, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 667, MPI_COMM_WORLD,
               &s);
         if (diff > maxdiff)
            maxdiff = diff;
      }
      printf ("mergesort: Total execution time is %f\n", maxdiff);
   }
   else
   {
      int i;
      MPI_Status s;
      MPI_Recv (&i, 1, MPI_INT, 0, 666, MPI_COMM_WORLD, &s);
      MPI_Send (&diff, 1, MPI_DOUBLE, 0, 667, MPI_COMM_WORLD);
   }
   end = MPI_Wtime ();
#endif

   MPI_Finalize ();                 /* Deinitialize MPI and stop timing */

   if (realproc == 0)               /* Ensure the array is sorted */
   {
      for (i = 0; i < (APPROXTOTAL-1); i++)
         if (a[i] > a[i+1])
	 {
	    printf ("mergesort: Test failed (unsorted)\n");
	    break;
	 }
      if (i == (APPROXTOTAL-1))
	 printf ("mergesort: Test succeeded!\n");
   }

   return 0;
}

/* EOF */
