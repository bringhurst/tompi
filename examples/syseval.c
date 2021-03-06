/* System evaluator.
 * Determines network performance according to the model
 *   time = node distance * (init. time * ceil(size/blocking factor)
 *                           + size * byte time)
 * where init. time, blocking factor, and byte time are the parameters
 * generated by the system evaluator.
 * All tests are done node-to-node.  (This program should be run with
 * two processors.)
 */

#include <stdio.h>
#include "mpi.h"

#ifdef UseSsend
#ifndef MPI_Send
#define MPI_Send MPI_Ssend
#endif
#endif

#ifndef WARMUP
#define WARMUP 100
#endif

#define NTESTS 7
#define MAXSIZE 1000000  /* Biggest end */
#define MAXLIST 100      /* Biggest it_median */
struct
{
   int begin;
   int end;
   int step;
   int it_average;
   int it_median;
} tests[NTESTS] = {
#if defined(test_values)
   test_values
#else
   {0, 2, 1, 1000, 11},
   {10, 11, 1, 1000, 11},
   {100, 101, 1, 1000, 11},
   {1000, 1001, 1, 500, 11},
   {10000, 10001, 1, 100, 11},
   {100000, 100001, 1, 50, 11},
   {1000000, 1000001, 1, 50, 11}
#endif
};
/* From PUPPET:
   {   0,  2048,  128,  20, 61},
   {2048,  8192,  256,  20, 51},
   {8192, 32768,  512,  20, 31}
 */

#define p0 if (proc == 0)
#define p1 if (proc == 1)
#define sync MPI_Barrier (MPI_COMM_WORLD)

typedef struct
{
   int size;
   double time;
} datum;

double median (double *a, int n);

main (int argc, char *argv[])
{
   int test, size, totalsize, i, j, it_avg, it_med;
   datum *data;
   int nproc, proc;
   char *buf;
   MPI_Status status;
   double begintime, endtime;
   double list[MAXLIST];
   int err;
   FILE *f;
   char statsfile[256];

#ifdef initial_code
   initial_code
#endif

   buf = (char *) malloc (MAXSIZE);
   if (buf == NULL)
   {
      perror ("Couldn't allocate data");
      exit (1);
   }

   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &nproc);
   if (nproc < 2)
   {
      fprintf (stderr, "syseval: Two processors required.\n");
      MPI_Finalize ();
      exit (1);
   }
   MPI_Comm_rank (MPI_COMM_WORLD, &proc);

   p0
   {
      for (totalsize = 0, test = 0; test < NTESTS; test++)
	 for (size = tests[test].begin; size < tests[test].end;
	      size += tests[test].step)
	    totalsize++;
      data = (datum *) malloc (sizeof (datum) * totalsize);
      if (data == NULL)
      {
         fprintf (stderr, "syseval: Out of memory!\n");
         err = 1;
      }
      else
         err = 0;
   }
   MPI_Bcast (&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (err)
   {
      MPI_Finalize ();
      exit (1);
   }

   /* Warm up */
   p0
   {
      printf ("Warming up...\n");
      fflush (stdout);
   }
   p0 for (i = 0; i < WARMUP; i++)
   {
      MPI_Send (buf, MAXSIZE, MPI_BYTE, 1, 100, MPI_COMM_WORLD);
      MPI_Recv (buf, MAXSIZE, MPI_BYTE, 1, 101, MPI_COMM_WORLD, &status);
   }
   p1 for (i = 0; i < WARMUP; i++)
   {
      MPI_Recv (buf, MAXSIZE, MPI_BYTE, 0, 100, MPI_COMM_WORLD, &status);
      MPI_Send (buf, MAXSIZE, MPI_BYTE, 0, 101, MPI_COMM_WORLD);
   }

   /* Collect statistics */
   for (totalsize = 0, test = 0; test < NTESTS; test++)
   {
      p0
      {
	 printf ("Starting test %d/%d (size %d-%d).\n", test+1, NTESTS,
	    tests[test].begin, tests[test].end);
	 fflush (stdout);
      }
      /* sync; */ /* originally in j loop */
      for (size = tests[test].begin; size < tests[test].end;
           size += tests[test].step)
      {
	 it_avg = tests[test].it_average;
	 it_med = tests[test].it_median;
	 p0
	 {
	    for (j = 0; j < it_med; j++)
	    {
	       begintime = MPI_Wtime ();
	       for (i = 0; i < it_avg; i++)
	       {
		  MPI_Send (buf, size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
		  MPI_Recv (buf, size, MPI_BYTE, 1, 1, MPI_COMM_WORLD,
		     &status);
	       }
	       endtime = MPI_Wtime ();
	       list[j] = (endtime - begintime) / 2.0 / it_avg;
	    }
            data[totalsize].size = size;
	    data[totalsize].time = median (list, it_med);
	    totalsize++;
	 }
	 p1
	    for (j = 0; j < it_med; j++)
	       for (i = 0; i < it_avg; i++)
	       {
		  MPI_Recv (buf, size, MPI_BYTE, 0, 0, MPI_COMM_WORLD,
		     &status);
		  MPI_Send (buf, size, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
	       }
      }
   }

   p0
   {
      printf ("Writing data...\n");
      fflush (stdout);
      sprintf (statsfile, "%s.stats", argv[0]);
      f = fopen (statsfile, "w");
      if (f == NULL)
      {
         perror ("syseval: Couldn't open stats file");
	 MPI_Finalize ();
	 exit (1);
      }
      fprintf (f, "Size (byt) Time (sec)\n");
      for (test = 0; test < totalsize; test++)
	 fprintf (f, "%10d %e\n", data[test].size, data[test].time);
      fclose (f);
   }

   MPI_Finalize ();

   return 0;
}

int compdbl (double *x, double *y)
{
   return *x - *y;
}
double median (double *a, int n)
{
   qsort (a, n, sizeof (double), compdbl);
   return a[n/2 + 1];
}

