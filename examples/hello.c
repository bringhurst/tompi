#include "mpi.h"

int main (int argc, char *argv[])
{
   int rank, size;
   MPI_Request req;
   MPI_Status status;

#ifdef initial_code
   initial_code
#endif

   MPII_Stack_size = 10000;

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   printf ("Hello world from process %d/%d\n", rank+1, size);
   if (rank == 0)
   {
      int i;
      MPI_Request *req =
         (MPI_Request *) malloc (sizeof (MPI_Request) * (size-1));
      MPI_Status *statuses =
         (MPI_Status *) malloc (sizeof (MPI_Status) * (size-1));
      for (i = 1; i < size; i++)
         MPI_Issend (&rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &(req[i-1]));
      MPI_Waitall (size-1, req, statuses);
   }
   else
   {
      MPI_Recv (&rank, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
      printf ("Process %d (source %d, tag %d) also says hello\n", rank+1,
            status.MPI_SOURCE, status.MPI_TAG);
   }
   MPI_Finalize ();

   return 0;
}

