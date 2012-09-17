#include "mpi.h"

int main (int argc, char *argv[])
{
   int rank, size;
   MPI_Request req;
   MPI_Status status;

#ifdef initial_code
   initial_code
#endif

#ifdef MPII_STACK_SIZE
   MPII_Stack_size = 10000;
#endif

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   printf ("Process %d sleeping...\n", rank);
   sleep (1);
   printf ("Process %d done sleeping.\n", rank);
   MPI_Finalize ();

   return 0;
}

