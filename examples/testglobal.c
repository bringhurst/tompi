#include <mpi.h>

int p, np;

void f ()
{
    static int i = 0;
    i++;
    printf ("[%d/%d] f: This should be one: %d\n", p+1, np, i);
}

void main (int argc, char *argv[])
{
    MPI_Init (&argc, &argv);
    printf ("1\n");
    MPI_Comm_rank (MPI_COMM_WORLD, &p);
    printf ("2\n");
    MPI_Comm_size (MPI_COMM_WORLD, &np);
    printf ("3\n");
    printf ("[%d/%d] Calling f...\n", p+1, np);
    f ();
    MPI_Finalize ();
}

