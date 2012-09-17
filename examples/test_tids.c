/* Tests if the thread system satisfies the intervalic thread ids property.
 * Pass the number of threads to spawn.
 * If so, the output should look like: 2 3 4 5 6   or some other interval.
 */

#include <stdio.h>
#include "mpii.h"

Thread *ids;
int done = 1;
Mutex mutex;
Cond cond;

void *fun (int i)
{
    printf ("hi this is thread %d\n", i);
    ids[i] = thread_id ();
    lock (mutex);
    done++;
    notify (cond);
    unlock (mutex);
}

main (int argc, char *argv[])
{
    int i, n;

    init_threads ();
    new_mutex (mutex);
    new_cond (cond);

    if (argc != 2)
    {
        fprintf (stderr, "Usage: test_tids nthread     e.g. 5\n");
        exit (1);
    }

    n = atoi (argv[1]);
    if (n <= 0)
    {
        fprintf (stderr, "nthread must be at least 1\n");
        exit (1);
    }

    ids = (Thread *) malloc (sizeof (Thread) * n);
    if (ids == NULL)
    {
        perror ("malloc");
        exit (1);
    }

    for (i = 0; i < n; i++)
        ids[i] = -1;

    for (i = 0; i < n; i++)
    {
        Thread tid;
        if (spawn_thread (fun, (void *) i, tid))
        {
            i++;
            break;
        }
    }

    lock (mutex);
    while (done < i)
        wait (cond, mutex);
    unlock (mutex);

    printf ("%d", thread_id ());
    for (i = 1; i < n; i++)
        printf (" %d", ids[i]);
    printf ("\n");
}
