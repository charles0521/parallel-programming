#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include <ctime>
#define MAX 1.0
#define MIN -1.0

void count_pi(long long int number_of_tosses, long long int &number_in_circle)
{
    double distance_squared;

    unsigned int seed = (unsigned int)87;
    for (int i = 0; i < number_of_tosses; ++i)
    {
        double x = (double)(MAX - MIN) * rand_r(&seed) / (RAND_MAX + 1.0) + MIN;
        double y = (double)(MAX - MIN) * rand_r(&seed) / (RAND_MAX + 1.0) + MIN;
        distance_squared = x * x + y * y;
        if (distance_squared <= 1)
            number_in_circle++;
    }

}

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status status;

    // calculate parameter
    long long int local_tosses = tosses / world_size;
    long long int local_count = 0;

    // pass config
    unsigned int dest=0;
    int tag = 0;

    // count pi
    count_pi(local_tosses, local_count);

    // TODO: binary tree redunction
    int tree_height = 0, temp=world_size;
    int base = 1;
    while(temp > 1)
    {
        temp /= 2;
        tree_height++;
    }
    for(int i=0; i<tree_height; ++i)
    {
        base *= 2;

        MPI_Barrier(MPI_COMM_WORLD);
        if(world_rank % base == 0)
        {
            long long int buf;
            MPI_Recv(&buf, 1, MPI_UNSIGNED_LONG, world_rank+base/2, tag, MPI_COMM_WORLD, &status);
            local_count += buf;
        }
        else if( (world_rank + base/2) % base == 0)
        {
            MPI_Send(&local_count, 1, MPI_UNSIGNED_LONG, world_rank-(base/2), tag, MPI_COMM_WORLD);
            // return 0;
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0)
    {
        // TODO: PI result
        pi_result = 4 * (local_count) /(( double ) tosses);
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
