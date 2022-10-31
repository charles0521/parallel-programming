#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <pthread.h>

using namespace std;

#define MAX 1.0
#define MIN -1.0

pthread_mutex_t mutexsum;

typedef struct thread_arg
{
    int thread_id;
    long long int part;
    long long int *number_in_circle;
} arg;

void *count_pi(void *th_arg)
{
    double distance_squared;
    arg* data = (arg*)th_arg;
    long long int number_of_tosses = data->part;
    long long int *number_in_circle = data->number_in_circle;

    long long int local_number_in_circle = 0;
    unsigned int seed = (unsigned int)pthread_self();
    
    for (int i = 0; i < number_of_tosses; ++i)
    {
        double x = (double)(MAX - MIN) * rand_r(&seed) / (RAND_MAX + 1.0) + MIN;
        double y = (double)(MAX - MIN) * rand_r(&seed) / (RAND_MAX + 1.0) + MIN;
        distance_squared = x * x + y * y;
        if (distance_squared <= 1)
            local_number_in_circle++;
    }
    
    // critical section
    pthread_mutex_lock(&mutexsum);
    *number_in_circle += local_number_in_circle;
    pthread_mutex_unlock(&mutexsum);

    pthread_exit((void *)0);

}

int main(int argc, char *argv[])
{

    srand(time(NULL));

    pthread_mutex_init(&mutexsum, NULL);

    // 設定 pthread 性質是要能 join
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // 每個 thread 都可以存取的 PI
    // 因為不同 thread 都要能存取，故用指標
    long long int *number_in_circle = (long long int *)malloc(sizeof(long long int));
    *number_in_circle = 0;

    if (argc != 3)
    {
        throw out_of_range("number of arguments error");
    }

    // get argv
    const long long int number_of_tosses = atoll(argv[2]);
    const int number_of_thread = atoi(argv[1]);

    // each thread toss
    long long int part = number_of_tosses / number_of_thread;

    // def number of thread
    pthread_t thread[number_of_thread];

    // set thread arg
    arg thread_arg[number_of_thread];

    // create thread
    // auto begin = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < number_of_thread; ++i)
    {
        thread_arg[i].thread_id = i;
        thread_arg[i].part = part;
        thread_arg[i].number_in_circle = number_in_circle;

        pthread_create(&thread[i], &attr, count_pi, (void *)&thread_arg[i]);
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    // printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
    
    // destroy attr
    pthread_attr_destroy(&attr);

    //wait all of thread
    void *status;
    for(int i=0; i<number_of_thread; ++i)
    {
        pthread_join(thread[i], &status);
    }

    double pi_estimate  = 4 * (*number_in_circle) /(( double ) number_of_tosses);

    cout << pi_estimate << endl;

    pthread_mutex_destroy(&mutexsum);
    return 0;
}