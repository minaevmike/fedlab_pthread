#include <QCoreApplication>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include <cmath>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
#include <QDateTime>
#include <QDebug>
typedef std::vector<std::vector<double> > matrix;
typedef std::vector<double> row;
std::vector<std::vector<double> > phi, I;
std::vector<int> condition;
static int THREADS, X_SIZE, Y_SIZE;
pthread_barrier_t bar, syn;
pthread_mutex_t lock;
static double R = 10, C = 1e-6, dt, T;
long long clock_time() {
    struct timespec tp;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
}
void print_matrix(const matrix &v) {

    for(int i = 0; i < v.size(); i++) {
        for(int j = 0; j < v[i].size(); j++) {
            std::cout << v[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_matrix(const row &v) {
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << std::endl;
    }
}

row solve_progon(const matrix &A,const row &B) {
    row P, Q, ans;
    int N = B.size();
    P.resize(N);
    Q.resize(N);
    ans.resize(N);
    P[0] = - A[0][1] / A[0][0];
    Q[0] = B[0]/A[0][0];
    //Pr9moy xod
    for(int i = 1; i < N; i++) {
        P[i] = A[i][i + 1]/(-A[i][i] - A[i][i - 1]*P[i-1]);
        Q[i] = (A[i][i - 1] * Q[i - 1] - B[i]) / (-A[i][i] - A[i][i-1] * P[i - 1]);
    }

    //Obratniy xod
    ans[N - 1] = (A[N-1][N-2] * Q[N-2] - B[N-1]) / (-A[N-1][N-1] - A[N-1][N-2] * P[N-2]);
    for(int i = N - 2; i >= 0; --i) {
        ans[i] = P[i] * ans[i+1] + Q[i];
    }
    return ans;
}

row check_solution(const matrix &A , const row &X, const row &B) {
    row diff;
    diff.resize(B.size());
    int size = B.size();
    for(int i = 0; i < size; i++) {
        double temp = 0;
        for(int j = 0; j < size; j++) {
            temp += A[i][j] * X[j];
        }
        diff[i] = temp - B[i];
    }
    return diff;
}

row solve(const row &left,const row &mid,const row &right, int j) {
    int size = left.size();
    row B(size);
    matrix A(size, row(size));
    //std::cout << left.size() << " " << mid.size() << " " << right.size() << std::endl;
    for(int i = 0; i < A.size(); i++) {
        if(i - 1 >=0)
            A[i][i-1] = -1;
        A[i][i] = 2;

        if(i + 1 < A.size())
            A[i][i+1] = -1;
        B[i] = (left[i] - 2 * mid[i] + right[i]) * (dt / (R * C) - 1) + I[j][i] * dt / C ;
    }
    //print_matrix(A);
    row s = solve_progon(A, B);
    row check = check_solution(A, s, B);
    for (int i = 0; i < check.size(); i++) {
        if (fabs(check[i]) > 1e-10)
            std::cout << "ERROR TOO BIG" << check[i] << std::endl;
    }
    return s;
}

void * calc(void *thread) {
    long t = (long) thread;
    int start_index = t * (X_SIZE / THREADS);
    int end_index = (t != THREADS - 1)?(t + 1) * (X_SIZE / THREADS) - 1: X_SIZE - 1;
    matrix local, next;
    row zeros;
    for (int i = 0; i < Y_SIZE; i++) {
        zeros.push_back(0);
    }
    double cur_time = 0;
    std::cout << "I am thread " << t << ". My start index: " << start_index << ". My end index: " << end_index << std::endl;
    while (cur_time < T) {
        /*if (start_index == 0) {
            local.push_back(zeros);
        } else {
            local.push_back(phi[start_index - 1]);
        }
        //pthread_mutex_lock(&lock);
        for (int i = start_index; i <= end_index; i++) {
            local.push_back(phi[i]);
        }
        //pthread_mutex_unlock(&lock);
        if (end_index == X_SIZE - 1) {
            local.push_back(zeros);
        }
        else {
            local.push_back(phi[end_index + 1]);
        }*/
        //print_matrix(local);
        //std::cout << local.size() << std::endl;
        for (int i = start_index; i <= end_index; i ++) {
            if (i == start_index && start_index == 0)
                next.push_back(solve(zeros, phi[i], phi[i + 1], i - start_index));
            else if (i == end_index && end_index == X_SIZE - 1)
                next.push_back(solve(phi[i - 1], phi[i], zeros, i - start_index));
            else
                next.push_back(solve(phi[i - 1], phi[i], phi[i + 1], i - start_index));
        }
        cur_time += dt;
        pthread_barrier_wait(&bar);
        for (int i = start_index; i <=end_index; i++) {
            phi[i] = next[i - start_index];
        }
        next.clear();
        pthread_barrier_wait(&syn);
    }
    pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    qint64 temp = QDateTime::currentMSecsSinceEpoch();
    int time_to_modulate;
    int num_threads;
    int x_grid, y_grid;
    const char* short_options = "tpxy";
    int c;
    T = 2;
    dt = 1e-5;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    x_grid = atoi(argv[1]);
    y_grid = atoi(argv[2]);
    num_threads = atoi(argv[3]);
    THREADS = num_threads;
    time_to_modulate = 2;
    X_SIZE = x_grid;
    Y_SIZE = y_grid;
    phi.resize(x_grid);
    for(int i = 0; i < phi.size(); i++) {
        phi[i].resize(y_grid);
        for (int j = 0; j < phi[i].size(); j++)
            phi[i][j] = (double) rand() / RAND_MAX;
    }
    I.resize(x_grid);
    for(int i = 0; i < I.size(); i++) {
        I[i].resize(y_grid);
    }
    I[1][1] = 10;
    pthread_barrier_init(&bar, NULL, THREADS);
    pthread_barrier_init(&syn, NULL, THREADS);
    pthread_mutex_init(&lock, NULL);
    /*std::vector<std::vector<double> > A;
    std::vector<double> B, X;
    int size = atoi(argv[1]);;
    A.resize(size);
    for(int i = 0; i < A.size(); i++) {
        A[i].resize(size);
    }
    B.resize(size);
    */
    pthread_barrier_init(&bar, NULL, THREADS);
    pthread_t *threads = new pthread_t[THREADS];
    unsigned long long start = clock_time();
    for (long i = 0; i < THREADS; i++) {
        if (pthread_create(&threads[i], NULL, calc, (void *)i) != 0) {
            std::cout << "Can't create thread " << i << std::endl;
        }
    }
    for (int i = 0; i < THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    std::cout << "It takes " << (double)(clock_time() - start) / 1e9 << std::endl;
    temp =  QDateTime::currentMSecsSinceEpoch() - temp;
    qDebug()<< temp;
    return a.exec();
}
