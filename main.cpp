#include <iostream>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include <cmath>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
#include <utility>
#include <fstream>
#ifndef NOGNUPLOT
    std::string data_file_name = "data.out";
    std::ofstream out(data_file_name.c_str(), std::ofstream::out | std::ofstream::trunc);
    std::vector<std::pair<int,int> > gnuplot_points;
#endif
typedef std::vector<std::vector<double> > matrix;
typedef std::vector<double> row;
std::vector<std::vector<double> > phi, I;
std::vector<int> condition;
static int THREADS, X_SIZE, Y_SIZE;
pthread_barrier_t bar, syn, val;
pthread_mutex_t lock;
static double R = 1, C = 1e-1, dt, T;
long long clock_time() {
    return time(NULL);
}
void print_matrix(const matrix &v) {
	int size = v.size();
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
    double pt, bt;
	row B(size);
    matrix A(size, row(size));
    for(int i = 0; i < A.size(); i++) {
        if(i - 1 >=0)
            A[i][i-1] = 1;
        A[i][i] = -2;
        if(i + 1 < A.size())
            A[i][i+1] = 1;
        if(i == 0)
            pt = 0;
        else
            pt = mid[i-1];
        if (i == A.size() - 1)
            bt = 0;
        else 
            bt = mid[i+1];
        B[i] = (-left[i] + 2 * mid[i] - right[i]) * dt /(R * C) + (pt - 2 * mid[i] + bt) + I[j][i] * dt/C;//(left[i] - 2 * mid[i] + right[i]) * (dt / (R * C) - 1) + I[j][i] * dt / C ;
    }
    row s = solve_progon(A, B);
    return s;
}

void * calc(void *thread) {
    long t = (long) thread;
    int start_index = t * (X_SIZE / THREADS); 
    int end_index = (t != THREADS - 1)?(t + 1) * (X_SIZE / THREADS) - 1: X_SIZE - 1;
    matrix local(phi.size()), next(phi.size());
    row zeros;
    for (int i = 0; i < Y_SIZE; i++) {
        zeros.push_back(0);
    }
    double cur_time = 0, a;
    int x, y;
	unsigned long long iters = 0;
    std::cout << "I am thread " << t << ". My start index: " << start_index << ". My end index: " << end_index << std::endl;
	row left_t, mid_t, right_t;
    unsigned long long start = clock_time();
    while (cur_time < T) {
    //for(int i = 0; i < 3; i ++) {
        pthread_barrier_wait(&syn);
#ifndef NOGNUPLOT
        if (t == 0) {
            out << cur_time << " ";
            for (int i = 0; i < gnuplot_points.size(); i++) {
                x = gnuplot_points[i].first;
                y = gnuplot_points[i].second;
                if (x >= 0 && x < X_SIZE && y >=0 && y < Y_SIZE)
                    out << phi[x][y] << " ";
            }
        }
        out << std::endl;
#endif
		cur_time += dt;
        for (int i = start_index; i <= end_index; i ++) {
			iters ++;
			mid_t = phi[i];
            if (i == start_index && start_index == 0) {
				left_t = zeros;
				right_t = phi[i+1];
			}
            else if (i == end_index && end_index == X_SIZE - 1) {
				left_t = phi[i - 1];
				right_t = zeros;
			}
            else {
				left_t = phi[i - 1];
				right_t = phi[i + 1];
			}
			next[i] = solve(left_t, mid_t, right_t, i);
        }
        pthread_barrier_wait(&bar);
        for (int i = start_index; i <=end_index; i++) {
            phi[i] = next[i];
        }
    }
    pthread_exit(NULL);
}

int main(int argc, char **argv) {
    int time_to_modulate;
    int num_threads;
    int x_grid, y_grid;
    const char* short_options = "tpxy";
    int c;
    T = atoi(argv[4]);
    dt = 1e-5;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    x_grid = atoi(argv[1]);
    y_grid = atoi(argv[2]);
    num_threads = atoi(argv[3]);
#ifndef NOGNUPLOT
    int print_points = atoi(argv[5]);
    std::cout << "Enter " << print_points << " nodes to print in format <m> <n>" << std::endl;
    int x, y;
    for(int i = 0; i < print_points; i++) {
        std::cin >> x >> y;
        gnuplot_points.push_back(std::make_pair(x, y));
    }
#endif
    THREADS = num_threads;
    time_to_modulate = 2;
    X_SIZE = x_grid;
    Y_SIZE = y_grid;
    phi.resize(x_grid);
    for(int i = 0; i < phi.size(); i++) {
        phi[i].resize(y_grid);
    }
    I.resize(x_grid);
    for(int i = 0; i < I.size(); i++) {
        I[i].resize(y_grid);
		for(int j = 0; j < y_grid; j++)
			I[i][j] = 0;// rand() % 6 - 3;
    }
   // I[1][1] = 1;
   // I[0][1] = -1;
    pthread_barrier_init(&bar, NULL, THREADS);
    pthread_barrier_init(&val, NULL, THREADS);
    pthread_barrier_init(&syn, NULL, THREADS);
    pthread_mutex_init(&lock, NULL);
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
    std::cout << "It takes " << (clock_time() - start) << std::endl;
#ifndef NOGNUPLOT
    std::ofstream gnuplot("cmd", std::ofstream::trunc | std::ofstream::out);
    gnuplot << "plot ";
    for (int i = 0; i < gnuplot_points.size(); i++) {
        gnuplot << "\"" << data_file_name << "\" using 1:" << i + 2 << "with lines title \"phi" << gnuplot_points[i].first << "," << gnuplot_points[i].second << "\"";
        if(i != gnuplot_points.size() - 1)
            gnuplot << ",";
    }
    gnuplot.close();
    int res = system("gnuplot -persist cmd");
#endif
    return 0;
}

