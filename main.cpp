#include <iostream>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include <cmath>
std::vector<std::vector<double> > phi;
std::vector<int> condition;

void print_matrix(std::vector<std::vector<double> > v) {
    for(int i = 0; i < v.size(); i++) {
        for(int j = 0; j < v.size(); j++) {
            std::cout << v[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_matrix(std::vector<double> v) {
    for(int i = 0; i < v.size(); i++){ 
        std::cout << v[i] << std::endl;
    }
}

std::vector<double> solve_progon(std::vector<std::vector<double> > A, std::vector <double> B) {
    std::vector<double> P, Q, ans;
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

std::vector<double> check_solution(std::vector<std::vector<double> >A , std::vector<double> X, std::vector<double> B) {
    std::vector<double> diff;
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

int main(int argc, char **argv) {
    int time_to_modulate;
    int num_threads;
    int x_grid, y_grid;
    const char* short_options = "tpxy";
    int c;
    double C = 1e-6;
    double R = 10;
    double delta_t = 1e-5;
    x_grid = 2;
    y_grid = 2;
    num_threads = 2;
    time_to_modulate = 2;
    phi.resize(x_grid);
    condition.resize(x_grid);
    for(int i = 0; i < phi.size(); i++) {
        phi[i].resize(y_grid);
    }
    std::vector<std::vector<double> > A;
    std::vector<double> B, X;
    int size = atoi(argv[1]);;
    A.resize(size);
    for(int i = 0; i < A.size(); i++) {
        A[i].resize(size);
    }
    B.resize(size);
    for(int i = 0; i < A.size(); i++) {
        if(i - 1 >=0)
            A[i][i-1] = (double)rand() / RAND_MAX;
        A[i][i] = (double)rand() / RAND_MAX;
        
        if(i + 1 < A.size())
            A[i][i+1] = (double)rand() / RAND_MAX;
        B[i] = (double)rand() / RAND_MAX;
        
    }
    //print_matrix(A);
    //print_matrix(B);
    X = solve_progon(A, B);
    std::vector<double> check = check_solution(A,X,B);
    for (int i = 0; i < check.size(); i++) {
        if (fabs(check[i]) > 1e-10)
            std::cout << "ERROR TOO BIG" << check[i] << std::endl;
    }
    //print_matrix(check_solution(A, X, B));
    /*while (1) {
        const struct option long_options[] = {
            {"delta_t", required_argument, 0, 0},
            {"threads", required_argument, 0, 0},
            {"x_size", required_argument, 0, 0},
            {"y_size", required_argument, 0, 0},
            {NULL, 0, NULL, 0}
        };
        int index = 0;
        c = getopt_long(argc, argv, "t:p:x:y:", long_options, &index);
        if (c == -1)
            break;
        std::cout << c << " " << optarg << std::endl;
        switch (c) {
            case 't':
                delta_t = atoi(optarg);
                break;
            case 'p':
                num_threads = atoi(optarg);
                break;
            case 'x':
                x_grid = atoi(optarg);
                break;
            case 'y':
                y_grid = atoi(optarg);
                break;
        }
    }
    std::cout << delta_t << " " << num_threads << " " << x_grid << " " << y_grid<< std::endl;*/

    return 0;
}

