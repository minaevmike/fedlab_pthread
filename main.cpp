#include <iostream>
#include <vector>
#include <getopt.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int delta_t;
    int num_threads;
    int x_grid, y_grid;
    const char* short_options = "tpxy";
    int c;
    while (1) {
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
    std::cout << delta_t << " " << num_threads << " " << x_grid << " " << y_grid<< std::endl;
    return 0;
}

