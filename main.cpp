
#include <iostream>
#include <cmath>

template<size_t rows, size_t cols>

double lambda(double x, double y);

double lambda(double x, double y);

using namespace std;

template<size_t rows, size_t cols>
void solve_matrix_x(double (&u)[rows][cols], double (&u_prev)[rows][cols], int n, double h, int cur_index,
                    char direction, double *x, double *y) {
    auto *alpha = new double[n + 1];
    auto *betta = new double[n + 1];

    alpha[0] = 0;
    betta[0] = 500;

    auto *a = new double[n];
    auto *b = new double[n];
    auto *c = new double[n];
    auto *f = new double[n];

    // init a, b, c, f
    for (int i = 1; i < n; i++) {
        if (direction == 'u') {
            a[i] = -lambda(x[cur_index - 1], y[i]) / (2 * h * h);
            b[i] = -lambda(x[cur_index + 1], y[i]) / (2 * h * h);
            f[i] = u_prev[cur_index][i] / h +
                   (lambda(x[cur_index], y[i + 1]) * (u_prev[cur_index][i + 1] - u_prev[cur_index][i])
                    - lambda(x[cur_index], y[i - 1]) * (u_prev[cur_index][i] - u_prev[cur_index][i - 1]))
                   / (2 * h * h);
        } else {
            a[i] = -lambda(x[i], y[cur_index - 1]) / (2 * h * h);
            b[i] = -lambda(x[i], y[cur_index + 1]) / (2 * h * h);
            f[i] = u_prev[cur_index][i] / h +
                   (lambda(x[i + 1], y[cur_index]) * (u_prev[i + 1][cur_index] - u_prev[cur_index][i])
                    - lambda(x[i - 1], y[cur_index]) * (u_prev[i][cur_index] - u_prev[i - 1][cur_index]))
                   / (2 * h * h);
        }
        // h here equals tau
        c[i] = 1.0 / h - a[i] - b[i];
        cout << a[i] << " " << b[i] << " " << c[i] << " " << f[i] << endl;


    }

    // forward sweep
    for (int i = 1; i < n; i++) {
        alpha[i] = -b[i] / (c[i] + a[i] * alpha[i - 1]);
        betta[i] = f[i] - a[i] * betta[i - 1] / (c[i] + a[i] * alpha[i - 1]);

    }

    // back substitution
    for (int i = n - 1; i >= 1; i--) {
        if (direction == 'u') {
            u[cur_index][i] = alpha[i - 1] * u[cur_index][i + 1] + betta[i - 1];
        } else {
            u[i][cur_index] = alpha[i - 1] * u[i + 1][cur_index] + betta[i - 1];
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            u_prev[i][j] = u[i][j];
        }
    }
}

double lambda(double x, double y) {
    if (x >= 0.25 && x <= 0.65 && y >= 0.1 && y <= 0.25) {
        return 0.01;
    }
    return 0.0001;
}


template<size_t rows, size_t cols>
void solveEqProgonka(int n_t, int n_x, int n_y, double (&u)[rows][cols], double (&u_prev)[rows][cols],
                     double *x, double *y, double hx, double hy) {

    for (int t = 1; t < n_t; t++) {
        // y(j) direction
        for (int i = 1; i < n_x; i++) {
            solve_matrix_x(u, u_prev, n_y, hy, i, 'u', x, y);
        }

        // x(i) direction
        for (int j = 1; j < n_y; j++) {
            solve_matrix_x(u, u_prev, n_x, hx, j, 'r', x, y);
        }
        // print results on current time step
        for (int i = 0; i <= n_x; i++) {
            for (int j = 0; j <= n_y; j++) {
                cout << u[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
};

double u_left_border(double y) {
    return 600;
}

double u_right_border(double y) {
    return 1200;
}

double u_bottom_border(double x) {
    return 600 * (1 + x);
}

double u_upper_border(double x) {
    return 600 * (1 + x * x * x);
}

int main() {

    const int n_t = 1000;
    const int n_y = 100;
    const int n_x = 200;
    const double x_min = 0;
    const double x_max = 1;
    const double y_min = 0;
    const double y_max = 0.5;
    const double t_min = 0;
    const double t_max = 10;

    double hx = (x_max - x_min) / n_x;
    double hy = (y_max - y_min) / n_y;
    double ht = (t_max - t_min) / n_t;
    double *x = new double[n_x + 1];
    double *y = new double[n_y + 1];
    double u[n_x][n_y] = {{}};
    double u_prev[n_x][n_y] = {{}};

    double a[n_x + 1][n_y + 1];
    double b[n_x + 1][n_y + 1];
    double c[n_x + 1][n_y + 1];
    double f[n_x + 1][n_y + 1];


    x[0] = 0;
    for (int i = 1; i <= n_x; i++)
        x[i] = x[i - 1] + hx;
    y[0] = 0;
    for (int i = 1; i <= n_y; i++)
        y[i] = y[i - 1] + hy;

    // T(x,y,0) = 300
    for (int i = 0; i <= n_x; i++) {
        for (int j = 0; j <= n_y; j++) {
            u[i][j] = 300;
        }
    }
    // assign border values of equation
    for (int i = 0; i <= n_x; i++) {
        u[i][0] = u_bottom_border(x[i]);
        u[i][n_y] = u_upper_border(x[i]);
    }
    for (int j = 0; j <= n_y; j++) {
        u[0][j] = u_left_border(y[j]);
        u[n_x][j] = u_right_border(y[j]);

    }
    // copy u to u_prev
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            u_prev[i][j] = u[i][j];
        }
    }

    solveEqProgonka(n_t, n_x, n_y, u, u_prev, x, y, hx, hy);

    return 0;
}