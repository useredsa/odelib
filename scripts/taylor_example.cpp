#include "methods/fixed_step.hpp"
#include "problems/taylor_example.hpp"
#include "Eigen/Dense"
#include <cmath>
using namespace fixed_step;
using namespace taylor_example;
using namespace std;
using namespace Eigen;

double step_size = 0.01;
double max_t = 3;
int num_it = (int) (max_t / step_size);

inline Vectord<1> analytical_sol(double t) {
    return Vectord<1>{(t+1)*(t+1) - 0.5 * exp(t)};
}

int main() {
    Vectord<1> x0{0.5};
    auto sol = fixed_step_ode_solver(
        0,
        x0,
        Taylor<derivative, 5>(),
        step_size,
        num_it
    );
    double err = 0.0;
    double t = 0;
    for(Vectord<1>& x : sol) {
        err = max(err, (x - analytical_sol(t)).norm());
        t += step_size;
    }
    cout << err << endl;
}

