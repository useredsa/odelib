#include "methods/fixed_step.hpp"
#include "problems/rigid_1d.hpp"
using namespace fixed_step;
using namespace rigid_1d;
using namespace std;

double step_size = 0.0001;
double max_t = 8;
int niters = int(max_t / step_size);
analytical_sol as;

int main() {
    Vectord<1> x0{-1};
    auto sol = fixed_step_ode_solver(
        0,
        x0,
        RK4<derivative>(),
        step_size,
        niters
    );
    double err = 0.0;
    double h = 0.0;
    for(auto x : sol) {
      err = max(err, (as(h)-x).norm());
      h += step_size;
    }
    cout << err << endl;
}

