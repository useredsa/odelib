#include "methods/fixed_step.hpp"
#include "problems/arenstorf.hpp"
using namespace fixed_step;
using namespace arenstorf;
using namespace std;

double step_size = 1e-6;
int num_it = 40e6;

void print(Vector4d& x) {
    cout << "p(" << x[0] << ", " << x[1] << ") v(" << x[2] << ", " << x[3] << ")\n";
}

int main() {
    // cin >> step_size >> num_it;
    Vector4d x0{0.994, 0, 0, -2.001585106};
    auto sol = fixed_step_ode_solver(
        0,
        x0,
        RK4<derivative>(),
        step_size,
        num_it
    );
    print(sol.back());
}

