#include "methods/fixed_step.hpp"
#include "problems/arenstorf.hpp"
using namespace fixed_step;
using namespace arenstorf;
using namespace std;

double K_STEP = 1e-6;
int K_NUMIT = 40e6;

void print(Vectord<K_DIM>& x) {
    cout << "p(" << x[0] << ", " << x[1] << ") v(" << x[2] << ", " << x[3] << ")\n";
}

int main() {
    // cin >> K_STEP >> K_NUMIT;
    Vectord<K_DIM> x0{0.994, 0, 0, -2.001585106}; 
    auto sol = fixed_step_ode_solver(
        0,
        x0,
        RK4<derivative, K_DIM>(),
        K_STEP,
        K_NUMIT
    );
    // for (auto x : sol) {
    //     cout << "p(" << x[0] << ", " << x[1] << ") v(" << x[2] << ", " << x[3] << ")\n";
    // }
    print(sol.back());
}

