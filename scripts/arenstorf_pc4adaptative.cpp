#include "methods/fixed_step.hpp"
#include "methods/multistep.hpp"
#include "methods/adaptative_multistep.hpp"
#include "problems/arenstorf.hpp"

using namespace fixed_step;
using namespace adaptative_multistep;
using namespace multistep;
using namespace arenstorf;
using namespace std;

double step_size = 1e-6;
int num_it = 40e6;

void print(Vector4d& x) {
    cout << "p(" << x[0] << ", " << x[1] << ") v(" << x[2] << ", " << x[3] << ")\n";
}

int main() {
    // cin >> step_size >> num_it;
    double t0 = 0;
    Vector4d x0{0.994, 0, 0, -2.001585106}; 
    struct PredictorCorrector4<
        derivative,
        AdamsBashforth4<derivative>,
        RK4<derivative>
        > pc4;
    pc4.h = 1;
    pc4.tol = 1e-8;
    auto [t, x] = adaptative_multistep_ode_solver(
        t0,
        x0,
        pc4,
        18.0
    );
    // print(sol.second.back());
    for(size_t i=0; i<t.size(); ++i){
        cout << x[i][0] << " " << x[i][1] << endl;
    }
}

