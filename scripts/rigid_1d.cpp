#include "methods/fixed_step.hpp"
#include "methods/implicit.hpp"
#include "problems/rigid_1d.hpp"
#include "utils/general.hpp"
using namespace implicit;
using namespace fixed_step;
using namespace rigid_1d;
using namespace std;

double step_size = 0.21;
double max_t = 8;
int niters = int(max_t / step_size);

int main() {
    cin >> step_size;
    Vectord<1> x0{-1};
    vector<Vectord<1>, Eigen::aligned_allocator<Vectord<1>>> as(niters+1);
    vector<double> t(niters+1);
    for (int i = 0; i < niters+1; ++i) {
        as[i] = analytical_sol()(step_size*i);
        t[i] = step_size*i;
    }
    {
        auto sol = fixed_step_ode_solver(
            0,
            x0,
            Euler<derivative>(),
            step_size,
            niters
        );
        double err = abs_diff_vectors(sol, as);
        cout << err << " " << abs_diff_func<1, analytical_sol>(t, sol) << endl;
    }

    {
        auto sol = newton_ode_solver(
            0,
            x0,
            BackwardsEuler<derivative>(),
            step_size,
            niters
        );
        double err = abs_diff_vectors(sol, as);
        cout << err << " " << abs_diff_func<1, analytical_sol>(t, sol) << endl
    }

    {
        auto sol = secant_ode_solver(
            0,
            x0,
            BackwardsEuler<derivative>(),
            step_size,
            niters
        );
        double err = abs_diff_vectors(sol, as);
        cout << err << " " << abs_diff_func<1, analytical_sol>(t, sol) << endl;
    }
}

