// #include "methods/fixed_step.hpp"
// #include "problems/two_bodies.hpp"
// using namespace fixed_step;
// using namespace two_bodies;
// using namespace std;

// double step_size = 1;
// int years = 2;
// int hours = years * 365 * 24;
// int num_it = (int) (hours / step_size);

// void print(Vector4d& x) {
//     cout << "p(" << x[0] << ", " << x[1] << ") v(" << x[2] << ", " << x[3] << ")\n";
// }

// int main() {
//     Vector4d x0{152.100533, 0, 0, 0.105444};
//     auto sol = fixed_step_ode_solver(
//         0,
//         x0,
//         RK4<derivative>(),
//         step_size,
//         num_it
//     );
//     for(auto x : sol)
//         cout << x[0] << " " << x[1] << endl;
// }

