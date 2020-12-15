// #include "methods/fixed_step.hpp"
// #include "methods/adaptive.hpp"
// #include "problems/arenstorf.hpp"
// using namespace fixed_step;
// using namespace adaptive;
// using namespace arenstorf;
// using namespace std;

// int num_it = 18e2;
// double tolerance = 1e-8;
// double min_step = 1e-12;
// double max_step = 0.1;

// void print(Vector4d& x) {
//     cout << "p(" << x[0] << ", " << x[1] << ") v(" << x[2] << ", " << x[3] << ")\n";
// }

// int main() {
//     // cin >> step_size >> num_it;
//     Vector4d x0{0.994, 0, 0, -2.001585106};
//     // RichardsonExtrapolation<derivative, RK4<derivative>>
//     RKFelhberg<derivative>
//         method(tolerance, min_step, max_step);
//     auto sol = adaptive_step_ode_solver(
//         0,
//         x0,
//         method,
//         num_it
//     );
//     for (auto v : sol.second) {
//         cout << v[0] << " " << v[1] << "\n";
//     }
// }
