# C++ Numerical Solvers for Ordinary Differential Equations

Welcome to our numerical library for ODEs.
This repository contains C++ implementations of common ODE numerical methods.

We chose C++ because it is a language that lets us
write high-level code using vector operations,
while performing **just as fast as plain C style** for loops
(pending upload of tests).
ODE Solvers are generally difficult to optimize for the compiler
because to be calculated,
each point depends on the previous ones and a derivative function that is
generally unknown in other languages
during the compilation process of the method,
making it difficult to vectorize.
We make use of the C++ template system to pass all possible information
in compilation time.
The result is a generation of **code specific for each intial value problem**
and **incredibly fast**.

## Target Audience

We developed this library for the subject
*Numerical Methods for Differential Equations*
which we taked in our Mathematics and Computer Engineering degrees in
the [University of Murcia].

Our target audience are Computer Engineering and/or Mathematics students
that are also learning and implementing this algorithms.
This may have influenced some of our decisions.
For example, the project is very easy to download and test because
our main dependency comes in the form of header files,
avoiding the necessity to install external libraries in your system,
and our build system is a [`Makefile`].

[University of Murcia]: https://www.um.es/

## Compiling and Running

You can find the headers used for the scripts under [include],
and their source files under [src].
We call *scripts* small programs that make use of this functionality
to apply methods to specific problems.
Scripts are files that contain a `main` function and
can be compiled to an executable.

We use the [make] utility to build this project
because it is readily available in most systems.
You can compile this project using the following commands.

```sh
make                   # compiles all scripts
make bin/<script_name> # compiles an specific script
make clean             # removes all the files generated during compilation
make cleandep          # deletes all dependency files
make info              # prints makefile debugging information
```

The [`Makefile`] automatically detects dependencies in C++ source files.
Therefore, you can add more headers or scripts to this project
and compile with the same commmands.
Aftter compiling, run any program with `./bin/<script_name>`.

[`Makefile`]: /Makefile
[include]: /include
[src]: /src
[make]: https://www.gnu.org/software/make/

## Dependencies

### Eigen

[Eigen] is a C++ template library for linear algebra:
matrices, vectors, numerical solvers, and related algorithms.
We use Eigen to code high level algorithms
using their powerful system of vectors and matrices
while maintaining peak performance.

Eigen is licensed under the [MPL2].
We include some parts of Eigen 3.3.8 under [lib] as is.

[Eigen]: http://eigen.tuxfamily.org/index.php?title=Main_Page
[lib]: /lib
[MPL2]: https://www.mozilla.org/en-US/MPL/2.0/

## LICENSE

The libraries used come with their own software licenses and apply only
to that source code.

Our own work is licensed under the [Apache 2.0 License].
This is a permissive license —
meaning that it grants you (and anyone) permission to
use, copy, modify, distribute, and run this code.
Had there not been a license, you would not have that permission.
If you find useful our work, we will be glad to hear from you!

[Apache 2.0 License]: https://apache.org/licenses/LICENSE-2.0
