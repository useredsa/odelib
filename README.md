# C++ Numerical Solvers for Ordinary Differential Equations

Pending description.

## Target Audience

Pending description.

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

Eigen is a C++ template library for linear algebra:
matrices, vectors, numerical solvers, and related algorithms.
We use Eigen to code high level algorithms
using their powerful system of vectors and matrices
while maintaining peak performance.

Eigen is licensed under the [MPL2].
We include some parts of Eigen 3.3.8 under [lib] as is.

[lib]: /lib
[MPL2]: https://www.mozilla.org/en-US/MPL/2.0/
