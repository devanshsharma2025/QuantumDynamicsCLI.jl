# QuantumDynamicsCLI

<!-- | **Documentation** | **Build Status** | **Citation** |
|:-----------------:|:---------:|:-------------:|
|[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amartyabose.github.io/QuantumDynamics.jl/dev/)|[![Run tests](https://github.com/amartyabose/QuantumDynamics.jl/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/amartyabose/QuantumDynamics.jl/actions/workflows/test.yml)|[![DOI](https://img.shields.io/badge/DOI-10.1063/5.0151483-blue.svg)](https://doi.org/10.1063/5.0151483)| -->

## What is QuantumDynamicsCLI?
Simulating the dynamics of quantum systems is a challenging task with a
multitude of complicated computational methods. The
[![QuantumDynamics.jl](https://github.com/amartyabose/QuantumDynamics.jl)]
package attempts to provide modular open-source implementations of an
increasingly growing number of these methods, while remaining a flexible
platform for further development. However, owing primarily to its exceptionally
flexible nature, the usage of QuantumDynamics.jl happens through short Julia
scripts. This means that for the most common simulation jobs, one needs to
effectively rewrite the same code multiple times increasing the chances of
errors. As a means to making some of the common types of simulations more
facile, we now offer the QuantumDynamicsCLI.jl package which installs the qdsim
application as a sister code of the QuantumDynamics.jl package. As the framework
grows, so will this application grow to accommodate the new methods and their
most common use cases.

## Installation
QuantumDynamicsCLI.jl is a registered package. Installation is a simple procedure. It can be done either through the Pkg REPL:
```bash
~ julia
```

```
julia> ]
pkg> add QuantumDynamicsCLI
```

or by using the `Pkg` package manager in a script as follows:
```julia
julia> using Pkg
julia> Pkg.add("QuantumDynamicsCLI")
```

After the package gets built, an executable called `qdsim` will be placed in `$HOME/.julia/bin` along with the code completions for the shell in `$HOME/.julia/completions`. Please add `$HOME/.julia/bin` to your path and source the correct completions file.

While QuantumDynamicsCLI.jl builds on top of the QuantumDynamics.jl package, separate installation of that package is unnecessary. Just installing QuantumDynamicsCLI.jl would install QuantumDynamics.jl as a dependency.