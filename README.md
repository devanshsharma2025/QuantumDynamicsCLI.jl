# QuantumDynamicsCLI

## What is QuantumDynamicsCLI?

Simulating the dynamics of quantum systems is a challenging task with a
multitude of complicated computational methods. The
[QuantumDynamics.jl](https://github.com/amartyabose/QuantumDynamics.jl) package
provides modular open-source implementations of an increasingly growing number
of these methods, while remaining a flexible platform for further development.
However, owing primarily to its exceptionally flexible nature, the usage of
QuantumDynamics.jl happens through short Julia scripts. This means that for the
most common simulation jobs, in addition to writing data post-processing
scripts, one also needs to effectively rewrite the same simulation code multiple
times increasing the chances of errors. 

QuantumDynamicsCLI.jl installs a command-line application, `qdsim`, that
leverages the QuantumDynamics.jl package to provide the user with an easy-to-use
interface for doing simulations without requiring scripting. After
QuantumDynamicsCLI.jl has been installed and built, `qdsim` is usually placed in
`$HOME/.julia/bin`. Please add this to your `PATH` environment variable.

## Basic Usage

`qdsim` follows a modular structure. The general syntax is as follows:

```bash
> qdsim <component_name> <command_name> <arguments>
```

Currently, two `command`s are supported:
- `simulate`: gives access to various techniques for simulating the dynamics
- `post`: provides post-processing tools for the output

The most important sub-command of `simulate` is `run`, and for `post` is `get-observable`.

## Examples

A series of examples are provided in the [examples](./examples/) folder. Below we list the broad ideas covered:

1. Path Integral Simulation of Dynamics:
    1. [Spin-Boson Problem](./examples/01-Spin-Boson/)
2. Hierarchical Equations of Motion:
    1. [Fenna-Matthews-Olson Complex](./examples/02-Ishizaki-Fleming-FMO/)