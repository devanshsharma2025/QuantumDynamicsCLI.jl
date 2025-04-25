# Inputs to qdsim

A simulation depends upon the specification of the system. This is done in the `system` TOML file, which specifies three different things:
- the `units` in use
- the Hamiltonian for the problem being simulated:
    - the description of the `system`
    - the description of the `bath` or environment

Along with this `system` TOML file, a `simulation` TOML file also needs to be prepared that provides the details of the simulation.

## System File
Each of the three sections in the `system` file has a dedicated function for parsing the data. Before giving a full example of a `system` input file, we discuss the parameters accepted by the various sections.

### Specifying the Units
```@docs
QuantumDynamicsCLI.ParseInput.parse_unit
```

### Specifying the Hamiltonian
The basic problem under study has a general form, $$\hat{H} = \hat{H}_0 + \hat{H}_\text{env}$$, where $\hat{H}_0$ forms the system Hamiltonian and the $\hat{H}_\text{env}$ is the environment or bath Hamiltonian.

#### System Hamiltonian
```@docs
QuantumDynamicsCLI.ParseInput.parse_system
```

#### Bath Hamiltonian
```@docs
QuantumDynamicsCLI.ParseInput.parse_bath
```

```@docs
QuantumDynamicsCLI.ParseInput.get_bath
```

## Simulation File
There are broadly three types of simulations that are currently supported:
- `dynamics` simulations that simulate the non-equilibrium dynamics of the given problem
- `equilibrium_rho` simulations that simulate the equilibrium density at the given temperature
- `complex_corr` simulations for calculating equilibrium correlation functions of various flavors

### Dynamics Simulations
Various methods of simulation are supported:
- Path Integral Methods:
    - Quasi-adiabatic Propagator Path Integrals
    - TEMPO
- HEOM
- Generalized Quantum Master Equation
- Multichromophore Incoherest Forster Theory
- Bloch-Redfield Master Equation

All of these methods require some core common parameters and then more specfic method-dependent parameters. The core parameters of all the dynamics methods are:
- `dt`: for the time-step
- `nsteps`: for the number of steps that the dynamics needs to be calculated for

## Example
### Spin-Boson Example
Let us say we are trying to simulate a typical Spin-Boson problem, where all parameters are specified in atomic units.

```toml
[system]
Hamiltonian = "Hamiltonian"

[baths]
beta = 5.0
[[baths.bath]]
type = "ohmic"
xi = 0.1
omegac = 7.5
svec = [1.0, -1.0]
```

Notice that the `[units]` section is completely skipped over because the default values specify atomic units.