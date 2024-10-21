# The Ishizaki-Fleming Fenna-Matthews-Olson Complex example

Here we use the Hierarchical equations of motion to generate the dynamics for the FMO model studied by Ishizaki and Fleming in 2009.

1. Run the HEOM simulation of the dynamics:
    ```bash
    qdsim simulate run system.toml simulate.toml
    ```
2. Finally, use the dynamics obtained to generate observables:
    ```bash
    qdsim post get-observable system.toml observables.toml
    ```

Notice that now, unlike in the TTM-based simulation of the [spin-boson](../01-Spin-Boson/) problem, we do not generate the dynamical maps and then use them to propagate the initial conditions. We directly obtain the dynamics which is subsequently used to obtain the observables.