# Simple Spin-Boson example

In this folder, we have all the files necessary for running a typical spin-boson calculation using the TNPI-TTM method.

1. Run the path integral calculation for the dynamical maps:
    ```bash
    qdsim simulate run system.toml simulate.toml
    ```
2. Propate the initial density matrix localized on the first site:
    ```bash
    qdsim simulate propagate-using-tmats system.toml propagate.toml
    ```
3. Finally, use the dynamics obtained to generate observables:
    ```bash
    qdsim post get-observable system.toml observables.toml
    ```