
# LinearResponse.jl

[![image](https://github.com/JuliaStellarDynamics/LinearResponse.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/LinearResponse.jl/)

**LinearResponse.jl** is a package written in Julia to perform Linear Response calculations. The software is described in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract). To see an example of the results from the software, see [this repository](https://github.com/michael-petersen/LinearResponse-paper/tree/main), which reproduces all the figures from Petersen et al. (2024).

---
## Installation

To fully use **LinearResponse.jl**, you will need **OrbitalElements.jl**, **FiniteHilbertTransform.jl**, and **AstroBasis.jl**. This installation assumes that you have not installed any of the four libraries; if you have already installed any of the others (globally, not locally) you may skip that portion of the installation.

The libraries under the JuliaStellarDynamics organisation are currently unregistered[^1]. To add them to your julia[^2] registry, follow these steps:

1. **Add Packages:** Use the package manager (`julia` at the command line, and then `]` in the interpreter) and execute the following command inside julia:
    ```julia
    add "git@github.com:JuliaStellarDynamics/OrbitalElements.jl.git"
    add "git@github.com:JuliaStellarDynamics/AstroBasis.jl.git"
    add "git@github.com:JuliaStellarDynamics/FiniteHilbertTransform.jl.git"
    add "git@github.com:JuliaStellarDynamics/LinearResponse.jl.git"
    ```

2. **Verify Version:** Confirm the current version with `status LinearResponse` in the julia package manager.

3. **Import Package:** Import the package in your julia environment with `import LinearResponse`.


---
## A first example

To reproduce the Plummer radial orbit instability calculation, see the example script in `examples/PlummerE/runExamplePlummer.jl`. This script will compute the location of the unstable radial orbit instability mode, using a simplified version of the calculation from [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract) (`n1max=1` instead of `n1max=10`, which results in a factor of 10 speedup). The outputs will all be cached, so re-running the example is inexpensive. This script will take approximately one minute to run.

Run the first example code with the following command:
```
$ julia examples/PlummerE/runExamplePlummer.jl
```

This example will first install some required libraries (`Plots`). These installations might take a few minutes when first called.

The resulting plot will be created with the name `examples/PlummerE/ROIdeterminant.png`.

![`Plummer ROI demonstration`](examples/PlummerE/ROIdeterminant.png)

In this image and test, we are highlighting two key results:

1. **Measurement of an unstable mode**, which is the pole located at $\omega=0.0+0.043i$. This mode location is verified against $N$-body simulations in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract).

2. **False poles in the lower half-plane**. All poles in the lower half-plane for this model are false poles owing to approximation of the function for linear response. We include them in this example as a caution against interpreting poles without validating via convergence tests.

An extension of this script, where the growth rate is computed for a range of radial anisotropy parameters, is Figure 1 in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract).

### Other examples

The `examples` directory also includes several other basic calculations for the Isochrone, Plummer, and Zang disc models. For example, to (very nearly) reproduce the Isochrone damped dipole mode calculation from Fouvry & Prunet (2022), see the example driver script in `examples/IsochroneE/runlinearresponseIsochroneISO.jl`. Note that this mode is not converged: adjusting parameters (in particular `n1max`) will result in different pole locations. The example `examples/IsochroneA/runlinearresponseIsochroneISO.jl` is precisely the calculation performed by [Fouvry & Prunet (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract) -- that is, with the analytic isochrone relations. This example will take substantially more computational effort, and can take up to twenty minutes to complete.

---

## Authors

Mathieu Roule - @MathieuRoule - roule@iap.fr

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

---

[^1]: For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

[^2]:If you are new to `julia`, install the latest version by running this in your terminal: `$ curl -fsSL https://install.julialang.org | sh`. If you are on Windows or run into problems with `curl`-based installation, please visit [this website](https://julialang.org/downloads/).
