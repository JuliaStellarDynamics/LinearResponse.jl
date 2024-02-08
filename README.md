
# LinearResponse.jl

[![image](https://github.com/JuliaStellarDynamics/LinearResponse.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/LinearResponse.jl/)

**LinearResponse.jl** is a package written in Julia to perform Linear Response calculations. The software is described in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract). To see an example of the results from the software, see [this repository](https://github.com/JuliaStellarDynamics/LinearResponse-paper/tree/main), which reproduces all the figures from Petersen et al. (2024).

---
## Installation

To fully use **LinearResponse.jl**, you will need **OrbitalElements.jl**, **FiniteHilbertTransform.jl**, and **AstroBasis.jl**. This installation assumes that you have not installed any of the four libraries; if you have already installed any of the others you may skip that portion of the installation.

The libraries under the JuliaStellarDynamics organisation are currently unregistered[^1]. To add them to your julia[^2] registry, follow these steps:

1. **Add Packages:** Use the package manager and execute the following command inside julia:
    ```julia
    add "git@github.com:JuliaStellarDynamics/OrbitalElements.jl.git"
    add "git@github.com:JuliaStellarDynamics/AstroBasis.jl.git"
    add "git@github.com:JuliaStellarDynamics/FiniteHilbertTransform.jl.git"
    add "git@github.com:JuliaStellarDynamics/LinearResponse.jl.git"
    ```
or at the command line
    ```$ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/OrbitalElements.jl.git")'
    $ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/AstroBasis.jl.git")'
    $ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git")'
    $ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/LinearResponse.jl.git")'
    ```

2. **Verify Version:** Confirm the current version with `status LinearResponse` in the julia package manager.

3. **Import Package:** Import the package in your julia environment with `import LinearResponse`.

## Working from Source

Alternatively, work directly from the codebase:

1. **Activate Environment:** In the main directory of the package, enter the Julia environment using `julia`.

2. **Access Package Manager:** Inside the Julia environment, open the package manager with `]`.

3. **Activate Project:** Activate the project using `activate .`. For added safety, resolve dependencies using `resolve` to check for updates.

4. **Return to Julia Interpreter:** Exit the package manager with `[backspace]`. You are now equipped with the latest package version.

5. **Import Package:** Import the package by typing `using LinearResponse` in the Julia interpreter.

Alternately, you may clone the repository wherever you want and create a local environment (or project) by running:
```
$ git clone https://github.com/JuliaStellarDynamics/LinearResponse.jl.git
$ cd LinearResponse.jl
$ julia --project=. -e 'using Pkg; Pkg.precompile()'
```

Note: If you are using a new Julia interpreter, you might need to download additional packages. Use the following command:
```julia
using(Pkg)
Pkg.instantiate()
```

---
### A first example

To reproduce the Plummer radial orbit instability calculation, see the example driver script in `examples/PlummerE/runlinearresponsePlummerROI.jl`. This script will compute the location of the unstable radial orbit instability mode, using a simplified version of the calculation from PRFPT24 (`n1max=1` instead of `n1max=10`, which results in a factor of 10 speedup). The outputs will all be cached, so re-running the example is inexpensive. This script will take approximately one minute to run.

Run the first example code with the following command:
```
$ julia examples/PlummerE/runlinearresponsePlummerROI.jl
```

This example will first install some required libraries (`Plots`). These installations might take a few minutes when first called.

The resulting plot will be created with the name `examples/PlummerE/ROIdeterminant.png`.

![`Plummer ROI demonstration`](examples/PlummerE/ROIdeterminant.png)

An extension of this script, where the growth rate is computed for a range of radial anisotropy parameters, is Figure 1 in Petersen et al. (2024).

---
### A second example

To (very nearly) reproduce the Isochrone damped dipole mode calculation from Fouvry & Prunet (2022), see the example driver script in `examples/IsochroneE/runlinearresponseIsochroneISO.jl`. Note that this mode is not converged: adjusting parameters (in particular `n1max`) will result in different pole locations. The example `examples/IsochroneA/runlinearresponseIsochroneISO.jl` is precisely the calculation performed by Fouvry & Prunet (2022) -- that is, with the analytic isochrone relations. This example will take substantially more computational effort, and can take up to twenty minutes to complete.

---

### Authors

Mathieu Roule - @MathieuRoule - roule@iap.fr

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

---

[^1]: For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

[^2]:If you are new to `julia`, install the latest version by running this in your terminal: `$ curl -fsSL https://install.julialang.org | sh`. If you are on Windows or run into problems with `curl`-based installation, please visit [this website](https://julialang.org/downloads/).
