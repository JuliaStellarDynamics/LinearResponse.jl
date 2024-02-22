
# LinearResponse.jl

[![image](https://github.com/JuliaStellarDynamics/LinearResponse.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/LinearResponse.jl/)

**LinearResponse.jl** is a package written in Julia to perform Linear Response calculations. The software is described in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract). To see an example of the results from the software, see [this repository](https://github.com/michael-petersen/LinearResponse-paper/tree/main), which reproduces all the figures from Petersen et al. (2024).

---
## Installation

To fully use **LinearResponse.jl**, you will need **OrbitalElements.jl**, **FiniteHilbertTransform.jl**, and **AstroBasis.jl**. This installation assumes that you have not installed any of the four libraries; if you have already installed any of the others (globally, not locally) you may skip that portion of the installation.

The libraries under the JuliaStellarDynamics organisation are currently unregistered[^1]. To add them to your julia[^2] registry, follow these steps:

You may add all the packages[^3] at once with this command:
    ```
    julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git");
    Pkg.add(url="https://github.com/JuliaStellarDynamics/OrbitalElements.jl.git");
    Pkg.add(url="https://github.com/JuliaStellarDynamics/AstroBasis.jl.git");
    Pkg.add(url="https://github.com/JuliaStellarDynamics/LinearResponse.jl.git")'
    ```

You can confirm the current version with `status LinearResponse` in the julia package manager.

---
## Quickstart

To reproduce the Plummer radial orbit instability calculation, see the example script in `examples/PlummerE/runExamplePlummer.jl`. Download the file by running:
```
wget https://github.com/JuliaStellarDynamics/LinearResponse.jl/blob/main/examples/PlummerE/runExamplePlummer.jl
```
This script will compute the location of the unstable radial orbit instability mode, using a simplified version of the calculation from [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract) (`n1max=1` instead of `n1max=10`, which results in a factor of 10 speedup). The outputs will all be cached (appearing as several files with the `.h5` extension in the folder where the script is run), so re-running the example is inexpensive. This script will take approximately one minute to run.

Run the first example code with the following command:
```
$ julia path/to/runExamplePlummer.jl
```

This example will first install some required libraries (`Plots`). These installations might take a few minutes when first called.

The resulting plot will be created with the name `ROIdeterminant.png`.

![`Plummer ROI demonstration`](examples/PlummerE/ROIdeterminant.png)

In this image and test, we are highlighting two key results:

1. **Measurement of an unstable mode**, which is the pole located at $\omega=0.0+0.043i$. This mode location is verified against $N$-body simulations in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract).

2. **False poles in the lower half-plane**. All poles in the lower half-plane for this model are false poles owing to approximation of the function for linear response. We include them in this example as a caution against interpreting poles without validating via convergence tests.

An extension of this script, where the growth rate is computed for a range of radial anisotropy parameters, is Figure 1 in [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract).

---
### Other examples

The `examples` directory also includes several other basic calculations for the Isochrone, Plummer, and Zang disc models. For example, to (very nearly) reproduce the Isochrone damped dipole mode calculation from Fouvry & Prunet (2022), see the example driver script in `examples/IsochroneE/runlinearresponseIsochroneISO.jl`. Note that this mode is not converged: adjusting parameters (in particular `n1max`) will result in different pole locations. The example `examples/IsochroneA/runlinearresponseIsochroneISO.jl` is precisely the calculation performed by [Fouvry & Prunet (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract) -- that is, with the analytic isochrone relations. This example will take substantially more computational effort, and can take up to twenty minutes to complete.

---
### Interactive notebook

If you prefer interactive Jupyter notebooks, you will need to install `IJulia` following these [instructions](https://julialang.github.io/IJulia.jl/stable/manual/installation/).

The interactive introduction example is then given in `examples/PlummerE/runExamplePlummer.ipynb`.


---
## Uninstall


First start by removing the packages from the environment by running
```
julia -e 'using Pkg; Pkg.rm("OrbitalElements"); Pkg.rm("AstroBasis");Pkg.rm("FiniteHilbertTransform");Pkg.rm("LinearResponse");'
```

If you worked in a test environment (that you do not want to keep) you can also simply erase the folder using `rm -r /path/to/my_env`.

Then to fully erase the package (installed in ~/.julia), run
```
julia -e 'using Pkg; using Dates; Pkg.gc(collect_delay=Day(0));'
```

It will erase all the packages which are not known in any of your "active" (i.e., for which the Manifest.toml file is reachable) project/environments, in particular `LinearResponse`.

---

## Authors

Mike Petersen -  [@michael-petersen](https://github.com/michael-petersen) - michael.petersen@roe.ac.uk

Mathieu Roule -  [@MathieuRoule](https://github.com/MathieuRoule) - roule@iap.fr

---

[^1]: For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

[^2]: If you are new to `julia`, install by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/). To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

[^3]: **Note on working with environments.** By default packages are added to the default environment at ~/.julia/environments/v1.#. It is however easy to create other, independent, projects. If you want to install the packages in a different/test environment, first create a folder to host the environment files (Project.toml and Manifest.toml which will be created later on). Then, for every command line invoking Julia, use `julia --project=/path/to/my_env` instead of `julia` alone. Note that packages will always be cloned in ~/.julia/packages but only accessible in your project's context.
