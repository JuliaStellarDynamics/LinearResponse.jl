# Installation

To fully use **LinearResponse.jl**, you will need **OrbitalElements.jl**, **FiniteHilbertTransform.jl**, and **AstroBasis.jl**. This installation assumes that you have not installed any of the four libraries; if you have already installed any of the others (globally, not locally) you may skip that portion of the installation.

The libraries under the JuliaStellarDynamics organisation are currently unregistered[^1]. To add them to your julia[^2] registry, follow these steps:

You may add all the packages[^3][^4] at once with this command:
    ```
    julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/OrbitalElements.jl.git"); Pkg.add(url="https://github.com/JuliaStellarDynamics/AstroBasis.jl.git"); Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git"); Pkg.add(url="https://github.com/JuliaStellarDynamics/LinearResponse.jl.git")'
    ```

You can confirm the current version with `status LinearResponse` in the julia package manager.

[^1]: For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

[^2]: If you are new to `julia`, install by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/). To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

[^3]: **Note on working with environments.** By default packages are added to the default environment at ~/.julia/environments/v1.#. It is however easy to create other, independent, projects. If you want to install the `LinearResponse` package in a different/test environment, first create a folder to host the environment files (Project.toml and Manifest.toml which will be created later on). Then, for every command line invoking Julia, use `julia --project=/path/to/my_env` instead of `julia` alone. Note that packages will always be cloned in ~/.julia/packages but only accessible in your project's context.

[^4]: You may also uninstall the libraries by removing the packages, by running
```
julia -e 'using Pkg; Pkg.rm("OrbitalElements"); Pkg.rm("AstroBasis");Pkg.rm("FiniteHilbertTransform");Pkg.rm("LinearResponse");'
```