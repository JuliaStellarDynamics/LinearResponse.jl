
# LinearResponse.jl
## Version 0.9 (missing documentation)

`LinearResponse.jl` is a package written in Julia to perform Linear Response calculations. The software is described in Petersen, Roule, Fouvry, Pichon, Tep (2024; PRFPT24). To see an example of the results from the software, see the repository https://github.com/michael-petersen/LinearResponse-paper/tree/main, which reproduces all the figures from PRFPT24.

-----------------------------
### Quick installation and activation

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`dev .`). To be extra safe, you can `resolve` to check for updates[^1]. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import the exports by typing `using LinearResponse` into the Julia interpreter[^2].

Note that to fully use this library, you will also require `OrbitalElements.jl`, `AstroBasis.jl`, and `FiniteHilbertTransform.jl`.

-----------------------------
### A first example

To reproduce the Plummer radial orbit instability calculation, see the example driver script in `examples/PlummerE/runlinearresponsePlummerROI.jl`. This script will compute the location of the unstable radial orbit instability mode, using a simplified version of the calculation from PRFPT24 (`n1max=1` instead of `n1max=10`, which results in a factor of 10 speedup). The outputs will all be cached, so re-running the example is inexpensive. An extension of this script, where the growth rate is computed for a range of radial anisotropy parameters, is Figure 1 in PRFPT24.

-----------------------------
### A second example

To (very nearly) reproduce the Isochrone damped dipole mode calculation from Fouvry & Prunet (2022), see the example driver script in `examples/IsochroneE/runlinearresponseIsochroneISO.jl`. Note that this mode is not converged: adjusting parameters (in particular `n1max`) will result in different pole locations. The example `examples/IsochroneA/runlinearresponseIsochroneISO.jl` is precisely the calculation performed by Fouvry & Prunet (2022) -- that is, with the analytic isochrone relations.

----------------------------

### Authors

Mathieu Roule - @MathieuRoule - roule@iap.fr

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

----------------------------

### Development notes

Function arguments: should have the evaluation locations listed first (i.e. radius, eccentricity), followed by model parameters, followed by distribution function parameters. The last argument is always the `ResponseParameter` structure.

----------------------------

[^1]: As `LinearResponse` is (currently) unregistered, if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Short version: when in the package manager, `add "git@github.com:michael-petersen/JuliaLinearResponse.git"`. If you are getting an error about git keys, you will need to register your private key using the julia shell prompt (access with `;`), and then pointing at your private key: `ssh-add ~/.ssh/id_rsa`.

[^2]: You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`. If you want to access specific elements listed below, I recommend `import LinearResponse`.
