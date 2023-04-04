
# LinearResponse.jl

`LinearResponse.jl` is a package written in Julia to perform Linear Response calculations.

-----------------------------
### Quick activate

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`dev .`). To be extra safe, you can `resolve` to check for updates[^1]. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import the exports by typing `using LinearResponse` into the Julia interpreter[^2].

-----------------------------
### Constructing W matrices

-----------------------------
### Constructing G functions

-----------------------------
### Computing Xi coefficients

-----------------------------

### Style notes

Function arguments: should have the evaluation locations listed first (i.e. radius, eccentricity), followed by model parameters, followed by distribution function parameters. The last argument is always the `ResponseParameter` structure.

----------------------------

### Authors

Mathieu Roule - @MathieuRoule - roule@iap.fr

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

----------------------------

[^1]: As `LinearResponse` is (currently) unregistered, if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Short version: when in the package manager, `add "git@github.com:michael-petersen/JuliaLinearResponse.git"`. If you are getting an error about git keys, you will need to register your private key using the julia shell prompt (access with `;`), and then pointing at your private key: `ssh-add ~/.ssh/id_rsa`.

[^2]: You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`. If you want to access specific elements listed below, I recommend `import LinearResponse`.
