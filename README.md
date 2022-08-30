
# CallAResponse.jl

`CallAResponse.jl` is a package written in Julia to perform Linear Response calculations.

-----------------------------
### Quick activate

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`). To be extra safe, you can `resolve` to check for updates. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import the exports by typing `using CallAResponse` into the Julia interpreter. You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`. If you want to access specific elements listed below, I recommend `import CallAResponse`.

As `CallAResponse` is (currently) unregistered, if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Short version: when in the package manager, `add "git@github.com:michael-petersen/JuliaCallAResponse.git"`. If you are getting an error about git keys, you will need to register your private key using the julia shell prompt (access with `;`), and then pointing at your private key: `ssh-add ~/.ssh/id_rsa`.

-----------------------------
### Constructing W matrices

-----------------------------
### Constructing G functions

-----------------------------
### Computing Xi coefficients

-----------------------------

### Authors

Mike Petersen -  @michael-petersen - petersen@iap.fr
Mathieu Roule - @MathieuRoule - roule@iap.fr
