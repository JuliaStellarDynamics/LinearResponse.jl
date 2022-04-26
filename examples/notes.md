## How to find the radial orbit instability of Saha (1991)

We are looking for the radial orbit instability in the isochrone cluster by scanning along the Re[omega]=0 axis for increasing Im[omega]. The key resonance is (n1,n2) = (-1,2), so we'll pay special attention to that.

0. Decide on parameters by seting up resonance vectors with `Resonances.jl` and basis vectors with `Basis.jl`.
1. Construct W_{l,n1,n2,p,q}(u,v) matrices for all basis functions and all resonances using `WMat.jl` and `run_wmat.jl`. W matrices are saved in wmat/.
2. Compute G_{l,n1,n2,p,q}(u) functions (GFunc.jl is hardwired for ROI right now!) using `GFunc.jl` and `run_Gfunc.jl`. G functions are saved in gfunc/.
3. Compute AXi using `Xi.jl` and `run_Axi.jl`.
4. Specify frequencies to probe with the main wrapper `find_ROI.jl`.
