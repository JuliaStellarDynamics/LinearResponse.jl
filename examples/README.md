## How to compute linear response

One only needs to modify two files in order to compute the linear response for a given model:
1. Design the model: potential, distribution function, and parameters, by building or modifying `ModelParam*.jl`.
2. Edit `runlinearresonse*.jl` to include the list of frequencies you would like to evaluate.

Several examples are provided:
1. IsochroneA/IsochroneE: We are looking for the radial orbit instability (ROI) in the isochrone cluster by scanning along the Re[omega]=0 axis for increasing Im[omega]. The key resonance is (n1,n2) = (1,-2)+(-1,2) so we'll pay special attention to those.
2. Plummer: We are approaching the ROI in a Plummer sphere.
3. MestelZang: A 2d linear response example.
