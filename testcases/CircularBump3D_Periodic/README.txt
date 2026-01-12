Test case with different sectors to test periodic BCs enforcement.
The strategy used in the solver (vertex-centered) is to sum and average the residuals contribution 
due to advection fluxes and volumetric source terms of periodic node mates.

In principle this should avoid to have discontinuities at the coordinate-cut for 
full annulus simulations.
