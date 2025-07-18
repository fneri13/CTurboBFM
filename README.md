# CTurboBFM #
C++ Euler solver.


### What is this repository for? ###

* Resolution of inviscid fluid dynamics problems (2D, Axisymmetric and 3D)
* Simulation of turbomachinery flows through the concept of body force models (axisymmetric and full-annulus)


### How do I get set up? ###

* Build system based on Makefile tested for Mac OS and UNIX systems. Three commands (make clean, make debug, make release).


### Results Example ###

##### Schulz benchmark test case 3 #####
Test case number 3 for shock-capturing properties of numerical schemes, taken from Schulz Riemann problems.
The following picture shows the magnitude of the density gradient obtained with Roe scheme, 2nd MUSCL reconstruction and Van Leer limiter on a 400x400 grid:
![Roe + MUSCL + Van Leer](images/schulz_roe2.png)



### Notes ###
* The code has been written for Mac OS systems, so there are likely of portability issues for windows machines.

* The system of Euler equations is solved with explicit time-stepping methods (3rd or 4th order Runge-Kutta). This means that the time-step must be accurately restricted below certain limits, and a large number of iterations may be required to simulate a certain problem.

* Three possible time-step methods. Local for steady state analysis, global for CFL limited global time step, and fixed time step for unsteady simulations.

* Two convection schemes available (JST and Roe). Extension of Roe to second order with MUSCL reconstruction and flux limiters available.

* Three domain topologies available (2D, 3D, axisymmetric). 

* Possibility to restart full-annulus turbo simulations with the axisymmetric solution (meridional grids must coincide).







### Contribution guidelines ###

* Validate the modifications by means of detailed test cases (use google gtest)

### Authors and contacts ###

- **Francesco Neri**, TU Delft, `f.neri@tudelft.nl`