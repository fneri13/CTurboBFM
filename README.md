# CTurboBFM #
BFM solver for turbomachinery flows. It couples a 3D finite volume Euler solver with body force models needed
to mimic the blades effects on the flow field


### What is this repository for? ###

* Calculation of inviscid fluid dynamics problems (2D, Axisymmetric and 3D). Viscous and turbulent fluxes currently on development.
* Simulation of turbomachinery flows through the concept of body force models (axisymmetric and full-annulus grids)
* Complete description of the capabilities in the associated article (pending, under review).

### How do I get set up? ###

* Build system based on Makefile, tested for Mac OS. Three commands (make clean, make debug, make release).

* The Python [Unsflow/grid](link_pending) package is needed to generate grid files. 

### Requirements
- C++17 compiler (GCC/Clang)
- Python 3.x
- GoogleTest (for unit tests)
- Make
- [Unsflow/grid](link_pending)


### Running a test case

```bash
turbobfm input.ini
```

### Project structure

grids/       → mesh generation scripts

include/    → headers

python/     → post-processing scripts

src/        → solver and physics

test/       → unit & regression tests






### Notes ###
* The code has been written for Mac OS systems, so there can be portability issues for Windows machines.

* The system of Euler equations is solved with explicit time-stepping methods (3rd or 4th order Runge-Kutta). 
The time-step must therefore be accurately restricted below certain limits, and a large number of iterations may be required.

* Three possible time-step methods are available. Local for steady state analysis, globally CFL limited time-marching, and global fixed time step for unsteady simulations.

* Two advection schemes available (JST and Roe). Extension of Roe to second order with MUSCL reconstruction and flux limiters available.

* Three domain topologies available (2D, 3D, axisymmetric). 

* The following BFMs are available: Hall [1], Hall-Thollet [2], Chima [3], Correlations (under dev.), Lift/Drag [2] (under dev.).

* Possibility to restart full-annulus turbo simulations with the axisymmetric solution (meridional grids must coincide).

### Contribution guidelines ###

* Validate the modifications and update the test folder. The test folder contains unit and regression tests.

### Authors and contacts ###

- **Francesco Neri**, TU Delft, `f.neri@tudelft.nl`
- **Matteo Pini**, TU Delft, `m.pini@tudelft.nl`


### References ###

[1] Hall, David Kenneth. Analysis of civil aircraft propulsors with boundary layer ingestion. Diss. Massachusetts Institute of Technology, 2015.

[2] Thollet, William. "Body force modeling of fan-airframe interactions." ISAE-SUPAERO: Toulouse, France (2017).

[3] Chima, Rodrick V. A three-dimensional unsteady CFD model of compressor stability. Vol. 4241. 2006.

[4] Tom-Robin Teschner, [CFD University](https://cfd.university/)