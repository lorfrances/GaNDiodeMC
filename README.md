# GaNDiodeMC
The following Julia source code was developed for the ensemble Monte Carlo (MC) Boltzmann transport simulation of electrons in a GaN-based diode device. It can simulate a thin heterobarrier layer intended to absorb more phonons and increase the electron current. Results generated from this simulation have been included in the following manuscript:

> Franceschetti, L., Kaviany, M., and Shin, S., “Heterobarrier in-situ phonon recycling in semiconductor diodes,” under review.

The MC simulation features include: scattering with phonons and impurties, spherical-parabolic Γ band electrons, Runge-Kutta 2 equation of motion integration, self-consistent Poisson equation, among others.

## Structure
The simulation is initialized from _MonteCarlo.jl_, which calls the primary MC temporal loop in _MainLoops.jl_ (executes one timestep). Within that timestep, every particle is drifted and potentially scattered. At the end of the timestep, the Poisson equation is evaluated, and data is collected. This continues until the end of the simulation with the specified duration.
## Usage
At this time, a sample input file (INPUT.txt) has been included with a few customizable parameters. The program is run by executing _MonteCarlo.jl_from the command line.

More details coming soon...
