# BlackHolesDarkDress - N-body Simulations

Here, we provide code for setting up, running and analysing the results of our N-body simulations using `Gadget-2` (which is available [here](https://wwwmpa.mpa-garching.mpg.de/gadget/)).

Note that at the moment many of the files have just been copied verbatim from the computing cluster where they were being used, so they may not run smoothly anywhere else (incorrect paths etc.). However, this will be cleared up soon, along with explanations of how to run the different codes/scripts.

### Requirements

Requires:

* `pyGadgetReader` - available [here](http://ascl.net/1411.001)
* `pyGadgetIC` - available [here](https://github.com/ldocao/pygadgetic)
* `sympy` - available [here](http://www.sympy.org/en/index.html)

### Downloading the results snapshot files

To keep the size of this repo manageable, the `Gadget-2` snapshots are archived on FigShare: [doi:10.6084/m9.figshare.6300110](https://doi.org/10.6084/m9.figshare.6300110).

To download and unpack all of the N-body results into the `sims/` folder, simply run the script `GetAllSnapshots.sh`. But beware that the total size of the snapshots is around 25 GB and may take an hour or so to download. 

To download one example simulation (~1-2 GB), for example to use some of the analysis notebooks, simply run `GetExampleSnapshot.sh`.

Alternatively you can visit [doi:10.6084/m9.figshare.6300110](https://doi.org/10.6084/m9.figshare.6300110) and download the individual simulations.


### Code

* `CalcOrbitalElements.ipynb` - Calculate the semi-major axis and eccentricity of the PBH binary system during the course of the simulations
* `PlotSimulationResults.ipynb` - Read in simulation results and plot some properties of the system - such as separation, enclosed DM halo mass and angular momentum.
* `NbodyResults.ipynb` - Notebook for plotting the results of the N-body simulations (final semi-major axis and eccentricity, etc).
* `DistributionFunction.ipynb` - Notebook for generating the distribution function for a given DM halo, following the Eddington inversion formula.
* `GenerateICs_*.py` - A number of scripts for generating the Gadget-2 snapshots for the initial conditions. These rely on `PBH.py` and `eddington.py` to sample from the DM halo distribution function. *More details on how to use these initial conditions generators will be added soon.*
*  `animate.py` and `animate_zoom.py` - Scripts for generating movies from snapshot files.


### Gadget-2 files

The folder `run/` contains example Gadget-2 parameter files, while `submit_scripts` contains submission scripts for setting up the relevant files and submitting to a computing cluster. These have not been cleaned up, so they're mostly there to give an idea of how the parameter files work (how the parameter values are set, etc.).
