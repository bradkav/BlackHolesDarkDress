# BlackHolesDarkDress - Simulations



### Requirements

Requires:

* `pyGadgetReader` - available [here](http://ascl.net/1411.001)
* `pyGadgetIC` - available [here](https://github.com/ldocao/pygadgetic)
* `sympy` - available [here](http://www.sympy.org/en/index.html)

### Downloading the results

To keep the size of this repo manageable, the `Gadget-2` snapshots are archived on FigShare: [doi:10.6084/m9.figshare.6300110](https://doi.org/10.6084/m9.figshare.6300110).

To download and unpack all of the N-body results into the `sims/` folder, simply run the script `GetAllSnapshots.sh`. But beware that the total size of the snapshots is around 25 GB and may take an hour or so to download. 

To download one example simulation (~1-2 GB), simply run`GetExampleSnapshot.sh`.

Alternatively you can visit [doi:10.6084/m9.figshare.6300110](https://doi.org/10.6084/m9.figshare.6300110) and download the individual simulations.


### Code

* `CalcOrbitalElements.ipynb` - Calculate the semi-major axis and eccentricity of the PBH binary system during the course of the simulations
* `PlotSimulationResults.ipynb` - Read in simulation results and plot some properties of the system - such as separation, enclosed DM halo mass and angular momentum.
* `NbodyResults.ipynb` - Notebook for plotting the results of the N-body simulations (final semi-major axis and eccentricity, etc).