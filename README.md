# BlackHolesDarkDress

Code for studying primordial black hole (PBH) binaries, clothed in dark matter (DM) halos. We use N-body simulations (along with analytical estimates) to follow the evolution of PBH binaries formed in the early Universe. *Movies are available [here](movies/) (archived on [FigShare](https://doi.org/10.6084/m9.figshare.6298397)).*

<img src="movies/snapshot1.png" alt="alt text" width=275> <img src="movies/snapshot2.png" alt="alt text" width=275> <img src="movies/snapshot3.png" alt="alt text" width=275>

If you have any questions, comments, bug-reports etc., please contact Bradley Kavanagh at bradkav@gmail.com. 

**You can run `StripAll.sh` before commiting any changes, which will clear all output from `.ipynb` notebooks. This should make things play nicer with git.**


### To-do

* Move notebooks to their own folder

## Contents

Folders:

* `Nbody/` - Code for setting up the Gadget-2 simulations and analysing the results. Note that the snapshot files are stored on FigShare - [https://doi.org/10.6084/m9.figshare.6300110](https://doi.org/10.6084/m9.figshare.6300110) - and there are a couple of scripts in the `Nbody/` folder for downloading and unpacking them.
* `plots/` - Plots associated with the code and paper. Some of the plots appearing here, do not feature in the paper, but we provide them for extra information.
* `movies/` - Animations of selected N-body simulations.
* `data/` - ...


Summary of notebooks (so far):

* `RemappingExamples.ipynb` - Example of how the remapping works, including a plot of the final semi-major axis as a function of initial semi-major axis
* `LIGO_sensitivities.ipynb` - Notebook for estimating the LIGO sensitivity curves $S(z)$




* `LIGO_limits.ipynb` - Notebook calculating the upper limit on the merger rate from LIGO, for a range of BH masses. The key output results can be found in `data/LIGO_R90.txt`.

## Requirements

The sampling routines, for generating samples of binaries, require [`emcee`](http://dfm.io/emcee/current/).

## To be added

* Some of the notebooks for analysing the N-body simulation results are currently only compatible with the low-resolution simulations. This will be fixed soon. In addition, there will also be some extra notebooks for checking the stability of the DM halos in single-PBH runs.
* Script for plotting the actual LIGO limits!

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
