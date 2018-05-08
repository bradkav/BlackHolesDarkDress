# PBHdress

If you have any questions, comments, bug-reports etc., please contact Bradley Kavanagh at bradkav@gmail.com. 

**You can run `StripAll.sh` before commiting any changes, which will clear all output from `.ipynb` notebooks. This should make things play nicer with git.**


### To-do

* Move notebooks to their own folder

## Contents

Summary of notebooks (so far):

* `RemappingExamples.ipynb` - Example of how the remapping works, including a plot of the final semi-major axis as a function of initial semi-major axis
* `PBH-sampling.ipynb` - Very rough example notebook for how the sampling procedure (using emcee) works
* `LIGO_limits.ipynb` - Notebook calculating the upper limit on the merger rate from LIGO, for a range of BH masses. The key output results can be found in `data/LIGO_R90.txt`.


## Requirements

The sampling routines, for generating samples of binaries, require [`emcee`](http://dfm.io/emcee/current/).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
