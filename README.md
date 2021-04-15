# CPM_RD_matlab
A general puprose tool for modelling reaction-diffusion systems in MATLAB (with extensible support for other computational environments). Inspired by the very fun and useful dynamical systems analysis tool [xppaut](http://www.math.pitt.edu/~bard/xpp/xpp.html), we have aimed to make chemical kinetic systems more approachable by providing a tool that automatically derives sets of equations from reactions schemes and can perform some of the more tedious aspects of model simplification in an reliable and automated manner. 

Moreover, due to our reaction-based model specification format, combining models is often as simple as copy-pasting the contents of two files.

A longer term aim of this project is to provide a tool that separates model specification from model analysis, simulation building, and simulation analysis. 

Currently, we have developped a simple text-based [model specification paradigm](MODELS.md) for reaction-diffusion systems that is flxeible to be extended to consider other types of systems (e.g., reaction-diffusion-convection).  


*N.B. All code is written from MATLAB 2018+ (i.e., array braodcasting is used), but one may be able to run it on older versions with some elbow grease.*

## Overview

This library provides three main components that are useful for CPM simulations:
1. An extensible [model specification format](MODELS.md) and a parser that will generate model-specific files for running CPM simulations.
2. A core set of functions/scripts for running CPM simulations.
3. A basic system for switching between models, running simulations, and storing results in a systematic manner.

### Index
- [Deploying Models](#deploying-models)
  - [Simulations](#simulations)
  - [Analysis](#analysis)


## Installation
The repository should be downloaded:
```
git clone https://github.com/laurentmackay/CPM_RD_matlab.git
cd CPM_RD_matlab
```

Furthermore, its base directory should be added to the MATLAB path (all subsequent modifications to the path are handled algorithmically).


## Testing
As a basic example, try the entering the following commands into the MATLAB command line:

```
deploy_model('chem_Rx_Pax_Asheesh')
main_FVM
```

Eventually, a window should pop up with CPM cell that moves around. Check the `_chem_Rx_Pax_Asheesh` folder for .ode and .f90 files that can be used with xppaut and AUTO, respectively.

You may also have a look at the script [main_FVM.m](/protocols/CPM/main_FVM.m), to see if there are any simulation parameters you wish to change.


## General Usage

After creating a model specification file (e.g., in the file `model_file`) for your model, one should:

1. Deploy the model
2. Create a "virtual experiment" to help organize results on the file system using `set_experiment()`. [Optional]
3. Run a CPM simulation using `main_FVM.m`
4. Analyze the results using tools in [/protocols/CPM/analysis/](/protocols/CPM/analysis/)

### Deploying Models
A model specified in the file `model_file` can be used in simulations by first calling `deploy_model('model_file')`. This will parse the model specified in `model_file` and generate model-dependent functions necessary for running simulations. These generated files are stored in a directory called `_model_file` and that directory is temporarily added to the MATLAB-path. 

Only one model can be deployed at a time, and `deploy_model` will take care of modifying the MATLAB-path when a user switches between models. One may check the currently active model using the `active_model` global variable.

Moreover, `deploy_model` will also check if the model specification file has been updated since it was last deployed, and will only run the model through the interpreter when changes have been made. This behaviour can be overridden by passing a second boolean argument `true` or `1` to `deploy_model`.



### Experiments

In order to prevent overwriting of files and disroganized data, we use "experiments" to define where results should be saved. In order to set the name of your current experiment, use the command `set_experiment('experiment_name')`.

The current results directory can be obtained using `results_dir()`, and, in general, is given by `<CPM_RD_matlab>/_model_name/results/experiment_name` where `<CPM_RD_matlab>` is the path to the base directory of this library.

The function `ls_results()` can be used to check the results directory for `.mat` files, or one may also use `ls_results(ext)` to find all files that end in `.ext`.


