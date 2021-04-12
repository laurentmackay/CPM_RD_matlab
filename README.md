# CPM_RD_matlab
A general puprose tool for modelling reaction-diffusion systems in MATLAB (with extensible support for other computational environments). Inspired by the very fun and useful dynamical systems analysis tool [xppaut](http://www.math.pitt.edu/~bard/xpp/xpp.html), we have aimed to make chemical kinetic systems more approachable by providing a tool that automatically derives sets of equations from reactions schemes and can perform some of the more tedious aspects of model simplification in an reliable and automated manner.

A longer term aim of this project is to provide a tool that separates model specification from model analysis, simulation building, simulation analysis. Moreover, due to our reaction-based model specification format, combining models can be made as simple as copy-pasting the contents of two files.

Currently, we have developped a simple text-based [model specification paradigm](#model-specification) for reaction-diffusion systems that is flxeible to be extended to consider other types of systems (e.g., reaction-diffusion-convection).  


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
- [Note on Units](#units)


## Deploying Models
A model specified in the file `model_file` can be used in simulations by first calling `deploy_model('model_file')`. This will run  the model specified in `model_file` through an interpreter (i.e., `mk_rxn_files`) and pass it along to simulation-dependent functions that will generate appropriate files necessary for running simulations. These output files are stored in a directory called `_model_file` and that directory is temporarily added to the MATLAB-path. 

Only one model can be deployed at a time, and `deploy_model` will take care of modifying the MATLAB-path when a user switches between models.

Moreover, `deploy_model` will also check if the model specification file has been updated since it was last deployed, and will only run the model through the interpreter when changes have been made. This behaviour can be overridden by passing a second boolean argument `true` or `1` to `deploy_model`.

## Simulations
Upon deploying a model, simulation specific files will be generated. Currently, we produce files that allow the user to simulate the model inside the Cellular Potts Model framework as well as some files for producing bifurcations in AUTO07p.

Moving forward, we will provide a more flexible approach such that arbitrary simulation files can be generated from user-specified model.

## Analysis

I am just starting to think about how to manage this within our system. However, we should probably distingush between model analysis and analaysis of simulation results as they are generally quite different, but often inter-related, beasts.


# Units
Model specification is unit-agnostic and it is up to a specific simulation to interpret the numerical values specified by the user apprpriately.


