# CPM_RD_matlab
A general puprose tool for modelling reaction-diffusion systems in MATLAB (with extensible support for other computational environments).

The aim of this project is to provide a tool that separates model specification from model analysis, simulation building, simulation analysis.

Currently, we have developped a simple text-based [model specification paradigm] (#model) for reaction-diffusion systems that is flxeible to be extended to consider other types of systems (e.g., reaction-diffusion-convection).  


*N.B. All code is written from MATLAB 2018+ (i.e., array braodcasting is used), but you may be able to run it on older versions with some elbow grease.*


## <a name='model'></a>Model Specification

# CPM_RD_matlab

A general purpose 2D stochastic implementation of reaction-diffusion systems using the SSSA (SSA0_fun.m), with a Cellular-Potts Model component (CPM_step.m) to provide a dyanmic domain for the chemical reactions and diffusion to occur in. See CPM_ReadMe.pdf for more details on these algorithms.



## Auto-generation of SSA Code (mk_fun)
Due to the way MATLAB handles functions (pass by value), it is computationally inefficient to compute chemical reaction propensities inside of a function. Instead, we specify the update rules in a script `update_alpha_chem0.m`, and automatically generate a SSSA function with those rules inside of it. 

The template for the SSSA is kept in a script called `SSA0.m`. This template is a fully functional MATLAB script, but it is written in a "covenient" manner and is not the most efficient implementation as it is a script calling scripts. Therefore, the template/script can be converted into a function automatically using mk_fun('SSA0'). This will create a new/overwrite the file `SSA0_fun.m` where any scripts called inside `SSA0.m` will be inlined for computational efficiency.

The generation of the function is fast enough that it can be run everytime a simulation is started. See line 115 of main2.m

## Chemical Reactions
Currently, the chemical reactions are based on the model found in [1] (see Tang2018.pdf). However, we aim to make the chemical model modular so that it can be modified to arbitrary chemical models without touching the SSSA algorithm files.

Currently, all propensities are computed from a single script `update_alpha_chem0.m`, which is a first step in the modularization. However, updates to the simulation state array are still hardcoded in `SSA0.m`, so there is still a fair amount of work to do on that front.




## References:


[1] Kaixi Tang,Colton G. Boudreau, Claire M. Brown ,Anmar Khadra. Paxillin phosphorylation at serine 273 and its effects on Rac, Rho and adhesion dynamic. PLOS Comp. Bio. 2018.
