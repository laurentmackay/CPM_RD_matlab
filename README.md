# CPM_RD_matlab

A general purpose 2D stochastic implementation of reaction-diffusion systems using the SSSA (CPM_chem_func.m), with a Cellular-Potts Model component (CPM_step.m) to provide a dyanmic domain for the chemical reactions and diffusion to occur in. See CPM_ReadMe.pdf for more details on these algorithms.

N.B. All code is written from MATLAB 2018+ (i.e., array braodcasting is used), but you may be able to run it on older versions with some elbow grease. 

## Compatability with MEX 
The SSSA component (CPM_chem_func.m) is mex-able so it can be converted to high-performance code for computationally-intensive scenarios (open CPM_chem_func.prj in coder and use main.m to autdefine variable types).

## Chemical Reactions
Currently, the chemical reactions are based on the model found in [1] (see Tang2018.pdf). However, we aim to make the chemical model modular so that it can be modified to arbitrary chemical models without touching the SSSA algorithm files.




## References:


[1] Kaixi Tang,Colton G. Boudreau, Claire M. Brown ,Anmar Khadra. Paxillin phosphorylation at serine 273 and its effects on Rac, Rho and adhesion dynamic. PLOS Comp. Bio. 2018.
