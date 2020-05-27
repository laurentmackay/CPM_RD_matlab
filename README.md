# CPM_RD_matlab

A general purpose 2D stochastic implementation of reaction-diffusion systems using the SSSA, with a Cellular-Potts Model component to provide a dyanmic domain for the chemical reactions and diffusion to occur in. All code is written from MATLAB 2018+ (i.e., array braodcasting is used), and the SSSA component is mex-able so it can be converted to high-performance code for computationally-intensive scenarios.

Currently, the chemical reactions are based on the model found in [1]. However, we aim to make the chemical model modular so that it can be modified to arbitrary chemical models without touching the SSSA algorithm files.






References:
[1] Kaixi Tang,Colton G. Boudreau, Claire M. Brown ,Anmar Khadra. Paxillin phosphorylation at serine 273 and its effects on Rac, Rho and adhesion dynamic. PLOS Comp. Bio. 2018.
