# CPM_RD_matlab
A general puprose tool for modelling reaction-diffusion systems in MATLAB (with extensible support for other computational environments).

The aim of this project is to provide a tool that separates model specification from model analysis, simulation building, simulation analysis.

Currently, we have developped a simple text-based [model specification paradigm](#model-specification) for reaction-diffusion systems that is flxeible to be extended to consider other types of systems (e.g., reaction-diffusion-convection).  


*N.B. All code is written from MATLAB 2018+ (i.e., array braodcasting is used), but you may be able to run it on older versions with some elbow grease.*


## Model Specification
### Reaction Notation
Reaction-Diffusion systems are specified primarily using a text-based chemical reaction notation. For, example a bimolecular complexing reaction between chemical species `A` and `B` producing `C` with rate constant `k1` is written as:

```
A + B -> C; k1
```
Reversible pairs of reactions may be specified using two separate declartions as above (with products and reactants swapped), or one may use the shorthand notation,
```
A + B <-> C; k1, k2
```
where `k2` is the rate constant for the degradation of `C` into `A` and `B`.

Most generally, reactions are described by stoichiometirc coeffient and chemical species. Consider a reaction whose `n`th chemical species is denoted by `Xn` and has a stoichiometric coefficent of `rn` as a reactant and `pn` as a product. Such a reaction can be written as  
```
r1*X1 + r2*X2 + ... + rN*XN -> p1*X1 + p2*X2 + ... + pN*XN; forward_rate
```
or
```
r1*X1 + r2*X2 + ... + rN*XN <-> p1*X1 + p2*X2 + ... + pN*XN; forward_rate , backward_rate
```
depending on whether one wishes the reaction to be reversible or not.

#### Variable definitions
While we have used the name of a rate "constant", the value of the rate constants specified above may in fact be functions of the chemical species (e.g., when QSSA has already been applied to the model). In such a case, one may define the rate constants to be variable quantities by using simple algebraic expressions.

For example, we may model mutual inhibition between two molecules (`A` and `B`) that can switch betweem two isomers (e.g., through isomerization reactions `A0<->A1` and `B0<->B1`) using the following model specification:
```
A0 <-> A1; kA, delta_A
B0 <-> B1; kB, delta_B

kA = 1/(1+B^2)
kB = 1/(1+A^2)
```

Variable definitions may reference other variable definitions, but note that any variable name referenced in a variable definition must have already been defined in the file That is, we will not sort out the order of variable definitions for you.


