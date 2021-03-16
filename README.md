# CPM_RD_matlab
A general puprose tool for modelling reaction-diffusion systems in MATLAB (with extensible support for other computational environments). Inspired by the very fun and useful dynamical systems analysis tool [xppaut](http://www.math.pitt.edu/~bard/xpp/xpp.html), we have aimed to make chemical kinetic systems more approachable by providing a tool that automatically derives sets of equations from reactions schemes and can perform some of the more tedious aspects of model simplification in an reliable and automated manner.

A longer term aim of this project is to provide a tool that separates model specification from model analysis, simulation building, simulation analysis. Moreover, due to our reaction-based model specification format, combining models can be made as simple as copy-pasting the contents of two files.

Currently, we have developped a simple text-based [model specification paradigm](#model-specification) for reaction-diffusion systems that is flxeible to be extended to consider other types of systems (e.g., reaction-diffusion-convection).  


*N.B. All code is written from MATLAB 2018+ (i.e., array braodcasting is used), but one may be able to run it on older versions with some elbow grease.*

### GoTo
- [Model Specification](#model-specification)
 - [Reaction Notation](#reaction-notation)
 - [Variable Definitions](#variable-definitions)
 - [Diffusion Notation](#diffusion-notation)


## Model Specification
We have implemented a basic model specification scheme for reaction-diffusion systems. This model specification scheme serves as a basic template that can be expanded upon to model more complex types of systems.

Model are specified in simple line-oriented text files, where we are agnostic about file extensions. By line-oriented, we mean that declarations (e.g., a reaction declaration, variable declaration, or diffusion coefficient declaration) are separated from one another by newline characters. There is currently no limit on line lengths in the model specification files.


### Reaction Notation
Reaction-diffusion systems are specified primarily using a text-based chemical reaction notation. For, example a bimolecular complexing reaction between chemical species `A` and `B` producing `C` with rate constant `k1` is written as:

```
A + B -> C; k1
```
Reversible pairs of reactions may be specified using two separate declarations as above (with products and reactants swapped), or one may use the shorthand notation,
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
depending on whether one wishes the reaction to be reversible or not. In the absence of an explicit stoichiometric coefficient, a value of 1 is assumed.

#### Variable Definitions
While we have used the name of a rate "constant", the value of the rate constants specified above may in fact be functions of the chemical species (e.g., when QSSA has already been applied to the model). In such a case, one may define the rate constants to be variable quantities by using simple algebraic expressions.

For example, we may model mutual inhibition between two molecules (`A` and `B`) that can switch betweem two isoforms (e.g., through reversible isomerization reactions `A0<->A1` and `B0<->B1`) using the following model specification:
```
A0 <-> A1; kA, delta_A
B0 <-> B1; kB, delta_B

kA = 1/(1+B^2)
kB = 1/(1+A^2)
```

Variable definitions may reference other variable definitions, but note that any variable name referenced in a variable definition must have already been defined in the file. That is, we will not sort out the order of variable definitions for the user, nor will we notify the user when a model is defined in an inconsistent manner.

#### Parameter Values
We have borrowed heavily from xppaut for our parameter specification paradigm. Any line starting with "par" or "param" denotes a line of parameter value declarations. For example, in the isomerization example considered above we could specify parameters values using:
```
par delta_A=0.1, delta_B=2
```

Moreover, we have added in  wildcard approach for specifying a default parameter value that is used for parameters that are unspecified. For axmple, we can specify a default parameter value of 1 using:
```
par *=1
```
Unlike in xppaut, we allow parameter values to reference one another. However, similarly to variable declarations, they must be specified in appropraiate order. For example:
```
par delta_A = 0.1, delta_B = 1 + 10*delta_A 
```

#### Initial Conditions
Initial condtions for the molecule `X` are specified using by placing `(0)=` after the chemical species name and either specifying a numerical value of paremeter name. For example:
```
X(0)=1
```
or 
```
X(0)=X0
par X0=1
```
where in both cases we have specified an initial condition of 1 (concentration units, see note on [units](#units)). For convenience, multiple species can be given the same initial condition using a compound declaration as follows:
```
X(0)=Y(0)=Z(0)=1
```
We also allow for wildcard notation to define a default intial conditions:
```
*(0)=1
```

### Diffusion Notation
The diffusion coefficient of a chemical species `X` can be specified using
```
D(X)=1
```
where we have specified a diffusion coefficient of 1 (space units^2 / time units, see note on [units](#units)). For convenience, when multiple species have the same diffusion coefficient, their diffusion coefficients may all be specified using a compound declaration as follows:
```
D(X)=D(Y)=D(Z)=1
```
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


