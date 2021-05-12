


## Model Specification
We have implemented a basic model specification scheme for reaction-diffusion systems. This model specification scheme serves as a basic template that can be expanded upon to model more complex types of systems.

Model are specified in simple line-oriented text files, where we are agnostic about file extensions. By line-oriented, we mean that declarations (e.g., a reaction declaration, variable declaration, or diffusion coefficient declaration) are separated from one another by newline characters. There is currently no limit on line lengths in the model specification files.

### Index
- [Reaction Notation](#reaction-notation)
  - [Sources and Sinks](#sources-and-sinks)
  - [Fast Reactions](#fast-reactions)
- [Variable Definitions](#variable-definitions)
- [Parameter Values](#parameter-values)
- [Initial Conditions](#initial-conditions)
- [Diffusion Notation](#diffusion-notation)
- [Note on Units](#units)

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

#### Sources and Sinks
Production and destruction of matter can be specified using the empty set or zero (i.e., `{}` or `0`) as the sole reactant or product of a reaction, respectively. For example, the production of `A` with a constant rate `v` can be written as
```
{} -> A; v
```
or 
```
0 -> A; v
```
where destruction can also be specified by simply reversing the direction of the arrow or swapping the reactant with the product.

#### Fast Reactions

Reactions set to quasi-steady state can be specified using a reversible reaction with a single rate constant. This rate constant corresponds to the the association constant of the reaction (i.e., `forward_rate/backward_rate`).

### Diffusion Notation
The diffusion coefficient of a chemical species `X` can be specified using
```
D(X)=1
```
where we have specified a diffusion coefficient of 1 (space units^2 / time units, see note on [units](#units)). For convenience, when multiple species have the same diffusion coefficient, their diffusion coefficients may all be specified using a compound declaration as follows:
```
D(X)=D(Y)=D(Z)=1
```

### Variable Definitions
While we have used the name of a rate "constant", the value of the rate constants specified above may in fact be functions of the chemical species (e.g., when QSSA has already been applied to the model). In such a case, one may define the rate constants to be variable quantities by using simple algebraic expressions.

For example, we may model mutual inhibition between two molecules (`A` and `B`) that can switch betweem two isoforms (e.g., through reversible isomerization reactions `A0<->A1` and `B0<->B1`) using the following model specification:
```
A0 <-> A1; kA, delta_A
B0 <-> B1; kB, delta_B

kA = 1/(1+B1^2)
kB = 1/(1+A1^2)
```

Variable definitions may reference other variable definitions, but note that any variable name referenced in a variable definition must have already been defined in the file. That is, we will not sort out the order of variable definitions for the user, nor will we notify the user when a model is defined in an inconsistent manner.

### Parameter Values
We have borrowed heavily from xppaut for our parameter specification paradigm. Any line starting with "p" (e.g., "par" or "param") denotes a line of parameter value declarations. For example, in the isomerization example considered above we could specify parameters values using:
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

### Initial Conditions
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



# Units
Model specification is unit-agnostic and it is up to a specific simulation to interpret the numerical values specified by the user apprpriately.
