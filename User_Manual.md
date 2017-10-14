
# Chemical Reaction Rate Simulator

## Model Document

*Copyrights © G5 CS207, team member Thomas Culp, Brianna McHorse, Yujiao Chen, Yi Zhai*

## Introduction:
***

Chemical Reaction Rate Simulator (CRRS) is a suite of software for the calculation of chemical kinetics. It can simulate a system consisting of $N$ species undergoing $M$ **irreversible**, **elementary** reactions of the form 

\begin{align}
  \sum_{i=1}^{N}{\nu_{ij}^{\prime}\mathcal{S}_{i}} \longrightarrow 
  \sum_{i=1}^{N}{\nu_{ij}^{\prime\prime}\mathcal{S}_{i}}, \qquad j = 1, \ldots, M
\end{align}

the rate of change of specie $i$ (the reaction rate) can be written as 
\begin{align}
  f_{i} = \sum_{j=1}^{M}{\nu_{ij}\omega_{j}}, \qquad i = 1, \ldots, N
\end{align}

where the progress rate for each reaction is given by 
\begin{align}
  \omega_{j} = k_{j}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime}}}, \qquad j = 1, \ldots, M
\end{align}

and $k_{j}$ is the forward reaction rate coefficient.

Nomenclature Table:

| Symbol | Meaning |
|:--------:|:-------:|
| $\mathcal{S}_{i}$ | Chemical symbol of specie $i$ |
| $\nu_{ij}^{\prime}$ | Stoichiometric coefficients of reactants |
| $\nu_{ij}^{\prime\prime}$ | Stoichiometric coefficients of products |
| $N$                       | Number of species in system |
| $M$                       | Number of elementary reactions |
| $f_{i}$                   | Rate of consumption or formation of specie $i$ (reaction rate) |
| $\omega_{j}$              | Progress rate of reaction $j$ |
| $x_{i}$                   | Concentration of specie $i$ |
| $k_{j}$                   | Reaction rate coefficient for reaction $j$ |

### The CRRS has the following features:

* The CRRS can be applied to irreversible, elementary chemical reactions, and is extensible to reversible reactions and non-elementary reactions.


* The CRRS can simulate chemical reactions using 

   * Constant reaction rate coefficients
   
    \begin{align}
      &k_{\textrm{const}}   = k 
    \end{align}
   
   * Arrhenius reaction rate coefficients
   
    \begin{align}
      &k_{\textrm{arr}}     = A \exp\left(-\frac{E}{RT}\right) 
    \end{align}
   * Modified Arrhenius (Kooij) reaction rate coefficients
   
    \begin{align}
      &k_{\textrm{mod arr}} = A T^{b} \exp\left(-\frac{E}{RT}\right) 
    \end{align}
   
   It is also extensible to other types of reaction rate coefficient.





## Installation:
***

1. Please visit 
    https://github.com/CS207-G5/cs207-FinalProject

2. Download package CRRS.zip (*to be updated*) 

3. Unzip CRRS.zip to your local hard drive

4. Run Interface.py (*to be updated*) to simulate your chemical reaction

The running environment:
* Python 3.6

Dependency packages:
* Sys, numpy, math, cmath, number, xml, pytest

## Basic usage and examples:
***
There are two ways to interact with CRRS:
* **_read in .xml file_**
* **_direct manual input_**

Either way is acceptable to CRRS interface.

#### Direct manual input examples:

**Example 1**

Compute the reaction rates temperature $T = 1500$ of the following system of irreversible reactions:
\begin{align}
  2H_{2} + O_{2} \longrightarrow 2OH + H_{2} \\
  OH + HO_{2} \longrightarrow H_{2}O + O_{2} \\
  H_{2}O + O_{2} \longrightarrow HO_{2} + OH
\end{align}

The reaction 1 is a modified Arrhenius reaction with $A_{1} = 10^{8}$, $b_{1} = 0.5$, $E_{1} = 5\times 10^{4}$, reaction 2 has a constant reaction rate parameter $k = 10^{4}$, and reaction 3 is an Arrhenius reaction with $A_{3} = 10^{7}$ and $E_{3} = 10^{4}$.

Given the following species concentrations:

\begin{align}
  \mathbf{x} = 
  \begin{bmatrix}
    H_{2}  \\
    O_{2}  \\
    OH     \\
    HO_{2} \\
    H_{2}O
  \end{bmatrix} = 
  \begin{bmatrix}
    2.0 \\
    1.0 \\
    0.5 \\
    1.0 \\
    1.0
  \end{bmatrix}
\end{align}

The direct interaction:
```python
import numpy as np
import chemkin
import reaction_coeffs as rrc

M = 3 # Number of reaction
N = 5 # Number of species
concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
Temps = 1500.0
A1, b1, E1 = 1.0e+08, 0.5, 5.0e+04
A2, E2 = 1.0e+07, 1.0e+04
nu_react = np.array([[2.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
nu_prod = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [2.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])    
k = np.zeros(M)
# Get the reaction rate coefficients
k[0] = rrc.k_mod_arr(A1, b1, E1, T)
k[1] = rrc.k_const(1.0e+04)
k[2] = rrc.k_arr(A2, E2, T)
# Get progress rates for reactions
omega = chemkin.progress_rate(nu_react, concs, k)
# Get reaction rates for species
f = chemkin.reaction_rate(nu_react, nu_prod, omega)
# Print reaction rates to screen
print(f)
```

The reaction rates of the species returned by CRRS are: 
```python
[ -2.81117621e+08  -2.85597559e+08   5.66715180e+08   4.47993847e+06   -4.47993847e+06]
```


```python

```
