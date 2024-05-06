This repository hosts code to implement the mathematical model and analysis for the (draft) paper:

*A global mathematical model of climatic suitability for Plasmodium falciparum malaria.*

by Owen Brown, Jen A Flegg, Daniel J Weiss, Nick Golding.

This code adapts C++ code originally written by Oli Brady and Nick Golding for the paper [*Global temperature constraints on Aedes aegypti and Ae. albopictus persistence and competence for dengue virus transmission*](https://doi.org/10.1186/1756-3305-7-338) which in turn was adapted from earlier Fortran code written by David L Smith based on the model presented in the paper [*Modelling the global constraints of temperature on transmission of Plasmodium falciparum and P. vivax*](https://doi.org/10.1186/1756-3305-4-92).

### Installation

The modified C++ can be installed as an R package by calling the following from an R session (with the root of this repository as the working directory):
```r
install.packages("tempsuitcalc", type = "source", repos = NULL)
```

### Execution

This repository is an archive of experimental code, and is not intended to be reproducible with execution of a single command. The code to implement the analysis is in the directory `ProjectNew`. The three models (temperature, temperature+humidity, temperature+humidity+precipitation) can be implemented by running different combinations of the code in various files within that directory. Instructions on how to execute these can be found in `owen_readme.txt`. The `.R` files run the code, and the `.slurm` files are used to execute the code on a specific high-performance compute system (University of Melbourne's *Spartan*).

