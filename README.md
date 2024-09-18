This is a repository for the code for the '**Microbial Coexistence with Resource Toxicity**' Research Project.

**Files:**\
The 2 key files are:
1) _Nondimensionalised_fast_\
   This code creates 2-parameter bifurcation diagrams for the non-dimensionalised version of the 2-species model.
2) _n_species_\
   This code creates 2-parameter bifurcation diagrams for the non-dimensionalised version of the n-species model, however it creates one for each type of equilibrium. Note that _Nondimensionalised_fast_ creates 1 full diagram that describes the bifurcations for all of the equilibria types.

**Instructions on how to use these 2 files are found below, and in the code itself.**

The other files are:
1) (OLD)_2parameter_\
  This code was an old draft for the 2-species model, no longer used.
2) _Root Finder_\
  This code is exclusively for the Root Finder used in the 2 key files above.
3) _jacobian_n_species_\
  This code finds the Jacobian matrix for n-species, as used in _n_species_.
4) _nondimensionalised_\
  This code was an old, slower draft of _Nondimensionalised_fast_.


**Instructions for _Nondimensionalised_fast_**
1) Alter your equations for r_1 and r_2
2) Alter your equations for r_1' and r_2' (you will have to calculate these by hand)
3) Choose your parameters, and hide the two you will be plotting on the graph
4) Choose the length and width of the plots (do not set xmin or ymin to be 0, use a small value like 0.01 instead else it will break)
5) Choose the number of points you want to plot (e.g. 100x100)
6) Place the two parameters you want to plot into the for loops

**Instructions for _n_species_**
1) Change the number of species
2) Input your equatios for r(M) as a vector
3) Input your equations for r'(M) as a vector (again, you must calculate these by hand)
4) Change your values of I and P, or hide if you want to plot them on the graph
5) Choose the length and width of the plots (again, do not set xmin or ymin to be 0, use a small value like 0.01 instead else it will break)
6) Choose the number of points you want to plot (e.g. 100x100)
7) Place the parameter you want on the y-axis into the for loop
8) Place the parameter you want on the x-axis into the for loop
9) Change your values of A (or hide if you want to plot on the graph)
10) Change the structure of A based on the number of species
