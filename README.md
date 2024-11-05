This is a repository for the code for the '**Microbial Coexistence with Resource Toxicity**' Research Project.

**Files:**\
The 6 key files are:
1) _bifurcation_plots_\
   This is a live script where users can input their own growth functions and parameter values to get a bifurcation plot of metal influx i against competition strength gamma \
2) _one_species_\
   This code creates a 2-parameter bifurcation diagram for p and i, showing the number and types of each equilibria, for the 1-species model.
3) _two_species_ \
   This code creates 2-parameter bifurcation diagrams for p and i of two sorts for the 2-species model: one sort that shows which types of equilibria are stable, and one sort of separate subplots for each type of equilibria that show how many there are.
4) _three_species_simple_ \
   This code creates a 2-parameter bifurcation diagram for p and i that shows the number of non-extinction equilibria for the 3-species model.   
5) _Nondimensionalised_fast_\
   This code creates 2-parameter bifurcation diagrams for any two parameters of your choosing, showing the types of each equilibria, for the non-dimensionalised version of the 2-species model.
6) _n_species_\
   This code creates 2-parameter bifurcation diagrams for any two parameters of your choosing for the non-dimensionalised version of the n-species model, however it creates one subplot for each type of equilibrium. Note that _Nondimensionalised_fast_ creates 1 full diagram that describes the bifurcations for all of the equilibria types.

**Instructions on how to use these 6 files are found in the code itself.**

The other files are:
1) (OLD)_2parameter_\
  This code was an old draft for the 2-species model, no longer used.
2) _Root Finder_\
  This code is exclusively for the Root Finder used in the 2 key files above.
3) _jacobian_n_species_\
  This code finds the Jacobian matrix for n-species, as used in _n_species_.
4) _nondimensionalised_\
  This code was an old, slower draft of _Nondimensionalised_fast_.
