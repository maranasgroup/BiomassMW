# BiomassMW
Functions for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models.
Require COBRA toolbox and CPLEX

Main functions: 
1. computeFormulasFillMets
Compute the chemical formulas of the unknown metabolites given the formulas for a set of known metabolites using a set of reactions.
2. computeMetMWrange
Similar to computeFormulasFillMets but return the minimum and maximum possible MW of the target metabolite.

Other functions:
1. checkEleBalance
For converting chemical formulas into a matrix and checking the elemental balance of reactions 
2. checkSolFeas
For manually checking the feasibility of CPLEX solutions
3. convertMatrixFormulas
For converting a matrix into chemical formulas
4. MW
For calculating molecular weights
5. setCplexParam
For conveniently setting CPLEX parameters in Matlab
