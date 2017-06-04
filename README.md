# BiomassMW
MATLAB and Python functions for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models.
Require Cobrapy.

MATLAB functions (MatlabCobraToolbox):
Main functions: 
1. computeMetFormulae
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


Python functions (PythonCobrapy)
Main object: MetFormulae in metFormulae/MetFormulae.py. 
It is initialized with a cobra model. Use the following two methods to calculate chemical formulae or the range for biomass MW:
- .computeMetForm: Compute the chemical formulas of the unknown metabolites given the formulas for a set of known metabolites using a set of reactions.
- .computeMetRange: Compute the minimum and maximum possible MW of the target metabolite.