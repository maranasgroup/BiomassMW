# BiomassMW
MATLAB and Python functions for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models.  
Please see the following paper for more details:  
Siu H. J. Chan, Jingyi Cai, Lin Wang, Margaret N. Simons-Senftle, Costas D. Maranas; Standardizing biomass reactions and ensuring complete mass balance in genome-scale metabolic models, Bioinformatics, btx453  
**[Link](https://doi.org/10.1093/bioinformatics/btx453)**

**MATLAB** functions (`MatlabCobraToolbox/` Require **[COBRA Toolbox](https://github.com/opencobra/cobratoolbox)**)  
Main functions:  
1. `computeMetFormulae`  
Compute the chemical formulas of the unknown metabolites given the formulas for a set of known metabolites using a set of reactions.
2. `computeMetMWrange`  
Similar to computeFormulasFillMets but return the minimum and maximum possible MW of the target metabolite.

Other functions:  
1. `checkEleBalance`  
For converting chemical formulas into a matrix and checking the elemental balance of reactions 
2. `checkSolFeas`  
For manually checking the feasibility of CPLEX solutions  
3. `convertMatrixFormulas`  
For converting a matrix into chemical formulas  
4. `MW`  
For calculating molecular weights  
5. `setCplexParam`  
For conveniently setting CPLEX parameters in Matlab  


**Python** functions (`PythonCobrapy/` Require **[Cobrapy](https://github.com/opencobra/cobrapy)**)  
Main object: `MetFormulae` in `metFormulae/MetFormulae.py`.  
It is initialized with a cobra model. Use the following two methods to calculate chemical formulae or the range for biomass MW:  
- `.computeMetForm`: Compute the chemical formulas of the unknown metabolites given the formulas for a set of known metabolites using a set of reactions.  
- `.computeMetRange`: Compute the minimum and maximum possible MW of the target metabolite.
