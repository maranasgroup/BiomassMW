# BiomassMW/PythonCobrapy
Python functions for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models.
Require Cobrapy. See the help files for documentation

Main object: MetFormulae in metFormulae/MetFormulae.py. 
It is initialized with a cobra model. Use the following two methods to calculate chemical formulae or the range for biomass MW:
- .computeMetForm: Compute the chemical formulas of the unknown metabolites given the formulas for a set of known metabolites using a set of reactions.
- .computeMetRange: Compute the minimum and maximum possible MW of the target metabolite.