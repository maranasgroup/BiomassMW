import re
from cobra.core import Metabolite
from ..functions import formulaDict2Str

element_re = re.compile("([A-Z][a-z_]*)(\-?[0-9.]+[0-9.]?|(?=[A-Z])?)")
isnan = lambda x: x != x

class GenericFormula(object):
	'''GenericFormula(formula=None)
	Formula object allowing for generic elements (an uppercase letter followed by lowercase letters or _, e.g. Generic_element) and 'Charge' in the chemical formula.
	formula='', None or 'nan' input would all result in the formula being 'nan', the indeterminate state
	For metabolites with zero mass (photons), use 'Mass0' or any element with only 0 stoichiometry. Will be automatically changed into 'Mass0'
	The input formula can also be a dictionary with elements as keys and stoichiometries as values or a cobra.core.Metabolite object. 
	In the latter case, the formula and charge in the metabolite would be used.
	Make sure the '.charge' property in the cobra object 'Metabolite' has the correct charge. Or include all charges
	into the chemical formulae, e.g. 'HCharge1' for proton, 'C3H3O3Charge-1' for pyruvate. Set all .charge to None, float('nan') or 0 to ignore charge balance.
	'''
	def __init__(self, formula=None):
		if isinstance(formula, str) or (formula is None):
			#a string of chemical formula 
			self.formula = formula
		elif isinstance(formula, dict):
			#input is a dictionary with elements as keys and stoichiometries as values
			self.elements = formula
		elif isinstance(formula, Metabolite):
			self.formula = formula.formula
			if 'Charge' not in self.elements and self.__formula != 'nan':
				self.charge = formula.charge if isinstance(formula.charge, (int, long, float)) and not isnan(formula.charge) else 0
		else:
			raise ValueError, 'Input must be a string, dictionary or cobra.core.Metabolite object.'

	def __repr__(self):
		return '%s <GenericFormula object at 0x%s>' %(self.formula, id(self))

	@property
	def formula(self):
		'''Formula of the object'''
		return self.__formula

	@formula.setter
	def formula(self, text):
		'''Reset the formula'''
		self.__formula = 'nan' if text in ['', None] else text
		if text != 'nan':
			#reformat the input text to the default format if not undetermined
			self.__formula = formulaDict2Str(self.elements)

	@property
	def formulaWoCharge(self):
		'''Formula without charge'''
		return re.split('Charge',self.__formula)[0]

	@formulaWoCharge.setter
	def formulaWoCharge(self,text):
		'''Reset the formula but keep the charge'''
		charge = self.charge
		self.__formula = text
		self.updateElements({'Charge': charge})

	@property
	def elements(self):
		'''Return the element-stoichiometry dictionary.
		Modify from cobra.core.Metabolite.py'''
		if not isinstance(self.__formula,str):
			raise ValueError, 'Formula is not a string.'
		if self.__formula == 'nan':
			return {'nan': float('nan')}
		parsed = element_re.findall(self.__formula)
		s, formDict = '', {}
		for (element, count) in parsed:
			s = s + element + count
			if count == '':
				count = '1'
			formDict[element] = float(count) if element not in formDict else formDict[element] + float(count)
		if not s == self.__formula:
			raise ValueError, "Incorrect formula: %s" %self.__formula
		return formDict
	
	@elements.setter
	def elements(self, eleDict):
		'''Reset the formula by an element-stoichiometry dictionary.'''
		self.__formula = formulaDict2Str(eleDict)

	@property
	def charge(self):
		'''Return the charge of the formula'''
		if self.__formula == 'nan':
			return float('nan')
		return self.elements['Charge'] if "Charge" in self.elements else 0

	@charge.setter
	def charge(self, n):
		'''Reset the charge'''
		if self.__formula != 'nan':
			if isinstance(n, (int, long, float)) and n == n and n not in [float('inf'), float('-inf')]: 
				formDict = self.elements
				if n == 0:
					tmp = formDict.pop('Charge',None)
				else:
					formDict['Charge'] = n
				self.__formula = formulaDict2Str(formDict)	

	@property
	def mw(self):
		'''Return the molecular weight of the formula. Nan for indeterminate forms'''
		return float('nan') if self.__formula == 'nan' else sum([elements_and_molecular_weights[e] * self.elements[e] for e in self.elements if e in elements_and_molecular_weights])

	@property
	def unknown(self):
		'''Return True for indeterminate forms'''
		return self.__formula == 'nan'

	@property
	def generic(self):
		'''Return True if the formula contains any generic elements.
		Unknown metabolites are not counted as generic
		'''
		return False if self.__formula == 'nan' else any([e not in elements_and_molecular_weights and e != 'Charge' for e in self.elements])
	
	def updateElements(self, updateDict=None, add=True, **kwargs):
		'''Update the formula by adding the new dictionary updateDict into it.
		When add = True, if an element in updateDict is already in the formula, sum up the stoichiometries in the original formula and updateDict
		When add = False, replace the stoichiometry of an existing element
		kwargs allow for key-value argument inputs, e.g. updateElements(C=1,O=2)
		'''
		formDict = self.elements
		if 'nan' in formDict:
			print 'Unknwon formula cannot be updated using .updateElements'
		else:
			if isinstance(updateDict, dict):
				if add:
					#if an element is in both the original formula and the updating dictionary, add up the coefficient
					formDict.update({e: formDict[e] + m_ie if e in formDict else m_ie for e, m_ie in updateDict.iteritems()})
				else:
					#directly update
					formDict.update(updateDict)
			elif updateDict is not None:
				raise ValueError, 'Input must be a dictionary or keyword arguements (e.g. C=6)'
			formDict.update(kwargs)
			self.__formula = formulaDict2Str(formDict)

elements_and_molecular_weights = {
	'H':   1.007940,	'He':  4.002602,	'Li':  6.941000,	'Be':  9.012182,	'B':   10.811000,
	'C':   12.010700,	'N':   14.006700,	'O':   15.999400,	'F':   18.998403,	'Ne':  20.179700,
	'Na':  22.989770,	'Mg':  24.305000,	'Al':  26.981538,	'Si':  28.085500,	'P':   30.973761,
	'S':   32.065000,	'Cl':  35.453000,	'Ar':  39.948000,	'K':   39.098300,	'Ca':  40.078000,
	'Sc':  44.955910,	'Ti':  47.867000,	'V':   50.941500,	'Cr':  51.996100,	'Mn':  54.938049,
	'Fe':  55.845000,	'Co':  58.933200,	'Ni':  58.693400,	'Cu':  63.546000,	'Zn':  65.409000,
	'Ga':  69.723000,	'Ge':  72.640000,	'As':  74.921600,	'Se':  78.960000,	'Br':  79.904000,
	'Kr':  83.798000,	'Rb':  85.467800,	'Sr':  87.620000,	'Y':   88.905850,	'Zr':  91.224000,
	'Nb':  92.906380,	'Mo':  95.940000,	'Tc':  98.000000,	'Ru':  101.070000,	'Rh':  102.905500,
	'Pd':  106.420000,	'Ag':  107.868200,	'Cd':  112.411000,	'In':  114.818000,	'Sn':  118.710000,
	'Sb':  121.760000,	'Te':  127.600000,	'I':   126.904470,	'Xe':  131.293000,	'Cs':  132.905450,
	'Ba':  137.327000,	'La':  138.905500,	'Ce':  140.116000,	'Pr':  140.907650,	'Nd':  144.240000,
	'Pm':  145.000000,	'Sm':  150.360000,	'Eu':  151.964000,	'Gd':  157.250000,	'Tb':  158.925340,
	'Dy':  162.500000,	'Ho':  164.930320,	'Er':  167.259000,	'Tm':  168.934210,	'Yb':  173.040000,
	'Lu':  174.967000,	'Hf':  178.490000,	'Ta':  180.947900,	'W':   183.840000,	'Re':  186.207000,
	'Os':  190.230000,	'Ir':  192.217000,	'Pt':  195.078000,	'Au':  196.966550,	'Hg':  200.590000,
	'Tl':  204.383300,	'Pb':  207.200000,	'Bi':  208.980380,	'Po':  209.000000,	'At':  210.000000,
	'Rn':  222.000000,	'Fr':  223.000000,	'Ra':  226.000000,	'Ac':  227.000000,	'Th':  232.038100,
	'Pa':  231.035880,	'U':   238.028910,	'Np':  237.000000,	'Pu':  244.000000,	'Am':  243.000000,
	'Cm':  247.000000,	'Bk':  247.000000,	'Cf':  251.000000,	'Es':  252.000000,	'Fm':  257.000000,
	'Md':  258.000000,	'No':  259.000000,	'Lr':  262.000000,	'Rf':  261.000000,	'Db':  262.000000,
	'Sg':  266.000000,	'Bh':  264.000000,	'Hs':  277.000000,	'Mt':  268.000000,	'Ds':  281.000000,
	'Rg':  272.000000,	'Cn':  285.000000,	'Uuq': 289.000000,	'Uuh': 292.000000}