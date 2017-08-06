class DataObject(object):
	'''DataObject has a representation summarizing all public attributes
	'''
	def __repr__(self):
		print 'Attribute\tClass\tValue/Length'
		for k in dir(self):
			if k[0] != '_':
				val = eval("self." + k)
				s = ''
				if isinstance(val, (int, long, float, complex)):
					s = str(val)
				elif isinstance(val, (list, dict, tuple, set)):
					s = str(len(val))
				elif isinstance(val, str):
					s = val
				print '%s\t%s\t%s\t at 0x%x' %(k, type(val), s, id(val))
		print ''
		return "<%s at 0x%x>" % (self.__class__.__name__, id(self))

class PreprocessData(DataObject):
	'''Preprocess data:
	metKnown: {i: GenericFormula(i)} for all metabolite i in the input model
	metUnknown: list of all unknown metabolite objects in the input model
	metFill: {'X1Y2Z3': GenericFormula('X1Y2Z3')} for all filling metabolite (e.g. with formula 'X1Y2Z3')
	rxnKnown: list of reaction objects in the input model used for balancing
	ele: list of elements existing in the chemical formulae of metKnown
	eleConnect: list of connected elements in tuples
	metFillConnect: {eleConnent: metFill}, dictionary associating each filling metabolite metFill to the connected sets of elements eleConnect.
	'''
	def __init__(self):
		self.metKnown, self.metUnknown, self.metFill, self.rxnKnown, \
		self.ele, self.eleConnect, self.metFillConnect = (None for i in range(7))

class MinInconParsiInfo(DataObject):
	'''MinInconParsiInfo object summarizing the results of solving MIP.
		formulae: result formulae
		mwRange: (min MW, max MW) (computeMetRange only)
		rhs: {e: {i: RHS[i] for all metabolite i}} the RHS value in the MIP problem for each element solved. (computeMetRange only)
		infeas: infeasibility of each solve
		bound: bound used for total inconsistency or the relaxation value eps0 for each solve
		obj: objective function value for each solve
		solution: solution values for each type of variables (m, xp, xn, Ap, An)
		metModel: the MIP problems solved for each element/connected set of elements as individual cobra models
		final: the final status of the solution
		solStat: solution status for each element/connected set of elements
		solConstrain: 'minIncon' or 'minFill' indicating the solution used to constrain the minimal formulae problem. (computeMetForm only)
	'''
	def __init__(self):
		self.formulae, self.mwRange, self.rhs, self.infeas, self.bound, self.obj, self.solution, self.metModel, \
		self.final, self.solStat, self.solConstrain = (None for i in range(11))
		
class ConservedMoietyInfo(DataObject):
	'''ConservedMoietyInfo object containing the conserved moiety information:
		cm: list of conserved moieties, each being a metabolite-coefficient dictionary
		cmGeneric: list of entries in cm that are generic (no known metabolites involved)
		cmGenericDict: {fomrula: cmGeneric} dictionary. formula is the GenericFormula object with the defaulted or inputted name as the formula for the conserved moiety.
	'''
	def __init__(self):
		self.cm, self.cmGeneric, self.cmGenericDict = (None for i in range(3))
  
class ResultStructure(DataObject):
	'''ResultStructure object with the following attributes:
		formulae: the final formulae for unknown metabolites
		S_fill: a dictionary of stoichiometric coefficients for filling metabolites
		rxnBalance: {rxn:{e: balance}} dictionary for the elemental balance for each element e of each reaction rxn  
		model: a copy of the input cobra model with updated chemical formulae
		mipInfo: MinInconParsiInfo object summarizing the results of solving MIP.
		cmInfo: ConservedMoietyInfo object containing the conserved moiety information.
	'''
	def __init__(self):
		self.formulae, self.S_fill, self.rxnBalance, self.model, self.mipInfo, self.cmInfo = (None for i in range(6))
