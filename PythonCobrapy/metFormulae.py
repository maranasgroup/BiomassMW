import pulp
import re

#####How should I define the class for storing the results? Or is it needed?
# class Result(object):
# 	def __init__(self):
# 		#input data
# 		self.model, self.metKnown, self.metUnnown, self.rxnIn = [None] * 4
# 		#pre-processing data
# 		self.ele, self.eleConnect, self.feasTol, self.solver = [None] * 4
# 		#solutions and properties
# 		self.metForm, self.metFormDict, self.rxnBalance, self.infeas, self.bound, self.sol, self.LP = [None] * 7 

def computeMetFormulae(model, metKnown, rxns, metFill=['HCharge1'], findCM=0, nameCM=0,solver):
	useModelCharge = True
	#metabolites with known formulae
	if metKnown is None:
		#default all mets with chemical formulae
		metKnown = {model.metabolites.index(j):j for j in model.metabolites if bool(j.elements)}
	else:
		if (type(metKnown) is list or isinstance(metKnown,type(model.metabolites)))
			and all([isinstance(j, model.metabolites[0]) for j in metKnown]):
			#Dictlist of met objects
			#metK = {model.metabolites.index(j):j for j in metKnown if j in model.metabolites}
			metK = [j for j in metKnown if j in model.metabolites]
		elif all([isinstance(j, str) for j in metKnown]):
			#list of met IDs
			#metK = [model.metabolites.index(j) for j in metKnown if j in model.metabolites]
			#metK = {j: model.metabolites[j] for j in metK}
			metK = [model.metabolites[model.metabolites.index(j)] for j in metKnown if j in model.metabolites]
		else:
			raise ValueError, 'metKnown must be a list of met IDs or met objects.'
		if len(metK) < len(metKnown):
			raise ValueError, 'Some mets in metKnown are not in the model'
		if not all([bool(j.formula) for j in metK])
			print 'Some mets in metKnown do not have formula. Ignore them.'
			#metK = {j: metK[j] for j in metK if bool(metK[j].formula)}
			metK = [j for j in metK if bool(metK[j].formula)]
	#reactions to be balanced for minimum inconsistency
	if rxns is None:
		#all non-exchange reactions
		#rxns = {model.reactions.index(j):j for j in model.reactions if bool(j.reactants) and bool(j.products)}
		rxns = [j for j in model.reactions if bool(j.reactants) and bool(j.products)]
	else:
		if (type(rxns) is list or isinstance(rxns,type(model.reactions)))
			and all([isinstance(j, model.reactions[0]) for j in rxns]):
			#Dictlist of rxn objects
			#rxnK = {model.reactions.index(j):j for j in rxns if j in model.reactions}
			rxnK = [j in rxns if j in model.reactions]
		elif all([isinstance(j, str) for j in rxns]):
			#list of rxn IDs
			#rxnK = [model.reactions.index(j) for j in rxns if j in model.reactions]
			#rxnK = {j: model.reactions[j] for j in rxnK}
			rxnK = [model.reactions[model.reactions.index(j)] for j in rxns if j in model.reactions]
		else:
			raise ValueError, 'rxns must be a list of rxn IDs or rxn objects.'
		if len(rxnK) < len(rxns):
			raise ValueError, 'Some rxns are not in the model'
	#Transform the formulae into dictionaries. Modify from cobra.core.Metabolite.py
	element_re = re.compile("([A-Z][a-z_]*)(\-?[0-9.]+[0-9.]?|(?=[A-Z])?)")
	metKform = {} #known metabolites
	for j in metK:
		parsed = element_re.findall(j.formula)
		s, dictJ = '', {}
		for (element, count) in parsed:
			s = s + element + count
			if count == '':
				count = '1'
			dictJ[element] = float(count)
		if not s == j.formula:
			raise ValueError, "''%s' has incorrect formula '%s'." %(metK[j].id, metK[j].formula)
		if useModelCharge and "Charge" not in dictJ.keys():
			#eval("model.metabolites." + j.id + '.charge')
			chargeJ = j.charge
			if isinstance(chargeJ, (int, long, float, complex)):
				dicJ['Charge'] = chargeJ
		metKform[j] = dictJ
	metFform = {} #metabolites for filling inconsistency
	for j in metFill:
		parsed = element_re.findall(j)
		s, dictJ = '', {}
		for (element, count) in parsed:
			s = s + element + count
			if count == '':
				count = '1'
			dictJ[element] = float(count)
		if not s == j:
			raise ValueError, "''%s' in metFill has incorrect formula." %(j)
		metFform[j] = dictJ

	
	######## get feasibility tolerance from cobrapy?
	feasTol = 1e-8
	
	#digits rounded to for avoiding numerical issues
	digitRounded = 10

	#find elements that are connected because of metFill. They need to be optimized in the same problem.
	eleConnect, ele, metFillConnect = [], [], []
	ele = list(set([k for j in metKform for k in metKform[j]] + [k for j in metFform for k in metFform[j]]))
	eleSet = set(ele)
	while bool(eleSet):
		eleCur = set([eleSet.pop()])
		while True:
			metFillCur = [j for j in metFill if bool(eleCur.intersection(metFform[j].keys()))]
			eleNext = set([metFform[j].keys() for j in metFillCur]).union(eleCur) #in case of metFillCur = []
			if eleCur == eleNext:
				break
			eleCur = eleNext.copy()
		eleConnect.append(list(eleCur))
		metFillConnect.append(metFillCur)
		eleSet.difference(eleCur)
	#S matrix in dict
	S_ij = {j: {i: j.get_coefficients([i])[0] for i in j.reactants + j.products} for j in model.reactions} 
	metU = [j for j in model.metabolites if j not in metK] #unknown mets
	metUform, S_fill, infeas, bound, sol, LP = {}, {}, [], [], {}, [] #data to be stored
	for k in ['minIncon', 'minFill', 'minForm']: #pre-assignment
		metUform[k] = {i: {} for i in metU}
		S_fill[k] = {}
		sol[k] = {'m': {}, 'xp': {}, 'xn': {}, 'Ap': {}, 'An': {}}
	#Optimize for each entry in eleConnect
	for jEC in range(len(eleConnect)):
		infeasJ, boundJ = {}, {}
		LPj = pulp.LpProblem("MIP_ec" + str(jEC), pulp.LpMinimize)
		

		############### How to correctly set solvers in pulp?
		if solver != None:
			LPj.setSolver(solver)
		############### How to set other user-supplied parameters?


		m_ie = pulp.LpVariable.dicts("m", (eleConnect[jEC], metU), 0) #elemental stoichiometry
		if 'Charge' in eleConnect[jEC]: #allow negative charges
			for i in m_ie['Charge']:
				m_ie['Charge'][i].lowBound = None
		xp_je = pulp.LpVariable.dicts("xp", (eleConnect[jEC], rxnK), 0) #inconsistence variable, +ve part
	 	xn_je = pulp.LpVariable.dicts("xn", (eleConnect[jEC], rxnK), 0) #inconsistence variable, -ve part
		if bool(metFillConnect[jEC]): #add adjustment variable if filling mets associated with the current connected elements
			Ap_ij = pulp.LpVariable.dicts("Ap", (metFillConnect[jEC], rxnK), 0)
			An_ij = pulp.LpVariable.dicts("An", (metFillConnect[jEC], rxnK), 0)
		#objective function
		LPj += pulp.lpSum([xp_je[e][j] + xn_je[e][j] for e in eleConnect[jEC] for j in rxnK]), 'Minimum inconsistency'
		for j in rxnK:
			# S * m_unknown + x_pos - x_neg + m * A_pos - m * A_neg = - S * m_known
			for e in eleConnect[jEC]:
				constraint = pulp.lpSum([S_ij[j][i] * m_ie[e][i] for i in S_ij[j] if i in metU]) + xp_je[e][j] - xn_je[e][j]
				if bool(metFillConnect[jEC]):
					constraint += pulp.lpSum([metFform[i][e] * (Ap_ij[i] - An_ij[i]) for i in metFillConnect[jEC] if e in metFform[i]])  
				LPj += constraint == -sum([S_ij[j][i] * metKform[i][e] for i in S_ij[j] if e in metKform[i] and i in metK]), 'balance_%s_%s' %(j.id, e)
		# minimum inconsistency
		LPj.solve()
		infeasJ['minIncon'] = LPj.infeasibilityGap()
		boundJ['minIncon'] = {e: sum([pulp.value(xp_je[e][j]) + pulp.value(xn_je[e][j]) for j in rxnK]) for e in eleConnect[jEC]}
		
		if not infeasJ['minIncon'] <= feasTol: 
			#infeasible (should not happen)
			for i in metU:
				metUform['minIncon'][i].update({e: float('nan') for e in eleConnect[jEC]})
		else: 
			#feasible
			sol['minIncon']['m'].update({e: {i: pulp.value(m_ie[e][i]) for i in metU} for e in eleConnect[jEC]})
			sol['minIncon']['xp'].update({e: {j: pulp.value(xp_je[e][j]) for j in rxnK} for e in eleConnect[jEC]})
			sol['minIncon']['xn'].update({e: {j: pulp.value(xn_je[e][j]) for j in rxnK} for e in eleConnect[jEC]})
			sol['minIncon']['Ap'].update({i: {j: pulp.value(Ap_je[i][j]) for j in rxnK} for i in metFillConnect[jEC]})
			sol['minIncon']['An'].update({i: {j: pulp.value(An_je[i][j]) for j in rxnK} for i in metFillConnect[jEC]})
			sol['minIncon']['obj'] = pulp.value(LPj.objective)
			for i in metU:
				metUform['minIncon'][i].update({e: pulp.value(m_ie[e][i]) for e in eleConnect[jEC] if pulp.value(m_ie[e][i]) != 0})
			if bool(metFillConnect[jEC]):
				S_fill['minIncon'].update({i: {j: pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j]) for j in rxnK 
					if abs(pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j])) > 1e-6} for i in metFillConnect[jEC]
					if any([pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j]) for j in rxnK])})
			#Add constraint to fix the total inconsistency for each element
			for e in eleConnect[jEC]:
				#rounding to avoid infeasibility due to numerical issues
				LPj += pulp.lpSum([xp_je[e][j] + xn_je[e][j] for j in rxnK]) <= round(boundJ['minIncon'][e], digitRounded), 'incon_%s' %(e)
			
			#minimize total adjustment under minimum inconsistency
			LPj.setobjective(pulp.lpSum([Ap_ij[i][j] + An_ij[i][j] for i in metFillConnect[jEC] for j in rxnK]))
			eps0 = 1e-6
			while True:
				LPj.solve()
				if LPj.infeasibilityGap() <= feasTol or eps0 > 1e-4 + 1e-8:
					break
				eps0 = eps0 * 10
				for e in eleConnect[jEC]:
					if LPj.constraints['incon_%s' %(e)].sense < 0: # '<= constraint'
						LPj.constraints['incon_%s' %(e)].changeRHS(round(boundJ['minIncon'][e] * (1 + eps0), digitRounded))	
					else: # '>=' constraint (should not happen)
						LPj.constraints['incon_%s' %(e)].changeRHS(-round(boundJ['minIncon'][e] * (1 + eps0), digitRounded))	
			infeasJ['minFill'] = LPj.infeasibilityGap()
			boundJ['minFill'] = eps0

			#Optimize for a set of minimal formulae given the minimum inconsistency and adjustment
			#Choose a solution used to constrain the minimal formula problem
			if infeasJ['minFill'] <= feasTol:
				solChoice = 'minFill'
				#feasible, get the solution
				sol['minFill']['m'].update({e: {i: pulp.value(m_ie[e][i]) for i in metU} for e in eleConnect[jEC]})
				sol['minFill']['xp'].update({e: {j: pulp.value(xp_je[e][j]) for j in rxnK} for e in eleConnect[jEC]})
				sol['minFill']['xn'].update({e: {j: pulp.value(xn_je[e][j]) for j in rxnK} for e in eleConnect[jEC]})
				sol['minFill']['Ap'].update({i: {j: pulp.value(Ap_je[i][j]) for j in rxnK} for i in metFillConnect[jEC]})
				sol['minFill']['An'].update({i: {j: pulp.value(An_je[i][j]) for j in rxnK} for i in metFillConnect[jEC]})
				sol['minFill']['obj'] = pulp.value(LPj.objective)
				for i in metU:
					metUform['minFill'][i].update({e: pulp.value(m_ie[e][i]) for e in eleConnect[jEC] if pulp.value(m_ie[e][i]) != 0})
				if bool(metFillConnect[jEC]):
					S_fill['minFill'].update({i: {j: pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j]) for j in rxnK 
						if abs(pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j])) > 1e-6} for i in metFillConnect[jEC]
						if any([pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j]) for j in rxnK])})
			else: 
				#infeasible, should not happen
				solChoice = 'minIncon'
				for i in metU:
					metUform['minFill'][i].update({e: float('nan') for e in eleConnect[jEC]})
			
			#remove constraints on total inconsistency
			for e in eleConnect[jEC]:
				tmp = LPj.constraints.pop('incon_%s' %(e))
			#add positive and negative parts for charge variables
			if 'Charge' in eleConnect[jEC]:
				chargeP_i, chargeN_i = pulp.LpVariable('chargeP', metU, 0), pulp.LpVariable('chargeN', metU, 0)
				for i in metU:
					LPj += m_ie['Charge'][i] - chargeP_i[i] + chargeN_i[i] == 0
			#reset objective to sum(m_ie) + sum(chargeP_i + chargeN_i)
			LPj.setObjective(pulp.lpSum([m_ie[e][i] if e != 'Charge' else chargeP_i[i] + chargeN_i[i] for e in eleConnect[jEC] for i in metU]))
			#fix inconsistency and adjustment
			for j in rxnK:
				for e in eleConnect[jEC]:
					xp_je.lowBound, xp_je.upBound = [round(sol[solChoice]['xp'][e][j], digitRounded)] * 2
					xn_je.lowBound, xn_je.upBound = [round(sol[solChoice]['xn'][e][j], digitRounded)] * 2
				if bool(metFillCOnnect[jEC]):
					for i in metFillConnect[jEC]:
						Ap_ij.lowBound, Ap_ij.upBound = [round(sol[solChoice]['Ap'][i][j], digitRounded)] * 2
						An_ij.lowBound, An_ij.upBound = [round(sol[solChoice]['An'][i][j], digitRounded)] * 2
			#solve
			eps0 = 1e-10
			while True:
				LPj.solve()
				if LPj.infeasibilityGap() <= feasTol or eps0 > 1e-5 + 1e-8:
					break
				eps0 = eps0 * 10
				for j in rxnK:
					for e in eleConnect[jEC]:
						xp_je.lowBound, xp_je.upBound = round(sol[solChoice]['xp'][e][j] * (1 - eps0), digitRounded), round(sol[solChoice]['xp'][e][j] * (1 + eps0), digitRounded)
						xn_je.lowBound, xn_je.upBound = round(sol[solChoice]['xn'][e][j] * (1 - eps0), digitRounded), round(sol[solChoice]['xn'][e][j] * (1 + eps0), digitRounded)
					if bool(metFillCOnnect[jEC]):
						for i in metFillConnect[jEC]:
							Ap_ij.lowBound, Ap_ij.upBound = round(sol[solChoice]['Ap'][i][j] * (1 - eps0), digitRounded), round(sol[solChoice]['Ap'][i][j] * (1 + eps0), digitRounded)
							An_ij.lowBound, An_ij.upBound = round(sol[solChoice]['An'][i][j] * (1 - eps0), digitRounded), round(sol[solChoice]['An'][i][j] * (1 + eps0), digitRounded)
			infeasJ['minForm'] = LPj.infeasibilityGap()
			boundJ['minForm'] = eps0
			if infeasJ['minForm'] <= feasTol:
				#feasible, get the solution
				sol['minForm']['m'].update({e: {i: pulp.value(m_ie[e][i]) for i in metU} for e in eleConnect[jEC]})
				sol['minForm']['xp'].update({e: {j: pulp.value(xp_je[e][j]) for j in rxnK} for e in eleConnect[jEC]})
				sol['minForm']['xn'].update({e: {j: pulp.value(xn_je[e][j]) for j in rxnK} for e in eleConnect[jEC]})
				sol['minForm']['Ap'].update({i: {j: pulp.value(Ap_je[i][j]) for j in rxnK} for i in metFillConnect[jEC]})
				sol['minForm']['An'].update({i: {j: pulp.value(An_je[i][j]) for j in rxnK} for i in metFillConnect[jEC]})
				sol['minForm']['obj'] = pulp.value(LPj.objective)
				for i in metU:
					metUform['minForm'][i].update({e: pulp.value(m_ie[e][i]) for e in eleConnect[jEC] if pulp.value(m_ie[e][i]) != 0})
				if bool(metFillConnect[jEC]):
					S_fill['minForm'].update({i: {j: pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j]) for j in rxnK 
						if abs(pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j])) > 1e-6} for i in metFillConnect[jEC]
						if any([pulp.value(Ap_ij[i][j]) - pulp.value(An_ij[i][j]) for j in rxnK])})
			else: 
				#infeasible, should not happen
				for i in metU:
					metUform['minForm'][i].update({e: float('nan') for e in eleConnect[jEC]})
		#store data
		infeas.append(infeasJ)
		bound.append(boundJ)
		LP.append(LPj)

	result = {}
	#pre-processing data
	result['ele'], result['eleConnect'], result['metFillConnect'], result['metKown'], result['metUnkown'], result['rxn'], result['feasTol'] 
		= ele, eleConnect, metFillConnect, metK, metU, rxnK, feasTol
	#solution info
	result['infeas'], result['bound'], result['LP'] = infeas, bound, LP
	#Choose the final solution: 'minForm' is preferred. If not feasible, 'minFill', then 'minIncon'
	if all([k['minForm'] <= feasTol for k in infeas]):
		sol['stat'] = 'minForm'
	elif all([k['minFill'] <= feasTol for k in infeas]):
		sol['stat'] = 'minFill'
	elif all([k['minIncon'] <= feasTol for k in infeas]):
		sol['stat'] = 'minIncon'
	else:
		print 'No feasible solution can be found. Problematic' 
		sol['stat'] = 'infeasible'
		result['model'] = None
		result['sol'] = sol
		return result

	result['sol'] = sol
	#Find conserved moieties by computing extreme rays
	###########


	###########

	#main useful results
	result['metFormDict'] = {k: metKform.update(metUform[k]) for k in ['minIncon', 'minFill', 'minForm']}
	result['metForm'] = {k: {i: formulaDict2Str(result['metFormDict'][k][i],1) for i in result['metFormDict'][k]} for k in ['minIncon', 'minFill', 'minForm']}
	result['metFormWtCharge'] = {k: {i: formulaDict2Str(result['metFormDict'][k][i],0) for i in result['metFormDict'][k]} for k in ['minIncon', 'minFill', 'minForm']}
	result['rxnBalance'] = {k: {j: {e: sum([S_ij[j][i] * result['metFormDict'][k][i][e] for i in S_ij[j] if e in result['metFormDict'][k][i]]) 
		for e in ele} for j in model.reactions} for k in ['minIncon', 'minFill', 'minForm']}
	
	#create a new model if solution exists
	result['model'] = model.copy()
	for i in result['model'].metabolites:
		i.formula = result['metForm'][sol['stat']][i]
	if 'Charge' in ele:
		for i in result['model'].metabolites:
			i.charge = sol[sol['stat']]['m']['Charge'][i]
			
	return result
