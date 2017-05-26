import re
import pdb
from cobra.core import Model, Reaction, Metabolite
from cobra.util.solver import set_objective
from .objects.GenericFormula import GenericFormula as Formula
from .functions import solution_infeasibility, active_met_rxn, formulaDict2Str, num2alpha, extreme_rays_from_null_basis
from .objects.ResultStructure import DataObject, PreprocessData, MinInconParsiInfo, ConservedMoietyInfo, ResultStructure

try:
	import cdd
	cddImported = True
except ImportError:
	print 'pycddlib is not installed. Unable to use it for extreme ray calculations. Calculate the null space biasis only'
	cddImported = False

element_re = re.compile("([A-Z][a-z_]*)(\-?[0-9.]+[0-9.]?|(?=[A-Z])?)")

class MetFormulae(DataObject):
	#this class property determines whether conserved moiety calculation can be performed
	__cddImported = cddImported
	def __init__(self, cobra_model):
		#must be associated with a normal cobra model
		if isinstance(cobra_model, Model):
			self.model = cobra_model
		else:
			raise ValueError, 'MetFormulae must be initiated with a cobrapy model'
		#default maximum value for variables
		self.infinity = float('inf')
		#default minimum value for charge variables
		self.negative_infinity = - 10000 # -float('inf') seems to cause error for some solvers
		self.__solver = cobra_model.solver

	def computeMetForm(self, metKnown=None, rxns=None, metFill=['HCharge'], findCM='null', deadend=True, nameCM=0, **kwargs):
		if not hasattr(self, 'pre'):
			#Pre-porocessing
			self.preprocessing(metKnown, rxns, metFill, **kwargs)
		#Solve 'Minimum Inconsistency under Parsimony' to find a set of minimal formulae
		mipInfo = self.minInconParsi(**kwargs)
		#Find conserved moieties by calculating extreme rays
		if not self.__cddImported and findCM =='cdd':
			print 'pycddlib is not installed. Use rational null space basis.'
			findCM = 'null'
		cmInfo = None
		if findCM in ['null', 'cdd'] and mipInfo.final != 'infeasible':
			cmInfo = self.conserved_moieties(mipInfo, findCM, deadend, nameCM)
		#get the updated model and info on solutions
		metFormResults = self.getResultsFromMIPandCM(mipInfo, cmInfo)
		self.metFormResults = metFormResults
		return metFormResults

	def computeMetRange(self, metInterest, metKnown=None, rxns=None, **kwargs):
		if not hasattr(self, 'pre'):
			#Pre-porocessing
			self.preprocessing(metKnown, rxns, metFill=None, **kwargs)
		#Solve 'Minimum Inconsistency under Parsimony' to find a set of minimal formulae
		mipInfo = self.minInconParsi_mwRange(metInterest, **kwargs)
		self.metRangeResults = mipInfo
		return mipInfo

	def preprocessing(self, metKnown=None, rxns=None, metFill = ['HCharge'], **kwargs):
		model = self.model
		#metabolites with known formulae
		if metKnown is None:
			#default all mets with non-generic chemical formulae
			metK = {i: Formula(i) for i in model.metabolites}
			metK = {i: i2 for i, i2 in metK.iteritems() if not i2.unknown} #and not i2.generic}
		else:
			if (type(metKnown) is list or isinstance(metKnown,type(model.metabolites)))\
			and all([isinstance(i, type(model.metabolites[0])) for i in metKnown]):
				#Dictlist or list of met objects
				metK = {i: Formula(i) for i in metKnown if i in model.metabolites}
			elif all([isinstance(i, str) for i in metKnown]):
				#list of met IDs
				metK = {i: Formula(model.metabolites[model.metabolites.index(i)]) for i in metKnown if i in model.metabolites}
			else:
				raise	ValueError, 'metKnown must be a list of met IDs or met objects.'
			if len(metK) < len(metKnown):
				raise	ValueError, 'Some mets in metKnown are not in the model'
			if not all([bool(i.formula) for i in metK]):
				print 'Some mets in metKnown do not have formula. Ignore them.'
				metK = {i:i2 for i,i2 in metK.iteritems() if bool(i2.formula)}
		#unknown mets
		metU = [i for i in model.metabolites if i not in metK]
		#reactions to be balanced for minimum inconsistency
		if rxns is None:
			#all non-exchange reactions
			#rxns = {model.reactions.index(j):j for j in model.reactions if bool(j.reactants) and bool(j.products)}
			rxnK = [j for j in model.reactions if bool(j.reactants) and bool(j.products) and (j.upper_bound != 0 or j.lower_bound != 0)]
		else:
			if (type(rxns) is list or isinstance(rxns,type(model.reactions))) \
			and all([isinstance(j, Reaction) for j in rxns]):
				#Dictlist of rxn objects
				rxnK = [j for j in rxns if j in model.reactions]
			elif all([isinstance(j, str) for j in rxns]):
				#list of rxn IDs
				rxnK = [model.reactions[model.reactions.index(j)] for j in rxns if j in model.reactions]
			else:
				raise ValueError, 'rxns must be a list of rxn IDs or rxn objects.'
			if len(rxnK) < len(rxns):
				raise ValueError, 'Some rxns are not in the model'

		#metabolites for filling inconsistency
		metF = {metFill[i]: Formula(metFill[i]) for i in xrange(len(metFill))}

		#find elements that are connected because of metFill. They need to be optimized in the same problem.
		eleConnect, metFillConnect = [], {}
		eleSet = set([e for i in metK.values() + metF.values() for e in i.elements])
		ele = list(eleSet)
		while bool(eleSet):
			eleCur = set([eleSet.pop()])
			while True:
				metFillCur = [i for i in metF if bool(eleCur.intersection(metF[i].elements))]
				eleNext = set([e for i in metFillCur for e in metF[i].elements]).union(eleCur) #in case of metFillCur = []
				if eleCur == eleNext:
					break
				eleCur = eleNext.copy()
			eleConnect.append(tuple(eleCur))
			metFillConnect[tuple(eleCur)] = metFillCur
			eleSet.difference_update(eleCur)

		#feasibility tolerance
		feasTol = kwargs['tolerance_feasibility'] if 'tolerance_feasibility' in kwargs else 1e-8
		#digits rounded to for avoiding numerical issues
		digitRounded = 10

		pre = PreprocessData()
		pre.metKnown, pre.metUnknown, pre.metFill, pre.rxnKnown, pre.ele, pre.eleConnect, pre.metFillConnect \
		,pre.feasTol, pre.digitRounded = metK, metU, metF, rxnK, ele, eleConnect, metFillConnect, feasTol, digitRounded
		self.pre = pre
		return

	def minInconParsi(self, **kwargs):
		inf, neg_inf = self.infinity, self.negative_infinity
		pre = self.pre
		metK, metU, metF, rxnK, ele, eleConnect, metFillConnect, feasTol, digitRounded \
		= pre.metKnown, pre.metUnknown, pre.metFill, pre.rxnKnown, \
		pre.ele, pre.eleConnect, pre.metFillConnect, pre.feasTol, pre.digitRounded
		model = self.model
		kwargs['objective_sense'] = 'minimize' #always minimize
		#data to be stored
		infeas, bound, solution, metModel, solConstrain, solStat, obj = ({} for i in range(7))
		for k in ['minIncon', 'minFill', 'minForm']: #pre-assignment
			solution[k] = {'m': {}, 'xp': {}, 'xn': {}, 'Ap': {}, 'An': {}}
		
		#Optimize for each connected componenet in eleConnect
		for eCC in eleConnect:
			metModelJ = Model(', '.join(eCC))
			metModelJ.solver = self.__solver
			print 'eCC: %s\tsolver: %s' %(eCC, metModelJ.solver.__class__)
			infeasJ, boundJ, objJ = {}, {}, {}
			constraint = {e: {j: Metabolite(j.id + ',' + e) for j in rxnK} for e in eCC}
			m = {e: {i: Reaction('m_' + i.id + ',' + e) for i in metU} for e in eCC}
			xp = {e: {j: Reaction('xp_' + j.id + ',' + e) for j in rxnK} for e in eCC}
			xn = {e: {j: Reaction('xn_' + j.id + ',' + e) for j in rxnK} for e in eCC}
			#add adjustment variable if filling mets associated with the current connected elements
			Ap = {i: {j: Reaction('Ap_' + metF[i].formula + ',' + j.id) for j in rxnK} for i in metFillConnect[eCC]}
			An = {i: {j: Reaction('An_' + metF[i].formula + ',' + j.id) for j in rxnK} for i in metFillConnect[eCC]}
			# m = {e: {i: Reaction('m_' + i.id + ',' + e) for i in metU} for e in eCC}
			# xp = {e: {j: Reaction('xp_' + j.id + ',' + e, objective_coefficient=1) for j in rxnK} for e in eCC}
			# xn = {e: {j: Reaction('xn_' + j.id + ',' + e, objective_coefficient=0) for j in rxnK} for e in eCC}
			# #add adjustment variable if filling mets associated with the current connected elements
			# Ap = {i: {j: Reaction('Ap_' + metF[i].formula + ',' + j.id) for j in rxnK} for i in metFillConnect[eCC]}
			# An = {i: {j: Reaction('An_' + metF[i].formula + ',' + j.id) for j in rxnK} for i in metFillConnect[eCC]}
			
			for e in eCC:
				for j in rxnK:
					#RHS for each constraint: -sum(S_ij * m^known_ie)
					constraint[e][j]._bound = -sum([S_ij * metK[i].elements[e] for i, S_ij in j._metabolites.iteritems() if i in metK and e in metK[i].elements])
					constraint[e][j]._constraint_sense = 'E'
					#x_pos - x_neg for each constraint
					xp[e][j].add_metabolites({constraint[e][j]: 1})
					xn[e][j].add_metabolites({constraint[e][j]: -1})
					xp[e][j].lower_bound, xn[e][j].lower_bound, xp[e][j].upper_bound, xn[e][j].upper_bound = 0, 0, inf, inf
					#m * A_pos - m * A_neg for each constraint
					for i in metFillConnect[eCC]:
						if e in metF[i].elements:
							Ap[i][j].add_metabolites({constraint[e][j]:  metF[i].elements[e]})
							An[i][j].add_metabolites({constraint[e][j]: -metF[i].elements[e]})
				for i in metU:
					# S_ij x m_ie for each i
					m[e][i].add_metabolites({constraint[e][j]: j._metabolites[i] for j in list(i._reaction) if j in rxnK})
					m[e][i].upper_bound = inf
					m[e][i].lower_bound = neg_inf if e == 'Charge' else 0
			for i in metFillConnect[eCC]:
				for j in rxnK:
					Ap[i][j].lower_bound, An[i][j].lower_bound, Ap[i][j].upper_bound, An[i][j].upper_bound = 0, 0, inf, inf
			
			print 'Add reactions into the model'
			#add reactions into the model
			metModelJ.add_reactions([m[e][i] for e in eCC for i in metU])
			metModelJ.add_reactions([xp[e][j] for e in eCC for j in rxnK])
			metModelJ.add_reactions([xn[e][j] for e in eCC for j in rxnK])
			if bool(metFillConnect[eCC]):
				metModelJ.add_reactions([Ap[i][j] for i in metFillConnect[eCC] for j in rxnK])
				metModelJ.add_reactions([An[i][j] for i in metFillConnect[eCC] for j in rxnK])
                        # pdb.set_trace()
                        # the objective can be set as a dictionary of rxns {rxn_id: 1,...}
                        # check this repo https://github.com/opencobra/cobrapy/blob/devel/cobra/util/solver.py
                        # so the objective of the following section would be the summation of xp and summation of xn
                        # Note that the default value of objective_coefficient is 0, so you do not need to specify m and Ap, An
                        objective_dict = dict()
                        for j in rxnK:
                                objective_dict[j] = 1
                        set_objective(metModelJ, objective_dict)
			# for e in eCC:
			# 	for j in rxnK:
                        #                 print j
                        #                 new_xp = metModelJ.reactions.get_by_id('xp_' + j.id + ',' + e)
                        #                 new_xn = metModelJ.reactions.get_by_id('xn_' + j.id + ',' + e)
                        #                 new_xp.objective_coefficient(1)
                        #                 new_xn.objective_coefficient(1)
                        #                 # pdb.set_trace()
			# 		# xp[e][j].objective_coefficient = 1
                        #                 # xn[e][j].objective_coefficient = 1
                        #         pdb.set_trace()        
                        ## the default is 0, not nded to specify
			# 	for i in metU:
			# 		m[e][i].objective_coefficient = 0
			# # for i in metFillConnect[eCC]:
			# # 	for j in rxnK:
			# # 		Ap[i][j].objective_coefficient, An[i][j].objective_coefficient = 0, 0

			#Solve for minimum inconsistency
			print 'solve'
			sol = metModelJ.optimize(**kwargs)
			solStatJ = 'minIncon'
			infeasJ[solStatJ] = solution_infeasibility(metModelJ, sol)
			if sol.x_dict is None:
				boundJ[solStatJ] = {e: float('nan') for e in eCC}
			else:	
				boundJ[solStatJ] = {e: sum([sol.x_dict[j.id] for j in xp[e].values() + xn[e].values()]) for e in eCC}
			if not infeasJ[solStatJ] <= feasTol:
				#infeasible (should not happen)
				infeasJ['minFill'], infeasJ['minForm'] = inf, inf
				solConstrainJ = 'infeasible'
				solStatJ = 'infeasible'
				objJ['minIncon'], objJ['minFill'], objJ['minForm'] = (float('nan') for i in range(3))
			else:
				#Feasible. Store the solution
				solution[solStatJ]['m'].update({e: {i: sol.x_dict[m[e][i].id] for i in metU} for e in eCC})
				solution[solStatJ]['xp'].update({e: {j: sol.x_dict[xp[e][j].id] for j in rxnK} for e in eCC})
				solution[solStatJ]['xn'].update({e: {j: sol.x_dict[xn[e][j].id] for j in rxnK} for e in eCC})
				solution[solStatJ]['Ap'].update({i: {j: sol.x_dict[Ap[i][j].id] for j in rxnK} for i in metFillConnect[eCC]})
				solution[solStatJ]['An'].update({i: {j: sol.x_dict[An[i][j].id] for j in rxnK} for i in metFillConnect[eCC]})
				objJ[solStatJ] = sol.f

				solConstrainJ = 'minIncon'
				#minimize total adjustment if filling metabolites exist
				if bool(metFillConnect[eCC]):
					#Add constraint to fix the total inconsistency for each element
					constraint_minIncon = {e: Metabolite('minIncon_'+e) for e in eCC}
					for e in eCC:
						#rounding to avoid infeasibility due to numerical issues
						constraint_minIncon[e]._bound = round(boundJ[solStatJ][e], digitRounded)
						constraint_minIncon[e]._constraint_sense = 'L'
						for j in rxnK:
							xp[e][j].add_metabolites({constraint_minIncon[e]: 1})
							xn[e][j].add_metabolites({constraint_minIncon[e]: 1})
							xp[e][j].objective_coefficient = 0
							xn[e][j].objective_coefficient = 0
					for i in metFillConnect[eCC]:
						for j in rxnK:
							Ap[i][j].objective_coefficient = 1
							An[i][j].objective_coefficient = 1
					solStatJ = 'minFill'
					eps0 = 1e-6
					while True:
						sol = metModelJ.optimize(**kwargs)
						infeasJ[solStatJ] = solution_infeasibility(metModelJ, sol)
						if infeasJ[solStatJ] <= feasTol or eps0 > 1e-4 + 1e-8:
							break
						eps0 *= 10
						for e in eCC:
							#rounding to avoid infeasibility due to numerical issues
							constraint_minIncon[e]._bound = round(boundJ['minIncon'][e] * (1 + eps0), digitRounded)
					boundJ[solStatJ] = eps0
					if infeasJ[solStatJ] <= feasTol:
						#Feasible. Use this as the solution for constraining the minimal formula problem
						solConstrainJ = 'minFill'
						#Store the solution
						solution[solStatJ]['m'].update({e: {i: sol.x_dict[m[e][i].id] for i in metU} for e in eCC})
						solution[solStatJ]['xp'].update({e: {j: sol.x_dict[xp[e][j].id] for j in rxnK} for e in eCC})
						solution[solStatJ]['xn'].update({e: {j: sol.x_dict[xn[e][j].id] for j in rxnK} for e in eCC})
						solution[solStatJ]['Ap'].update({i: {j: sol.x_dict[Ap[i][j].id] for j in rxnK} for i in metFillConnect[eCC]})
						solution[solStatJ]['An'].update({i: {j: sol.x_dict[An[i][j].id] for j in rxnK} for i in metFillConnect[eCC]})
						objJ[solStatJ] = sol.f
					else:
						#infeasible, should not happen
						objJ[solStatJ] = float('nan')
					#prepare to compute minimal formulae
					eps0 = 1e-10
					for j in rxnK:
						for i in metFillConnect[eCC]:
							#reset the objective function
							Ap[i][j].objective_coefficient = 0
							An[i][j].objective_coefficient = 0
							#Fix the adjustment
							Ap[i][j].lower_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 - eps0), digitRounded)
							Ap[i][j].upper_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 + eps0), digitRounded)
							An[i][j].lower_bound = round(solution[solConstrainJ]['An'][i][j] * (1 - eps0), digitRounded)
							An[i][j].upper_bound = round(solution[solConstrainJ]['An'][i][j] * (1 + eps0), digitRounded)
						for e in eCC:
							#remove all coefficients on the constraints for total inconsistency
							tmp = xp[e][j]._metabolites.pop(constraint_minIncon[e])
							tmp = xn[e][j]._metabolites.pop(constraint_minIncon[e])
					for e in eCC:
						#then remove the constraints
						metModelJ.metabolites.remove(constraint_minIncon[e])
				else:
					#no filling metabolites
					infeasJ['minFill'], boundJ['minFill'] = 0, 0

				#compute minimal formulae
				eps0 = 1e-10
				for e in eCC:
					for j in rxnK:							
						#fix the inconsistency
						xp[e][j].lower_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 - eps0), digitRounded)
						xp[e][j].upper_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 + eps0), digitRounded)
						xn[e][j].lower_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 - eps0), digitRounded)
						xn[e][j].upper_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 + eps0), digitRounded)
					
					if e != 'Charge':
						for i in metU:
							#objective coefficients for minimal formulae
							m[e][i].objective_coefficient = 1
					else:	
						#add variables for the positive and negative part of charges
						chargePos = {i: Reaction('chargePos_' + i.id) for i in metU}
						chargeNeg = {i: Reaction('chargeNeg_' + i.id) for i in metU}
						constraint_chargeDecomp = {i: Metabolite('chargeDecomp_' + i.id) for i in metU}
						for i in metU:
							#constraint_chargeDecomp[i]: m_charge,i - chargePos_i + chargeNeg_i = 0
							m[e][i].add_metabolites({constraint_chargeDecomp[i]: 1})
							chargePos[i].add_metabolites({constraint_chargeDecomp[i]: -1})
							chargeNeg[i].add_metabolites({constraint_chargeDecomp[i]: 1})
							chargePos[i].objective_coefficient, chargePos[i].lower_bound, chargePos[i].upper_bound = 1, 0, inf
							chargeNeg[i].objective_coefficient, chargeNeg[i].lower_bound, chargeNeg[i].upper_bound = 1, 0, inf
							constraint_chargeDecomp[i]._bound, constraint_chargeDecomp[i]._constraint_sense = 0, 'E'
						metModelJ.add_metabolites(constraint_chargeDecomp.values())
						metModelJ.add_reactions(chargePos.values() + chargeNeg.values())
				#Solve for minimum formulae
				solStatPrev, solStatJ =solStatJ, 'minForm'
				while True:
					sol = metModelJ.optimize(**kwargs)
					infeasJ[solStatJ] = solution_infeasibility(metModelJ, sol)
					if infeasJ[solStatJ] <= feasTol or eps0 > 1e-5 + 1e-8:
						break
					eps0 *= 10
					#relax bounds, rounding to avoid infeasibility due to numerical issues
					for j in rxnK:
						for i in metFillConnect[eCC]:
							Ap[i][j].lower_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 - eps0), digitRounded)
							Ap[i][j].upper_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 + eps0), digitRounded)
							An[i][j].lower_bound = round(solution[solConstrainJ]['An'][i][j] * (1 - eps0), digitRounded)
							An[i][j].upper_bound = round(solution[solConstrainJ]['An'][i][j] * (1 + eps0), digitRounded)
						for e in eCC:
							xp[e][j].lower_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 - eps0), digitRounded)
							xp[e][j].upper_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 + eps0), digitRounded)
							xn[e][j].lower_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 - eps0), digitRounded)
							xn[e][j].upper_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 + eps0), digitRounded)
				boundJ[solStatJ] = eps0
				if infeasJ[solStatJ] <= feasTol:
					#Feasible. Store the solution
					solution[solStatJ]['m'].update({e: {i: sol.x_dict[m[e][i].id] for i in metU} for e in eCC})
					solution[solStatJ]['xp'].update({e: {j: sol.x_dict[xp[e][j].id] for j in rxnK} for e in eCC})
					solution[solStatJ]['xn'].update({e: {j: sol.x_dict[xn[e][j].id] for j in rxnK} for e in eCC})
					solution[solStatJ]['Ap'].update({i: {j: sol.x_dict[Ap[i][j].id] for j in rxnK} for i in metFillConnect[eCC]})
					solution[solStatJ]['An'].update({i: {j: sol.x_dict[An[i][j].id] for j in rxnK} for i in metFillConnect[eCC]})
					objJ[solStatJ] = sol.f
				else:
					#infeasible, should not happen
					solStatJ = solStatPrev
					objJ[solStatJ] = float('nan')

			#store data
			infeas[eCC] = infeasJ
			bound[eCC] = boundJ
			obj[eCC] = objJ
			metModel[eCC] = metModelJ
			solConstrain[eCC] = solConstrainJ
			solStat[eCC] = solStatJ

		#summarize the final solution state
		if any([k == 'infeasible' for k in solStat.values()]):
			print 'Failure: no feasible solution can be found.' 
			solFinal = 'infeasible'
		else:
			solMixed = True
			for stat in ['minForm', 'minFill', 'minIncon']:
				if all([k == stat for k in solStat.values()]):
					solFinal, solMixed = stat, False
					break
			if solMixed: 
				solFinal = 'mixed'

		#Get the resultant set of formulae. For each set of elements in eleConnect, choose the latest solution (minForm > minFill > minIncon), recorded in solStat.
		if solFinal != 'infeasible':
			formulae = {i: Formula('Mass0') for i in metU}
			for i in metU:
				for eCC in eleConnect:
					formulae[i].updateElements({e: round(solution[solStat[eCC]]['m'][e][i], digitRounded) for e in eCC})
		else:
			formulae = {i: Formula() for i in metU}

		formulae.update(metK)
		mipInfo = MinInconParsiInfo()
		mipInfo.formulae, mipInfo.infeas, mipInfo.bound, mipInfo.obj, mipInfo.solution, \
		mipInfo.metModel, mipInfo.final, mipInfo.solConstrain, mipInfo.solStat \
		= formulae, infeas, bound, obj, solution, metModel, solFinal, solConstrain, solStat
		
		return mipInfo

	def conserved_moieties(self, mipInfo, findCM='null', deadend=True, nameCM=0):
		#Find conserved moieties by computing extreme rays using the cdd library
		#Need to have pycddlib installed
		model = self.model
		if (not deadend) or nameCM == 1:
			activeMets, activeRxns = active_met_rxn(model)
		else:
			activeMets, activeRxns = model.metabolites, model.reactions
		if not deadend:
			cmMets, cmRxns = activeMets, activeRxns
		else:
			cmMets, cmRxns = model.metabolites, model.reactions
		
		#(matrix format: [[b_1, a_11, a_12, ..., a_1N], ..., [b_M, a_M1, a_M2, ..., a_MN]] for Ax + b >= 0 
		#where A = [a_ij], b = [b_1, ..., b_M])
		
		if findCM == 'null':
			#transpose of S
			S = [[j._metabolites[i] if i in j._metabolites else 0 for i in cmMets] for j in cmRxns]
			#This method calculates a rational basis for transpose(S) from the reduced row echelon form, usually find a subset of extreme rays, quite probably the whole set. 
			N = extreme_rays_from_null_basis(S)
			cm = [{cmMets[i]: N[i,k] for i in xrange(len(cmMets)) if N[i,k] != 0} for k in xrange(N.shape[1])]
		elif findCM == 'cdd':
			#transpose(S) >= 0 
			S = [[0] + [j._metabolites[i] if i in j._metabolites else 0 for i in cmMets] for j in cmRxns]
			# #transpose(S) <= 0 
			S += [[-i for i in j] for j in S]
			#all entries >= 0
			S += [[0] + [1 if i == j else 0 for i in xrange(len(cmMets))] for j in xrange(len(cmMets))]
			print 'Matrix size for cdd extreme ray calculation: %d x %d' %(len(S), len(S[0]))
			# The cdd library seems unable to cope with genome-scale models. The best is to call EFMtool. To be implemented.
			mat = cdd.Matrix(S, number_type='float')
			mat.rep_type = cdd.RepType.INEQUALITY
			poly = cdd.Polyhedron(mat)
			ext = poly.get_generators()
			cm = [{cmMets[i]: ext.__getitem__(k)[i+1] for i in xrange(len(cmMets)) if ext.__getitem__(k)[i+1] != 0} for k in xrange(ext.row_size)]

		#generic conserved moieties involing no known metabolites
		cmGeneric = [c for c in cm if all([i in self.pre.metUnknown for i in c])]
		cmGenericDict = {}
		NcmDefault = 0		
		for c in cmGeneric:
			#Use defaulted names for dead end metabolites if nameCM = 1, or always use defaulted names if nameCM = 0
			if nameCM == 0 or (any([not i in activeMets for i in c.keys()]) and nameCM == 1):
				#defaulted names
				NcmDefault += 1
				cmNameCur = 'Conserve_' + num2alpha(NcmDefault)
			else:
				print '\n\n'
				for i in c.keys():
					toPrint = self.pre.metKnown[i].formula if i in self.pre.metKnown else mipInfo.formulae[i].formula
					if toPrint == 'Mass0':
						toPrint = ''
					toPrint += formulaDict2Str({"Conserve": c[i]})
					print '%s\t%s\t%s' %(i.id, i.name, toPrint)			
				while True:
					cmNameCur = raw_input("\nEnter the formula for the conserved moiety: " \
					+ "(e.g. C2HRab_cd0.5Charge-1 -> {C:2, H:1, Rab_cd: 0.5, Charge: -1}, " \
					+ "hit return to use default name 'Conserve_xxx')\n")
					#check if the input is empty or a correct formula
					if cmNameCur == "" or ''.join([''.join(k) for k in element_re.findall(cmNameCur)]) == cmNameCur:
						break
					print 'Incorrect format of the input formula!\n'
				if cmNameCur == '':
					#empty string means using the default name
					NcmDefault += 1
					cmNameCur = 'Conserve_' + num2alpha(NcmDefault)
				
			cmGenericDict[Formula(cmNameCur)] = c

		cmInfo = ConservedMoietyInfo()
		cmInfo.cm, cmInfo.cmGeneric, cmInfo.cmGenericDict = cm, cmGeneric, cmGenericDict
  		return cmInfo
  
	def getResultsFromMIPandCM(self, mipInfo, cmInfo=None):
		
		metFormResults = ResultStructure()

		if mipInfo.final =='infeasible':
			metFormResults.mipInfo = mipInfo
			metFormResults.model, metFormResults.formulae, metFormResults.formulaeWoCharge,  \
			metFormResults.formulaeDict, metFormResults.rxnBalance, metFormResults.cmInfo = (None for k in range(6))
			self.metFormResults = metFormResults
			return metFormResults

		#Get S_fill. For each set of elements in eleConnect, choose the latest solution (minForm > minFill > minIncon), recorded in mipInfo.solStat.
		S_fill = {}
		for eCC in self.pre.eleConnect:
			for i in self.pre.metFillConnect[eCC]:
				d = {j: round(mipInfo.solution[mipInfo.solStat[eCC]]['Ap'][i][j] - mipInfo.solution[mipInfo.solStat[eCC]]['An'][i][j], self.pre.digitRounded) for j in self.pre.rxnKnown
					if abs(mipInfo.solution[mipInfo.solStat[eCC]]['Ap'][i][j] - mipInfo.solution[mipInfo.solStat[eCC]]['An'][i][j]) > 1e-6}
				if bool(d):
					S_fill[i] = d
		
		formulae = mipInfo.formulae.copy()	
		if cmInfo != None and bool(cmInfo.cmGenericDict):
			#include the unknown conserved moieties into the chemical formulae
  			for cmNameCur, c in cmInfo.cmGenericDict.iteritems():
				for i in c:
					formulae[i].updateElements({e: c[i]*m_ie for e, m_ie in cmNameCur.elements.iteritems()}, True)
					
		#elemental balance of each reaction
		rxnBalance = {j: {e: sum([j._metabolites[i] * formulae[i].elements[e] for i in j._metabolites if e in formulae[i].elements])
			for e in self.pre.ele} for j in self.model.reactions}
		rxnBalance = {j: {e: rxnBalance[j][e] for e in rxnBalance[j] if rxnBalance[j][e] != 0} for j in rxnBalance}
		#updated model
		model = self.model.copy()
		for i in model.metabolites:
			#the original metabolite object
			i0 = self.model.metabolites[self.model.metabolites.index(i.id)]
			i.formula, i.charge = formulae[i0].formulaWoCharge, formulae[i0].charge
		
		metFormResults.model, metFormResults.formulae, metFormResults.S_fill, metFormResults.rxnBalance, \
		metFormResults.mipInfo, metFormResults.cmInfo = model, formulae, S_fill, rxnBalance, mipInfo, cmInfo
		self.metFormResults = metFormResults
		return metFormResults

	def minInconParsi_mwRagne(self, **kwargs):
		inf, neg_inf = self.infinity, self.negative_infinity
		pre = self.pre
		metK, metU, metF, rxnK, ele, feasTol, digitRounded \
		= pre.metKnown, pre.metUnknown, pre.rxnKnown, \
		pre.ele, pre.eleConnect, pre.metFillConnect, pre.feasTol, pre.digitRounded
		model = self.model
		kwargs['objective_sense'] = 'minimize'
		infeas, bound, solution, metModel, solStat, obj = ({} for i in range(6))
		for k in ['minIncon', 'minFill', 'minForm']: #pre-assignment
			solution[k] = {'m': {}, 'xp': {}, 'xn': {}}

		for e in ele:
			metModelJ = Model('min/max ' + e)
			infeasJ, boundJ, objJ = {}, {}, {}
			constraint = {j: Metabolite(j.id + ',' + e) for j in rxnK}
			m = {i: Reaction('m_' + i.id + ',' + e) for i in metU}
			xp = {j: Reaction('xp_' + j.id + ',' + e) for j in rxnK}
			xn = {j: Reaction('xn_' + j.id + ',' + e) for j in rxnK}
			for j in rxnK:
				#RHS for each constraint: -sum(S_ij * m^known_ie)
				constraint[j]._bound = -sum([S_ij * metK[i].elements[e] for i, S_ij in j._metabolites.iteritems() if i in metK and e in metK[i].elements])
				constraint[j]._constraint_sense = 'E'
				#x_pos - x_neg for each constraint
				xp[j].add_metabolites({constraint[j]: 1})
				xn[j].add_metabolites({constraint[j]: -1})
				xp[j].lower_bound, xn[j].lower_bound, xp[j].upper_bound, xn[j].upper_bound = 0, 0, inf, inf
				xp[j].objective_coefficient, xn[j].objective_coefficient = 1, 1


	@property
	def solver(self):
		return self.__solver

	@solver.setter
	def solver(self, choice):
		self.model.solver = choice
		self.__solver = self.model.solver
	
