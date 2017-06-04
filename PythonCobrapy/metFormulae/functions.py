import numpy

def solution_infeasibility(model, sol):
	'''Return the infeasibility of the solution returned by model.optimize() for a cobra model 'model'
	'''
	if sol.fluxes is None:
		return float('inf')
	if False:
		#obsolete
		Svb = {i: sum([j._metabolites[i] * sol.fluxes[j.id] for j in list(i._reaction)]) - i._bound for i in model.metabolites}
		infeas = max([abs(Svb[i]) for i in model.metabolites if i._constraint_sense == 'E'] + [0])
		infeas = max([Svb[i] for i in model.metabolites if i._constraint_sense == 'L'] + [infeas])
		infeas = max([- Svb[i] for i in model.metabolites if i._constraint_sense == 'G'] + [infeas])
	else:
		infeas = max([sum([j._metabolites[i] * sol.fluxes[j.id] for j in list(i._reaction)]) - model.constraints[i.id].ub for i in model.metabolites])
		infeas = max([infeas] + [model.constraints[i.id].lb - sum([j._metabolites[i] * sol.fluxes[j.id] for j in list(i._reaction)]) for i in model.metabolites])
	return infeas

def active_met_rxn(model):
	'''Return [list of active metabolites, list of active reactions] of a COBRA model. Being active means being not dead end.
	'''
	deadMet, deadRxn = [], []
	while True:
		finish = True
		for i in model.metabolites:
			if i not in deadMet:
				rxn = i._reaction.difference(deadRxn)
				#only 1 or no active reaction or all active reactions have the same possible flux sign (all >= 0, or all <= 0)
				if len(rxn) <= 1 \
				or all([j._metabolites[i] * j.upper_bound >= 0 for j in list(rxn)] + [j._metabolites[i] * j.lower_bound >= 0 for j in list(rxn)]) \
				or all([j._metabolites[i] * j.upper_bound <= 0 for j in list(rxn)] + [j._metabolites[i] * j.lower_bound <= 0 for j in list(rxn)]):
					deadMet.append(i)
					deadRxn += [j for j in list(rxn)]
					finish = False
					break
		if finish:
			break
	#return the metabolites and reactions that are not dead end
	return [[i for i in model.metabolites if i not in deadMet], [j for j in model.reactions if j not in deadRxn]]

def formulaDict2Str(formDict, wtCharge=True, dMax=12):
	'''Convert an element-stoichiometry dictionary into a string of chemical formula.
	'nan' for indeterminate forms. 'Mass0' for empty dictionary or all stoichiometries = 0.
	The elements C, H, N, O, P, S are prioritized if present. 'Charge' is present at the end stating the charge of the metabolite
	'''
	if 'nan' in formDict or any([x != x or x in [float('infinity'), -float('infinity')] for x in formDict.values()]):
		return 'nan'
	if not bool(formDict) or not any(formDict.values()):
		#Empty dictionary or all m_ie being zero
		return 'Mass0'
	#prioritize elements
	elePrior = ['C', 'H', 'N', 'O', 'P', 'S']
	ele = [e for e in elePrior if e in formDict] + [e for e in formDict if e not in elePrior + ['Charge']]
	ele += ['Charge'] if wtCharge and 'Charge' in formDict else []
	s, sFormat = '', '%.' + str(dMax) + 'f'
	for e in ele:
		if formDict[e] != 0:
			if formDict[e] == 1 and e != 'Charge':
				#Add the element symbol only if m_ie = 1
				s += e
			elif abs(round(formDict[e]) - formDict[e]) < 1e-10:
				#integer (avoid all the checking below)
				s += e + '%d' %formDict[e]
			elif abs(formDict[e]) > 10^(-dMax):
				#round to the largest digit that gives a difference < 1e-10
				n, d = formDict[e], 0
				while abs(round(n, d) - n) > 1e-10 and d < dMax:
					d += 1
				n = round(n, d)
				#trim all the tailing 0s
				sE = sFormat %n
				if '.' in sE:
					k = len(sE)
					while sE[k - 1] == '0':
						k -= 1
					sE = sE[:k] if sE[k - 1] != '.' else sE[:k - 1]
				s += e + sE
  	return s	

def num2alpha(index, charSet='_abcdefghijklmnopqrstuvwxyz'):
	''' Transform an integer into the corresponding base-27 representation, with '_' equivalent to '0' and 'z' to the last digit.
	Used to automatically name conserved moieties. Called by computeMetForm
	'''
	N, s = len(charSet), ''
	k = index / N
	while k > 0:
		s = charSet[index - k * N] + s
		index, k = k, k / N
	s = charSet[index] + s
	return s

	
def extreme_rays_from_null_basis(A, tolerance=None):
	'''Calculate a subset of the extreme rays of the matrix A (in the format of list of list of rows)
	Return a numpy matrix.
	'''
	def rref(A, tolerance=None):
		#Find the reduced row echolon form of A. Follow the implementation in Metatool, rref.m, which is much faster than the original Matlab implementation
		#http://pinguin.biologie.uni-jena.de/bioinformatik/networks/metatool/metatool5.1/metatool5.1.html
		#A. von Kamp and S. Schuster: Metatool 5.0: fast and flexible elementary modes analysis. Bioinformatics 22 (15), 2006, 1930-1931.
		rows, cols = A.shape

		if rows == 0 or cols == 0:
			return A, [], 0

		if tolerance is None:
			tolerance = 1e-15 * max(rows, cols) * numpy.linalg.norm(A, ord=float('inf'))

		used = [False for j in range(cols)]
		r = 0
		for c in range(cols):
			#Find the pivot row
			m, pivot = abs(A[r:rows, c]).max(),  abs(A[r:rows, c]).argmax()
			pivot = r + pivot

			if m <= tolerance:
				#skip c, making sure the approximately zeros are actually zero
				A[r:rows, c] = 0
			else:
				#keep track of bound variables
				used[c] = True
				#swap current row and pivot row
				A[[pivot, r], c:cols] = A[[r, pivot], c:cols]
				#normalize pivot row
				A[r, c:cols] = A[r, c:cols] / A[r, c]
				#Eliminate the current column
				ridx = [k for k in range(rows) if k != r]
				A[ridx, c:cols] = A[ridx, c:cols] - A[ridx, c] * A[r, c:cols]

				#Check if done
				if r == rows - 1:
					break
				r = r + 1
		return A, [j for j in range(cols) if used[j]], tolerance

	if isinstance(A, (list,tuple)):
		if max([len(r) for r in A]) > min([len(r) for r in A]):
			raise ValueError, 'The rows in the matrix have different lengths.'

	#Follow the Matlab implementation. 
	R, pivcol, tolerance = rref(numpy.matrix(A), tolerance)
	r = len(pivcol) #rank
	n = R.shape[1]
	nopiv = [k for k in xrange(n) if k not in pivcol]
	N = numpy.matrix(numpy.zeros((n, n - r)))
	if n > r:
		N[nopiv, :] = numpy.eye(n - r)
		if r > 0:
			N[pivcol, :] = -R[range(r),:][:, nopiv]
	colpos = (N >= 0).all(0)
	colpos = [k for k in xrange(N.shape[1]) if colpos[0,k]]
	colneg = (N <= 0).all(0)
	colneg = [k for k in xrange(N.shape[1]) if colneg[0,k]]
	N = numpy.concatenate((N[:,colpos], -N[:,colneg]), 1)
	return N

	
