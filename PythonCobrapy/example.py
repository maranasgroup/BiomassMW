#example
import cobra
import metFormulae

model = cobra.io.load_matlab_model('B_theta_iAH991_before.mat')
met = metFormulae.MetFormulae(model)
met.computeMetForm(findCM='null',nameCM = 2)
cm = met.metFormResults.cmInfo.cmGenericDict
for cmName, cmDict in  cm.iteritems():
	print 'Conserved moiety: %s' %cmName.formula
	for i in cmDict:
		print '    %s\t%s' %(i.id, met.metFormResults.formulae[i])
print ' '
model2 = met.metFormResults.model
bm = model2.metabolites[model2.metabolites.index('biomass[e]')]
bm0 = met.model.metabolites[met.model.metabolites.index('biomass[e]')]
print 'Biomass formula: %s\nCharge: %.4f\nBiomass weight: %.4f g/mol' %(bm.formula, bm.charge, met.metFormResults.formulae[bm0].mw)

