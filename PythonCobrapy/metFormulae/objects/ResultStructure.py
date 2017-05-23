class DataObject(object):
	def __repr__(self):
		print 'Attribute\tClass\tValue/Length'
		for k, val in self.__dict__.iteritems():
			if k[0] != '_':
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
	def __init__(self):
		self.metKnown, self.metUnknown, self.metFill, self.rxnKnown, \
		self.ele, self.eleConnect, self.metFillConnect = (None for i in range(7))

class MinInconParsiInfo(DataObject):
	def __init__(self):
		self.formulae, self.infeas, self.bound, self.obj, self.solution, self.metModel, \
		self.final, self.solStat, self.solConstrain = (None for i in range(9))

class ConservedMoietyInfo(DataObject):
	def __init__(self):
		self.cm, self.cmGeneric, self.cmGenericDict = (None for i in range(3))
  
class ResultStructure(DataObject):
	def __init__(self):
		self.formulae, self.S_fill, self.rxnBalance, self.model, self.mipInfo, self.cmInfo = (None for i in range(6))
