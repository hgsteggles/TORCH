import numpy as np

class CornishData:
	def __init__(self, iden):
		self.dirname = "data/cornish_" + str(iden)
		self.star_data = np.genfromtxt(self.dirname + "/starpop-hm-w-a.txt", skip_header=1)
		self.exclusion_set = np.zeros(len(self.star_data[:,0]) + 1, dtype=bool)

		try:
			exclusion_list = np.fromfile(self.dirname + "/exclusion-list.txt", dtype=int, sep='\n')
		except IOError:
			exclusion_list = np.zeros(0)

		for i in exclusion_list:
			self.exclusion_set[i] = True

	def star_id_str(self, i):
		return "%04d" % (i,)
