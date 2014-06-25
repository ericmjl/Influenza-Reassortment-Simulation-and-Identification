"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

class Reconstructor(object):
	"""docstring for Reconstructor"""
	def __init__(self):
		"""
		The Reconstructor class allows us the flexibility to test different 
		reconstruction algorithms.

		ATTRIBUTES:
		NOTE: I am currently debating whether these attributes are actually 
		necessary or not. Therefore, this is still an alpha version of the 
		code. 

		-	DICTIONARY: sequences
				A dictionary that holds SeqRecord objects for each of the 
				segments.

		-	LIST: graphs
				A list of the graphs that are present.
		"""
		super(Reconstructor, self).__init__()
		
		self.sequences = dict()

		self.graphs = [] 
		