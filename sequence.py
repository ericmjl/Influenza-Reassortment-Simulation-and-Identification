"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from random import choice
from generate_id import generate_id

class Sequence(object):
	"""
	The Sequence object is the lowest level object in the pathogen simulator. 
	It provides a container for storing seed sequences for the pathogens present
	in the environment. 

	This can be subclassed to store seed sequences for other pathogens, rather 
	than using a generated sequence. 

	Note that when a virus replicates, the full sequence object is not copied 
	for each of its segments; rather, each segment only keeps track of the 
	mutations that have happened.
	"""
	def __init__(self, length=1000, sequence=None, id=None):
		"""
		Initialize the sequence with a random sequence of length 1000 if 
		sequence is not specified.

		Otherwise, initialize sequence with a sequence that is specified. 
		"""
		super(Sequence, self).__init__()

		self.sequence = None

		if sequence == None:
			self.sequence = self.generate_sequence(length)
		else:
			self.sequence = sequence

		if id == None:
			self.id = generate_id()
		else:
			self.id = id

	def __repr__(self):
		return self.id
		
	def generate_sequence(self, length):
		"""
		This method will generate a sequence, and set the Sequence object's 
		sequence to that sequence.
		"""
		sequence = ''
		for i in range(length):
			letter = choice(['A', 'T', 'G', 'C'])
			sequence += letter
		return sequence