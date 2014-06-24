"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from random import choice, random, randint, sample
from sequence import Sequence
from numpy.random import binomial

class Segment(object):
	"""
	The Segment class models a "genomic segment" of a pathogen.
	"""

	def __init__(self, segment_number, substitution_rate, sequence):
		"""
		Initialize a viral segment.

		INPUTS:
		-	INT: segment_number
				an integer that identifies the segment number, for 
				identification purposes only.
		
		-	FLOAT: substitution_rate
				a floating point number, usually a small number between 0 and 
				0.1, that specifies the expecte number of mutations that occur 
				per position per year.
		
		-	OBJECT: sequence
				a Sequence object that specifies the seed sequence of the 
				genomic segment.
		"""
		super(Segment, self).__init__()
		
		self.seed_sequence = sequence
		
		self.segment_number = segment_number

		self.mutations = dict()

		self.length = len(self.compute_sequence())

		self.substitution_rate = substitution_rate

	def __repr__(self):
		return 'Segment %s, sequence %s' % (self.segment_number, \
			self.seed_sequence)

	def compute_sequence(self):
		"""
		This method computes the segment's sequence by comparing the seed 
		sequence with the mutation dictionary.
		"""
		sequence = ''

		for i, letter in enumerate(self.seed_sequence.sequence):
			if i in self.mutations.keys():
				sequence += self.mutations[i]
			else:
				sequence += letter

		return sequence

	def mutate(self):
		"""
		This method uses the length of the segment and the segment's mutation 
		rate to identify the number of positions that will be mutated. It then
		chooses that many positions at random, and records the mutation in the
		segment's mutation dictionary.
		"""
		n = self.length
		p = float(self.substitution_rate)

		num_positions = binomial(n,p)

		def choose_positions(start, end, num_positions):
			"""
			This function chooses n positions at random within
			range(start, end)

			INPUTS:
			-	INT: start
					lower bound of the range of positions to choose from
			-	INT: end
					upper bound of the range of positions to choose from
			- 	INT: num_positions
					the number of positions to be mutated

			OUTPUTS:
			-	a list of positions within the bounds (start, end)
			"""
			return sample(range(start, end), num_positions)

		positions = choose_positions(0, len(self.compute_sequence()), \
			num_positions)

		def choose_new_letter(letter):
			"""
			This function chooses a new letter from ATGC that is
			different from the letter passed into the function.

			INPUTS:
			-	CHAR: letter
					the letter that will not be chosen from ATGC.
			"""
			possible_letters = set(['A', 'T', 'G', 'C'])
			new_letter = choice(list(
				possible_letters.difference(set(letter))))

			return new_letter

		for position in positions:
			if position in self.mutations.keys():
				letter = self.mutations[position]
			else:
				letter = self.seed_sequence.sequence[position]

			self.mutations[position] = choose_new_letter(letter)

			# Note: this mutational simulation process allows the virus to 
			# back-mutate. In this case, we consider the back-mutation to
			# remain a type of "mutation", rather than a reversion, because 
			# it is different from its  "parental" sequence. 



