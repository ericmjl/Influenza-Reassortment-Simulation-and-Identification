"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from generate_id import generate_id
from copy import copy, deepcopy
from random import choice
from numpy.random import normal

class Pathogen(object):
	"""The Pathogen class models a generic pathogen."""
	def __init__(self, segments, creation_time, progeny_size, \
		parent=dict(), convenient_id=None):
		"""
		Initialize the pathogen.

		ATTRIBUTES:
		-	LIST: segments
				a list of Segment objects that specify the genomic segments 
				that are present in the pathogen.

		-	INT: creation_time
				an integer number that specifies the time within the 
				simulation that the pathogen was created.

		- 	DICTIONARY: parent
				A dictionary that specifies the parent of the pathogen. 
				The (key, value) pairs are (parent, segments transmitted).
					-	If the dictionary is empty, then the pathogen is a seed pathogen.
					-	If the dictionary has one (key, value) pair, then the pathogen was replicated from one other pathogen.
					-	If the dictionary has two (key, value) pairs, then the pathogen was a reassortant of two parental pathogens.

		-	STRING: id
				a string generated form the generate_id function that uniquely 
				identifies the pathogen.

		-	TUPLE: progeny_size
				a two-tuple that describes the (mean, var) of the burst size 
				per pathogen. These are the Normal distribution parameters. In 
				practice, while a floating point number is drawn from the 
				Normal distribution, it will be rounded to an integer number.
		"""
		super(Pathogen, self).__init__()

		self.segments = segments

		self.creation_time = creation_time

		self.parent = parent

		self.id = generate_id()

		# self.convenient_id = self.id[0:5]

		self.progeny_size = progeny_size
		

	def __repr__(self):
		return self.id[0:5]


	def reassort_with(self, other, current_time): 
		"""
		This method takes in another pathogen and returns a progeny with 
		genomic segments that are randomly selected from segments from both
		pathogens. Reassortment with another pathogen is akin to 
		replication, but the step at which segments are added is different.

		INPUTS:
		-	OBJECT: other
				The other pathogen with which to reassort.

		-	INT: current_time
				The time at which the reassortment event happened.


		OUTPUTS:
		-	OBJECT: new_pathogen
				The reassortant progeny from the two pathogens.
		"""

		new_pathogen = copy(self)
		new_pathogen.creation_time = current_time
		new_pathogen.id = generate_id()
		new_pathogen.segments = []
		
		# Assign parent differently from the replicate() function.
		parent = dict()
		parent[self] = []
		parent[other] = []

		# Iterate over all of the segments.
		for segment in self.segments:
			i = choice([0, 1])
			segment_number = segment.segment_number
			segments_to_choose_from = (segment, \
				other.get_segment(segment_number))
			new_pathogen.segments.append(deepcopy(segments_to_choose_from[i]))
			if i == 0:
				parent[self].append(segments_to_choose_from[0].segment_number)
			if i == 1:
				parent[other].append(segments_to_choose_from[1].segment_number)

		# Clean the dictionary if the progeny is entirely derived from one 
		# parent.
		if parent[self] == []:
			del parent[self]
		if parent[other] == []:
			del parent[other]

		new_pathogen.parent = parent

		new_pathogen.mutate()

		# print new_pathogen.parent
		return new_pathogen

	def generate_progeny(self, current_time):
		"""
		This method takes calls on the replicate() function to generate n 
		progeny, drawn from the Normal distribution.

		INPUTS:
		-	INT: current_time:
				This number will be set to the "creation_time" of the new 
				pathogen, and is passed into the replicate() function.


		OUTPUTS:
		-	LIST: progeny
				A list of progeny pathogens that can be extended onto a 
				"master list" of pathogens.
		"""
		mean = self.progeny_size[0]
		var = self.progeny_size[1]
		num_progeny = normal(mean, var)

		progeny = []
		for i in range(num_progeny):
			progeny.append(self.replicate(current_time))

		return progeny


	def replicate(self, current_time): #, convenient_id=None):
		"""
		This method replicates the pathogen once. Mutation is guaranteed to be 
		called, but not guaranteed to happen, as it depends on the 
		substitution rate of each of the genomic segments.

		INPUTS:
		-	INT: current_time
				This number will be set to the "creation_time" of the new 
				pathogen.


		OUTPUTS:
		-	OBJECT: new_pathogen
				The replicated pathogen.
		"""
		new_pathogen = copy(self)
		new_pathogen.parent = dict()
		new_pathogen.parent[self] = self.get_segment_numbers()
		new_pathogen.creation_time = current_time
		new_pathogen.id = generate_id()
		new_pathogen.segments = deepcopy(self.segments)
		new_pathogen.mutate()

		# print new_pathogen.parent
		return new_pathogen

	def mutate(self):
		"""
		This method mutates each of the genomic segments according to 
		their specified substitution rates. It does not return anything.
		"""
		for segment in self.segments:
			segment.mutate()


	#################### PATHOGEN PROPERTY METHODS ############################
	
	def is_seed(self):
		"""
		This method returns a Boolean value telling us if the pathogen is 
		a "seed" pathogen. A "seed" pathogen is one that "seeded" the 
		population of pathogens present in the population.
		"""
		if len(self.parent) == 0:
			return True
		else:
			return False

	def is_reassorted(self):
		"""
		This method returns a Boolean value telling us if the pathogen is 
		a "reassorted" pathogen. A "reassorted" pathogen is one that has 
		genomic segments from two different parents.
		"""

		if len(self.parent) == 2:
			return True
		else:
			return False

	def mutations(self):
		"""
		This method prints the mutations present inside the pathogen. 
		"""
		for segment in self.segments:
			print segment.mutations

	def get_segment(self, segment_number):
		"""
		This method returns the given segment by segment number.

		INPUTS:
		-	INT: segment_number
				The integer number of the segment corresponding to the Segment 
				object's segment_number attribute

		OUTPUTS:
		-	SEGMENT OBJECT: seg
				The Segment object from the pathogen that has the 
				corresponding Segment number as specified in its 
				segment_number attribute.
		"""
		for segment in self.segments:
			if segment.segment_number == segment_number:
				return segment

	def get_segment_numbers(self):
		"""
		This method returns a list of all of the segment numbers for each of 
		the segments present in the pathogen.
		"""
		segment_numbers = []
		for segment in self.segments:
			segment_numbers.append(segment.segment_number)

		return segment_numbers

