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
		parent=None, convenient_id=None):
		"""
		Initialize the pathogen.

		ATTRIBUTES:
		-	LIST: segments
				a list of Segment objects that specify the genomic segments 
				that are present in the pathogen.

		-	INT: creation_time
				an integer number that specifies the time within the 
				simulation that the pathogen was created.

		- 	OBJECT: parent
				a Pathogen object that specifies the parent of the pathogen. 
				If the pathogen was created ab initio, then the default value 
				of Parent is None. If the pathogen is reassorted, then parent 
				is a two-tuple of two Pathogen objects.

		-	STRING: id
				a string generated form the generate_id function that uniquely 
				identifies the pathogen.

		-	TUPLE: progeny_size
				a two-tuple that describes the (mean, var) of the burst size 
				per pathogen. These are the Normal distribution parameters. In 
				practice, while a floating point number is drawn from the 
				Normal distribution, it will be rounded to an integer number.

		-	INT or STRING: convenient_id
				a string or an integer that provides a convenient 
				representation of the pathogen. This item is not involved in 
				any computation. 

				It is best practice to keep convenient_id to a short element, 
				such as a string of less than 5 characters, or an integer of 
				less than 5 digits in length.
		"""
		super(Pathogen, self).__init__()

		self.segments = segments

		self.creation_time = creation_time

		self.parent = parent

		self.id = generate_id()

		self.convenient_id = convenient_id

		self.progeny_size = progeny_size

	def __repr__(self):
		if self.convenient_id == None:
			return str(self.id)
		else:
			return str(self.convenient_id)

	def reassort_with(self, other_pathogen, current_time, convenient_id=None):
		"""
		This method takes in another pathogen and returns a progeny with 
		genomic segments that are randomly selected from segments from both
		pathogens. Reassortment with another pathogen is akin to 
		replication, but the step at which segments are added is different.

		INPUTS:
		-	OBJECT: other_pathogen
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
		new_pathogen.convenient_id = convenient_id
		new_pathogen.segments = []
		
		parent = set()
		for segment in zip(self.segments, other_pathogen.segments):
			i = choice([0, 1])
			new_pathogen.segments.append(deepcopy(segment[i]))

			if i == 0:
				parent.add(self)
			elif i == 1:
				parent.add(other_pathogen)
		new_pathogen.parent = tuple(parent)

		new_pathogen.mutate()

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


	def replicate(self, current_time, convenient_id=None):
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
		new_pathogen.parent = self
		new_pathogen.creation_time = current_time
		new_pathogen.convenient_id = convenient_id
		new_pathogen.id = generate_id()
		new_pathogen.segments = deepcopy(self.segments)
		new_pathogen.mutate()

		return new_pathogen

	def mutate(self):
		"""
		This method mutates each of the genomic segments according to 
		their specified substitution rates. It does not return anything.
		"""
		for segment in self.segments:
			segment.mutate()


	### PATHOGEN PROPERTY METHODS ###
	
	def is_seed(self):
		"""
		This method returns a Boolean value telling us if the pathogen is 
		a "seed" pathogen. A "seed" pathogen is one that "seeded" the 
		population of pathogens present in the population.
		"""
		if self.parent == None:
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