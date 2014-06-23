"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from generate_id import generate_id
from copy import copy, deepcopy
from random import choice

class Pathogen(object):
	"""The Pathogen class models a generic pathogen."""
	def __init__(self, segments, creation_time, parent=None, \
		convenient_id=None):
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

		-	INT or STRING: convenient_id
				a string or an integer that provides a convenient 
				representation of the virus. This item is not involved in any 
				computation. 

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
				The reassortant progeny from the two viruses.
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


	def replicate(self, current_time, convenient_id=None):
		"""
		This method replicates the pathogen. Mutation is guaranteed to be 
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
		THis method returns a Boolean value telling us if the pathogen is 
		a "reassorted" pathogen. A "reassorted" pathogen is one that has 
		genomic segments from two different parents.
		"""

		if len(self.parent) == 2:
			return True
		else:
			return False

	def mutations(self):
		for segment in self.segments:
			print segment.mutations