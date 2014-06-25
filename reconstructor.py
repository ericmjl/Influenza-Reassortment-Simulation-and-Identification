"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from Levenshtein import distance
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Reconstructor(object):
	"""
	The Reconstructor class holds the composable methods for reconstructing 
	networks from genetic information."""
	def __init__(self, segments):
		"""
		The Reconstructor class allows us the flexibility to test different 
		reconstruction algorithms.

		INPUTS:
		-	LIST: segments 
				A list of segment numbers. Currently, I assume that each 
				segment has an integer number associated with it. Therefore, 
				the elements in this list should be integer numbers.

		ATTRIBUTES:
		NOTE: I am currently debating whether these attributes are actually 
		necessary or not.

		-	DICTIONARY: sequences
				A dictionary that holds SeqRecord objects for each of the 
				segments.

		-	LIST: graphs
				A list of the graphs that are present.
		"""
		super(Reconstructor, self).__init__()
		
		self.sequences = dict()
		for segment in segments:
			self.sequences[segment] = []

		self.graphs = dict()

	def read_fasta_file(self, fasta_file):
		"""
		This method takes in a FASTA file, and adds the sequences to the self.
		sequences dictionary, split by segment number.

		INPUTS:
		-	STRING: fasta_file
				A string that describes the location of the FASTA file that 
				contains the sequences of each simulation run.
		"""


		

	def add_nodes_with_data(self, sequences):
		"""
		This method adds in nodes with the metadata attached to it.

		INPUTS:
		-	LIST: sequences
				A list of BioPython SeqRecord objects. Currently, I am 
				assuming that SeqRecord.id is of the format:

					virus_id|Segment_segmentnumber|Time_creationdate

				This can then be parsed by splitting the id by the '|' 
				character, and then by the '_' character.

				We are also assuming that the list of sequences is already 
				split by segments.
		"""

		pass 