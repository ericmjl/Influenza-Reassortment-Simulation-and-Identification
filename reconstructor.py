"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from Levenshtein import distance
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import networkx as nx

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
		
		self.segments = segments

		self.sequences = dict()
		self.graphs = dict()
		for segment in segments:
			self.sequences[segment] = []
			self.graphs[segment] = nx.DiGraph()

	def read_fasta_file(self, fasta_file, splitchar='|', pos_segment_num=1):
		"""
		This method takes in a FASTA file, and adds the sequences to the self.
		sequences dictionary, split by segment number.

		INPUTS:
		-	STRING: fasta_file
				A string that describes the location of the FASTA file that 
				contains the sequences of each simulation run.

		-	CHAR: splitchar
				A character that is used for splitting up the SeqRecord.id 
				attribute

		-	INT: pos_segment_num
				The position of the segment number in the id.
		"""

		sequences = [record for record in SeqIO.parse(fasta_file, 'fasta')]
		
		for sequence in sequences:
			segnum = int(sequence.id.split(split_char)[pos_segment_num])

			self.sequences[segnum].append(sequence)

	def add_nodes_with_data(self, id_attributes=None, splitchar='|'):
		"""
		This method adds in nodes with the metadata attached to it.

		INPUTS:
		-	DICT: id_attributes (optional)
				A dictionary of attributes and their positions that are stored 
				in the SeqRecord.id attribute. Currently required is the 
				following id_attributes dictionary:

					{'id':0, 
					 'segment_number':1, 
					 'creation_time':2}

				Those three fields must be guaranteed to be present. Otherwise, an Error message will be raised.

		-	CHAR: splitchar (optional)
				A character that is used for splitting up the SeqRecord.id 
				attribute
		"""
		####################### START HELPER FUNCTIONS #######################
		def reverse_dictionary(dictionary):
			"""
			This method takes in a dictionary and reverses the key, value 
			pairs. It is assumed that the keys and values are all unique.
			"""
			new_dictionary = dict()
			for k, v in dictionary.items():
				new_dictionary[v] = k

			return new_dictionary

		def add_attribute(node, attributes_dict):
			"""
			This method takes in a networkx graph's node, and adds attributes 
			to it from a dictionary of attributes. 
			"""

			for k, v in attributes.items():
				node[k] = v

		def make_attributes_dict(list_of_attributes, attribute_pos):
			"""
			This method takes in a dictionary of {positions:attribute_ids}, and returns an attribute dictionary of {attribute_ids:attribute_value}, using the list_of_attributes.
			"""
			attribute_dict = dict()
			for pos, attribute in attribute_pos.items():
				attribute_dict[attribute] = list_of_attributes[pos]

			return attribute_dict
	
		######################## END HELPER FUNCTIONS ########################

		required_attributes = {'id':0, 'segment_number':1, 'creation_time':2}

		if id_attributes = None:
			id_attributes = required_attributes

		attribute_pos = reverse_dictionary(id_attributes)

		for segment in self.sequences:
			G = self.graphs[segment]

			for sequence in self.sequences[segment]:
				# Get the list of attributes
				list_of_attributes = sequence.id.split(splitchar)

				# Make an attributes dictionary
				attribute_dict = make_attributes_dict(list_of_attributes, attribute_pos)

				# Add the node using the node name.
				node_name = attribute_dict['id'][0:5]
				G.add_node(node_name)

				# Add the attributes to the node
				node = G.node[node_name]
				add_attribute(node, attribute_dict)

				# Add the sequence to the node
				G.node[node_name]['sequence'] = str(sequence.seq)


	def add_edges_with_weight(self):
		"""
		This method adds in the edges to the graph, with the weight equal to 
		the Levenshtein distance between two sequences.
		"""

		for segment in self.graphs:
			G = self.graphs[segment]

			for node1 in G.nodes(data=True):
				for node2 in G.nodes(data=True):
					if node1[1]['creation_time'] < node2[1]['creation_time']:
						weight = distance(node1[1]['sequence'], \
							node2[1]['sequence'])
						G.add_edge(node1[0], node2[0], weight=weight, \
							segment=segment)

	def prune_graph(self):
		"""
		This method prunes the graph down, such that the only edges remaining 
		are the in_edges of minimum weight for each node.
		"""
		for segment in self.graphs:
			G = self.graphs[segment]

			for node in G.edges(data=True):
				in_edges =G.in_edges(node[0], data=True)

				if len(in_edges) != 0:
					min_weight = min([edge[2]['weight'] for edge in in_edges])

					for edge in in_edges:
						if edge[2]['weight'] != min_weight:
							G.remove_edge(edge[0], edge[1])

	def compose_segment_graphs(self):
		"""
		This method composes each of the segment transmission graphs into a 
		single MultiDiGraph.
		"""
		