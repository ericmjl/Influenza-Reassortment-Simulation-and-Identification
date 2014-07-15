"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from Levenshtein import distance
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from copy import copy
from numpy.random import binomial
from random import sample
from copy import deepcopy

import networkx as nx

class Reconstructor(object):
	"""
	The Reconstructor class holds the composable methods for reconstructing 
	networks from genetic information.
	"""
	def __init__(self, segments, sampling=False):
		"""
		The Reconstructor class allows us the flexibility to test different 
		reconstruction algorithms.

		INPUTS:
		-	LIST: segments 
				A list of segment numbers. Currently, I assume that each 
				segment has an integer number associated with it. Therefore, 
				the elements in this list should be integer numbers.

		-	BOOLEAN: sampling
				A boolean to tell us whether to do sampling from simulated 
				sequences or not. 

				-	TRUE: sample only 1/2 of the isolates in reconstruction.
				-	FALSE: use all of the isolates in reconstruction.

		ATTRIBUTES:
		NOTE: I am currently debating whether these attributes are actually 
		necessary or not.

		-	DICTIONARY: sequences
				A dictionary that holds SeqRecord objects for each of the 
				segments.

		-	LIST: graphs
				A list of the graphs that are present.

		-	NETWORKX GRAPH: _reassigned_source_graph
				The reassigned source graph that serves as the Null Model 
				reconstruction for statistics.
		"""
		super(Reconstructor, self).__init__()
		
		self.segments = segments

		self.sequences = dict()
		self.graphs = dict()
		for segment in segments:
			self.sequences[segment] = []
			self.graphs[segment] = nx.DiGraph()

		self.sampling = sampling

		# Property methods' corresponding attributes are stored here.
		self._reassigned_source_graph = None
		self._condensed_graph = None
		self._pruned_condensed_graph = None

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

		#################### BEGIN HELPER FUNCTIONS ###########################

		def get_isolate_ids(sequences):
			"""
			This method will get the ids from the FASTA files.
			"""
			ids = set()

			for sequence in sequences:
				ids.add(sequence.id.split(splitchar)[0])

			return ids

		def sample_isolates(ids):
			"""
			This method will sample approximately half of the isolates at 
			random.
			"""
			n = len(ids)
			p = 0.5

			num_sampled = binomial(n, p)

			sampled_isolates = sample(ids, num_sampled)

			return sampled_isolates

		#################### END HELPER FUNCTIONS #############################

		sequences = [record for record in SeqIO.parse(fasta_file, 'fasta')]

		if self.sampling == True:
			ids = get_isolate_ids(sequences)
			sampled_isolates = sample_isolates(ids)

		if self.sampling == False:
			sampled_isolates = get_isolate_ids(sequences)

		for sequence in sequences:
			segnum = int(sequence.id.split(splitchar)[pos_segment_num])

			seqid = sequence.id.split(splitchar)[0]

			if seqid in sampled_isolates:
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

				Those three fields must be guaranteed to be present. 
				Otherwise, an Error message will be raised.

		-	CHAR: splitchar (optional)
				A character that is used for splitting up the SeqRecord.id 
				attribute

		NOTE: Helper functions are defined here to aid in readability and keep 
		the code flat where the main logic takes place, (mostly) in accordance 
		with the Zen of Python.
		"""
		#################### START HELPER FUNCTIONS ###########################
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

			for k, v in attributes_dict.items():
				node[k] = v

		def make_attributes_dict(list_of_attributes, attribute_pos):
			"""
			This method takes in a dictionary of {positions:attribute_ids}, 
			and returns an attribute dictionary of 
			{attribute_ids:attribute_value}, using the list_of_attributes.
			"""
			attribute_dict = dict()
			for pos, attribute in attribute_pos.items():
				attribute_dict[attribute] = list_of_attributes[pos]

			return attribute_dict
		#################### END HELPER FUNCTIONS #############################

		#################### BEGIN IMPORTANT LOGIC ############################
		required_attributes = {'id':0, 'segment_number':1, 'creation_time':2}

		if id_attributes == None:
			id_attributes = required_attributes

		attribute_pos = reverse_dictionary(id_attributes)

		for segment in self.sequences:
			G = self.graphs[segment]

			for sequence in self.sequences[segment]:
				# Get the list of attributes
				list_of_attributes = sequence.id.split(splitchar)

				# Make an attributes dictionary
				attribute_dict = make_attributes_dict(list_of_attributes, \
					attribute_pos)

				# Add the node using the node name.
				node_name = attribute_dict['id'][0:5]
				G.add_node(node_name)

				# Add the attributes to the node
				node = G.node[node_name]
				add_attribute(node, attribute_dict)

				# Add the sequence to the node

				G.node[node_name]['sequence_%s' % str(segment)] \
				= str(sequence.seq)
		#################### END IMPORTANT LOGIC ##############################


	def add_edges_with_weight(self):
		"""
		This method adds in the edges to the graph, with the weight equal to 
		the Levenshtein distance between two sequences.
		"""

		self.add_nodes_with_data()

		for segment in self.graphs:
			G = self.graphs[segment]

			sequence = 'sequence_%s' % str(segment)

			for source in G.nodes(data=True):
				for sink in G.nodes(data=True):
					if int(source[1]['creation_time']) < \
					int(sink[1]['creation_time']):
						weight = distance(source[1][sequence], \
							sink[1][sequence])
						G.add_edge(source[0], sink[0], weight=weight, \
							segment=segment)
					if source[1]['creation_time'] >= sink[1]['creation_time']:
						pass

	def prune_graphs_by_weight(self):
		"""
		This method prunes the graph down, such that the only edges remaining 
		are the in_edges of minimum weight for each node.
		"""
		for segment in self.graphs:
			G = self.graphs[segment]

			for node in G.nodes(data=True):
				in_edges = G.in_edges(node[0], data=True)

				if len(in_edges) != 0:
					min_weight = min([edge[2]['weight'] for edge in in_edges])

					for edge in in_edges:
						if edge[2]['weight'] != min_weight:
							G.remove_edge(edge[0], edge[1])

	def composed_graph(self):
		"""
		This method composes each of the segment transmission graphs into a 
		single MultiDiGraph.

		OUTPUTS:
		-	NETWORKX MULTIDIGRAPH: composed 
				The composed graph.
		"""

		composed = nx.MultiDiGraph()

		for segment in self.graphs:
			G = self.graphs[segment]

			composed.add_nodes_from(G.nodes(data=True))
			composed.add_edges_from(G.edges(data=True))

		return composed

	@property
	def condensed_graph(self):
		"""
		This method condenses the composed segment graphs into a single 
		DiGraph that keeps track of the number of segments transmitted in each 
		edge, while also summing up the weights.

		OUTPUTS:
		-	NETWORKX DIGRPH: condensed
				A DiGraph constructed from the composed MultiDiGraph, in which 
				each edge captures the segments that have been transmitted 
				between two viruses.
		"""

		if self._condensed_graph == None:
			composed = self.composed_graph()

			# Initialize the new graph
			condensed = nx.DiGraph()
			condensed.add_nodes_from(composed.nodes(data=True))

			# Initialize the edges dictionary
			edges = dict()
			for edge in composed.edges(data=True):
				nodes = (edge[0], edge[1])
				edges[nodes] = dict()
				edges[nodes]['segments'] = []
				edges[nodes]['weight'] = 0

			# Condense the segments into a list, and sum up the weights.
			for edge in composed.edges(data=True):
				nodes = (edge[0], edge[1])
				edges[nodes]['segments'].append(edge[2]['segment'])
				edges[nodes]['weight'] += edge[2]['weight']

			# Add the edges into the comdensed_graph.
			for nodes, attributes in edges.items():
				condensed.add_edge(nodes[0], nodes[1], attr_dict=attributes)

			self._condensed_graph = condensed

		return self._condensed_graph

	@property
	def pruned_condensed_graph(self):
		"""
		This method prunes the condensed graph such that if a node has a "full 
		transmission" edge into it, we will remove edges that have fewer than 
		full transmissions going into it. This is the graph that will be used 
		for assessing accuracy under no sampling constraints.

		To keep the code logic readable and compact, two helper functions have 
		been defined.

		OUTPUTS:
		-	NETWORKX DIGRAPH: pruned 
				A DiGraph of the pruned condensed graph, in which the full 
				transmission edges are kept preferentially over the partial 
				(non-full) transmission edges. 
		"""
		#################### BEGIN HELPER FUNCTIONS ###########################
		def has_at_least_one_full_transmission(in_edges):
			"""
			This method checks whether full transmission edges exist amongst 
			all of the in_edges for a given node. If such an edge exists, 
			return True, otherwise return False.
			"""
			boolean = False
			for edge in in_edges:
				if set(edge[2]['segments']) == set(self.segments):
					boolean = True
					break
			return boolean

		def remove_non_full_transmission(graph, in_edges):
			"""
			This method removes any edges that do not involve all segments.
			"""
			for edge in in_edges:
				if set(edge[2]['segments']) != set(self.segments):
					graph.remove_edge(edge[0], edge[1])
		#################### END HELPER FUNCTIONS #############################

		#################### BEGIN IMPORTANT LOGIC ############################
		
		if self._pruned_condensed_graph == None:
			pruned = deepcopy(self.condensed_graph)

			for node in pruned.nodes(data=True):
				in_edges = pruned.in_edges(node[0], data=True)
				if has_at_least_one_full_transmission(in_edges): 
					remove_non_full_transmission(pruned, in_edges)

			self._pruned_condensed_graph = pruned

		return self._pruned_condensed_graph
		#################### END IMPORTANT LOGIC ##############################

	def is_full_transmission_edge(self, edge):
		"""
		This method will tell us whether an edge is a full transmission edge 
		or not.

		OUTPUTS:
		-	BOOLEAN: whether the edge is a full transmission edge or not.
		"""
		if edge[2]['segments'] == self.segments:
			return True
		else:
			return False

	def reassortants(self, reconstruction_type='reconstruction'):
		"""
		This method identifies the reassortant viruses that are present in the 
		reconstruction. 

		The reassortant viruses are the viruses that do not have full 
		transmissions going into it. We use the helper function defined at the 
		bottom of this file.

		INPUTS:
		-	STRING: reconstruction_type
				A string specifying the type of reconstructin that we want to 
				evaluate accuracy for. Currently, we allow:
					-	'reconstruction': the best reconstruction possible
					-	'reassigned_source': a shuffled version of the 
						reconstruction, in which the sources for each node 
						that has an in_edge is reassigned at random.

		OUTPUTS:
		-	LIST: reassortants
				A list of reassortants present in the graph.
		"""
		#################### BEGIN HELPER FUNCITON ############################
		def has_no_full_transmissions(in_edges):
			"""
			In contrast to the has_at_least_one_full_transmission() method 
			defined above, this function returns True if there are no full 
			transmissions present. As usual, this helper function is defined 
			to keep the main code readable and compact.

			OUTPUTS:
			-	BOOLEAN: boolean 
					A boolean variable that tells us whether a full 
					transmission exists within the set of the given in_edges.
			"""
			boolean = True
			for edge in in_edges:
				if set(edge[2]['segments']) == set(self.segments):
					boolean = False
					break
			return boolean
		#################### END HELPER FUNCITON ##############################

		reassortants = []

		if reconstruction_type == "reconstruction":
			graph = self.pruned_condensed_graph
		if reconstruction_type == "reassigned_source":
			graph = self.reassigned_source_graph

		for node in graph.nodes(data=True):
			in_edges = graph.in_edges(node[0], data=True)

			if len(in_edges) > 0 and has_no_full_transmissions(in_edges) == True:
				reassortants.append(node[0])

		return reassortants

	def reassortant_edges(self, reconstruction_type='reconstruction'):
		"""
		This method identifies the edges that are associated with the 
		reassortant viruses present in the reconstructions.

		INPUTS:
		-	STRING: reconstruction_type
				A string specifying the type of reconstructin that we want to 
				evaluate accuracy for. See reassortants() for the description.

		OUTPUTS:
		-	LIST: edges
				A list of edges that are associated with each of the 
				reassortant viruses.
		"""

		edges = []

		reassortants = self.reassortants(reconstruction_type)

		if reconstruction_type == "reconstruction":
			graph = self.pruned_condensed_graph
		if reconstruction_type == "reassigned_source":
			graph = self.reassigned_source_graph

		for node in reassortants:
			in_edges = graph.in_edges(node, data=True)
			edges.extend(in_edges)

		return edges

	@property
	def full_transmission_graph(self):
		"""
		This method will return only the full transmissions in a copy of the  
		condensed graph.
		"""
		full_graph = deepcopy(self.pruned_condensed_graph)
		for edge in full_graph.edges(data=True):
			if set(edge[2]['segments']) != set(self.segments):
				full_graph.remove_edge(edge[0], edge[1])

		return full_graph

	def segment_transmission_graph(self, segment):
		"""
		This method will return only the sepcified segment transmission graph 
		in the condensed graph. The way this is done is by removing edges that 
		do not contain the specified segment.
		"""

		seg_graph = deepcopy(self.pruned_condensed_graph)

		for edge in seg_graph.edges(data=True):
			if segment not in edge[2]['segments']:
				seg_graph.remove_edge(edge[0], edge[1])

		return seg_graph

	@property
	def reassigned_source_graph(self):
		"""
		This method takes in a graph, for which the edges are permuted. The 
		specific implementation here is that we iterate over all the edges, 
		and randomly select another source node that occurs before the sink 
		node.

		INPUTS:
		-	NETWORKX GRAPH: graph
				The graph on which the sources are to be randomly reassigned.

		OUTPUTS: 
		-	NETWORKX GRAPH: reassigned
				The graph where the sources are randomly reassigned.
		"""

		if self._reassigned_source_graph == None:

			new_graph = nx.DiGraph()

			edges = self.pruned_condensed_graph.edges(data=True)
			nodes = self.pruned_condensed_graph.nodes(data=True)

			#################### BEGIN HELPER FUNCTIONS #######################
			def get_new_source(graph, node):
				"""
				This method takes in a list of nodes and a sink node, and 
				returns a new source node in which the new node is chosen from 
				the list of nodes that have a creation_time prior to the sink 
				node.
				"""
				from random import choice 

				node_data = graph.node[node]

				node_list = graph.nodes(data=True)
				candidate_sources = []

				for candidate in node_list:
					if int(candidate[1]['creation_time']) < int(node_data['creation_time']):
						candidate_sources.append(candidate)

				new_source = choice(candidate_sources)
				return new_source[0] # return the node label only

			#################### END HELPER FUNCTIONS #########################

			for edge in edges:
				sink_node = edge[1]

				while True:
					new_source = get_new_source(self.pruned_condensed_graph, sink_node)
					new_edge = (new_source, sink_node)

					if new_edge not in new_graph.edges():
						break

				new_graph.add_edge(new_source, edge[1], attr_dict=edge[2])

			self._reassigned_source_graph = new_graph

		return self._reassigned_source_graph









