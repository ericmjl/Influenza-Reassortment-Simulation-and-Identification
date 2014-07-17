"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from random import randint
from Levenshtein import distance
from copy import deepcopy

import networkx as nx

class Simulator(object):
	"""The Simulator class defines how a simulation run happens."""
	def __init__(self):
		"""
		Initialize the simulator.

		ATTRIBUTES:
		-	INT: current_time
				The current time in the simulator. It is initialized to 0, and 
				can be reset to 0 using the reset() function.

		-	LIST: pathogens
				A list of pathogens currently present in the simulator. It is 
				initialized to an empty list, and can be rest to an empty list 
				using the reset() function.
		"""
		super(Simulator, self).__init__()
		
		self.current_time = 0

		self.pathogens = []

		self.node_draw_positions = None

		# All the properties that need to be pre-computed are initialized to None here first. They are computed once.
		self._reassortant_edges = None
		self._relabeled_transmission_graph = None
		self._full_transmission_graph = None
		self._full_transmission_paths = None
		self._transmission_paths = None
		self._transmission_graph = None

	def increment_timestep(self):
		"""
		This is the customizable part of the simulator. In the actual 
		simulation code, one has to subclass the Simulator and implement 
		increment_timestep. All other methods are available to the subclass.

		NOTE: self.current_time is set in the run() function, not in the 
		increment_timestep() function. Do not change that logic here.
		"""
		pass

	def run(self, timesteps):
		"""
		This method runs the simulation for as many timesteps are specified. 
		Nothing is returned.

		INPUTS:
		-	INT: timesteps
				The number of time steps that the simulation has to run.
		"""
		for i in range(1, timesteps):
			self.current_time = i
			self.increment_timestep()

	def reset(self):
		"""
		This resets the Simulator object by setting:
			-	self.current_time: 0
			-	self.pathogens: [] #empty list
		"""

		self.current_time = 0
		self.pathogens = []

	def add_pathogens(self, pathogens):
		"""
		This method takes in a list of pathogens and adds it to the current 
		list of pathogens in the simulation.

		INPUTS:
		-	ITERABLE: pathogens
				A list of Pathogen objects that will be added to self.pathogens
		"""
		self.pathogens.extend(pathogens)

	def write_sequences(self, outfile_name, folder_name=None):
		"""
		This method writes the sequences of every Pathogen object to a single 
		FASTA file.

		INPUTS:
		-	STRING: outfile_name
				The desired filename, without the ".fasta" extension.

		-	STRING: folder_name (optional)
				The folder in which the FASTA files are going to be stored
		"""

		sequences = []

		for pathogen in self.pathogens:
			creation_time = pathogen.creation_time

			for segment in pathogen.segments:
				segment_number = segment.segment_number
				sequence_string = Seq(segment.compute_sequence())


				id_tuple = (pathogen.id, segment_number, creation_time)
				sequence = SeqRecord(sequence_string, id="%s|%s|%s" % id_tuple)

				sequences.append(sequence)

		if folder_name == None:
			output_handle = open('%s.fasta' % outfile_name, 'w+')

		else:
			output_handle = open('%s/%s.fasta' % (folder_name, outfile_name), \
				'w+')

		SeqIO.write(sequences, output_handle, 'fasta')
		output_handle.close()

	@property
	def transmission_graph(self):
		"""
		This method creates the ground truth transmission graph in memory.

		NOTE: The data structure of the edges has to match the data structure 
		of the edges in the reconstruction. As of 27 June 2014, Reconstructor 
		graph edges have the following dictionary attributes:

			-	'segments' = [1, 4] (LIST: INT segment numbers)
			-	'weight' = 9 (INT Levenshtein distance between two isolates)

		NOTE: To draw to screen, call on draw_transmission_graph().

		NOTE: We have defined a few helper functions to get the sequences of 
		each segment in a pathogen.
		"""

		#################### BEGIN HELPER FUNCTIONS ###########################
		def edge_levenshtein_distance(parent, progeny, segments):
			"""
			This method computes the total Levenshtein distance over all 
			segments that were transmitted from parent to progeny. The parent 
			and progeny are Pathogen objects, and segments is a list of 
			integer numbers.
			"""	
			lev_dist = 0
			for segment in segments:
				lev_dist += segment_levdist(parent, progeny, segment)

			return lev_dist

		def segment_levdist(pathogen1, pathogen2, segment_number):
			"""
			This method returns the Levenshtein distance between two 
			pathogens' stipulated segment.
			"""
			# Get the sequence for each of the segments
			segment1 = pathogen1.get_segment(segment_number).compute_sequence()
			segment2 = pathogen2.get_segment(segment_number).compute_sequence()

			# Compute the Levenshtein distance
			lev_dist = distance(segment1, segment2)

			return lev_dist
		#################### END HELPER FUNCTIONS #############################

		#################### BEGIN MAIN LOGIC #################################
		if self._transmission_graph == None:
			transmission_graph = nx.DiGraph()

			for pathogen in self.pathogens:
				transmission_graph.add_node(pathogen, \
					creation_time=pathogen.creation_time)

				# Pass if the parent is empty - this means that the pathogen was a seed pathogen 
				if len(pathogen.parent) == 0:
					pass

				# Otherwise, add each edge with the weight.
				else:
					for parent, segments in pathogen.parent.items():
						if len(segments) != 0:
							weight = edge_levenshtein_distance(parent, pathogen, segments) 

							transmission_graph.add_edge(parent, pathogen, weight=weight, segments=segments)

			self._transmission_graph = transmission_graph

		return self._transmission_graph

		#################### END MAIN LOGIC ###################################

	def draw_transmission_graph(self, positions=False):
		"""
		This method draws the transmission graph to the screen using 
		matplotlib.

		INPUTS:
		-	BOOLEAN: positions
				If False: the circular layout will be drawn.

				If True: nodes will be restricted in the x-axis.

		"""

		# Step 1: Guarantee that transmission_graph is made.
		transmission_graph = deepcopy(self.relabeled_transmission_graph)

		# Step 2: Draw the graph according to the time-restricted layout or 
		# circular layout.
		if positions == False:
			nx.draw_circular(transmission_graph)

		if positions == True:
			if self.node_draw_positions == None:
				positions = dict()
				for pathogen in self.pathogens:
					positions[str(pathogen)] = (pathogen.creation_time, randint(0, 20))

				self.node_draw_positions = positions

			nx.draw(transmission_graph, pos=self.node_draw_positions)


	def write_transmission_graph(self, outfile_name, folder_name=None):
		"""
		This method writes the ground truth transmission network as a NetworkX 
		pickle file.

		INPUTS:
		-	STRING: outfile_name
				The desired filename, without the ".gpickle" extension.

		-	STRING: folder_name (optional)
				The folder in which the networkX gpickle files are going to be 
				stored
		"""
		transmission_graph = deepcopy(self.transmission_graph)

		if folder_name == None:
			output_handle = open('%s.gpickle' % outfile_name, 'w+')

		else:
			output_handle = open('%s/%s.gpickle' % (folder_name, outfile_name)\
			 , 'w+')

		nx.write_gpickle(transmission_graph, output_handle)

	def reassortants(self):
		"""
		This method will return the reassortant pathogens that are present in 
		the simulation. The reassortant pathogens are identifiable using their 
		is_reassortant() function.
		"""
		reassortants = []
		for pathogen in self.pathogens:
			if pathogen.is_reassorted():
				reassortants.append(pathogen)

		return reassortants
	
	@property
	def reassortant_edges(self):
		"""
		This method will return the edges that connect to reassortants as a 
		list.
		"""
		if self._reassortant_edges == None:
			edges = []

			for reassortant in self.reassortants():
				in_edges = \
				self.relabeled_transmission_graph.in_edges(str(reassortant), data=True)
				
				edges.extend(in_edges)
				# for edge in in_edges:
				# 	edges.append(edge)

			self._reassortant_edges = edges

		return self._reassortant_edges


	@property
	def relabeled_transmission_graph(self):
		"""
		This method will return a relabeled_transmission_graph, with the nodes 
		relabeled as strings rather than pathogen objects. 
		"""

		if self._relabeled_transmission_graph == None:

			# Call on.transmission_graph to guarantee that the graph is 
			# created. 
			transmission_graph = deepcopy(self.transmission_graph)

			# Create mapping from object to string
			mapping = dict()
			for node in transmission_graph.nodes():
				mapping[node] = str(node)

			# Relabel the transmission graph in a copy of the transmission graph
			relabeled = nx.relabel_nodes(transmission_graph, mapping)

			self._relabeled_transmission_graph = relabeled

		return self._relabeled_transmission_graph

	@property
	def full_transmission_graph(self):
		"""
		This method will generate the full transmission graph from the 
		relabeled transmission graph. This is done by iterating over the edges 
		present in the graph. If the progeny node in the edge is a 
		reassortant, then remove the edge. Otherwise, pass.
		"""
		if self._full_transmission_graph == None:

			# Call on relabel_transmission_graph() to guarantee that the graph is 
			# created. 
			full_graph = deepcopy(self.relabeled_transmission_graph)

			# Identify the reassortants, and then recast them as a list of 
			# strings, rather than a list of objects.
			reassortants = self.reassortants()
			reassortants = [str(item) for item in reassortants]

			for edge in full_graph.edges():
				progeny = edge[1]

				# This is the criteria for removal of an edge. If the progeny node 
				# is a reassortant, then the edge is not a full transmission edge. 
				# Therefore, the progeny node should not be in reassortant for the 
				# edge to be kept. Otherwise, the edge is removed.
				if progeny in reassortants:
					full_graph.remove_edge(edge[0], edge[1])
				else:
					pass

			self._full_transmission_graph = full_graph

		return self._full_transmission_graph

	@property
	def full_transmission_paths(self):
		"""
		This method will update and return the set of full transmission paths 
		between all nodes present in the simulation. The exact structure is a 
		disjoint set. We first create a singleton set for each node in the 
		transmission graph. We then iterate over each edge and union the sets 
		containing the nodes.

		NOTE: We are going to use the relabeled transmission graph in this 
		case, rather than the original transmission graph, as this will yield 
		a set of strings that can be compared with the reconstruction, which 
		also uses strings as node labels.

		OUTPUTS:
		-	LIST of SETS: paths
				A list of disjoint sets, for which in each set, a path exists 
				between each of the nodes.
		"""
		if self._full_transmission_paths == None:
			full_graph = self.full_transmission_graph

			paths = identify_paths(full_graph)

			self._full_transmission_paths = paths

		return self._full_transmission_paths

	def full_transmission_path_exists(self, node1, node2):
		"""
		This method will return True if a full transmission path exists 
		between two nodes.

		INPUTS:
		-	STRING: node1, node2
				The nodes that we are using for path identification are 
				strings. Therefore, node1 and node2 are strings. They are the 
				labels of the nodes present in the graph.

		OUTPUTS:
		-	BOOLEAN result that tells us if a path exists between the two 
			nodes.
		"""

		boolean = False
		paths = self.full_transmission_paths
		for path in paths:
			if node1 in path and node2 in path:
				boolean = True
				break

		return boolean


	def segment_transmission_graph(self, segment):
		"""
		This method will iterate over all of the edges in the relabeled 
		transmission graph, and if the edge's segment attribute does not 
		contian the segment specified, then the edge will be removed.

		INPUTS:
		-	INTEGER: segment 
				The segment number for which we want the segment transmission 
				graph.

		OUTPUTS:
		-	NETWORKX DIGRAPH: seg_graph
				The segment transmission graph.
		"""
		seg_graph = deepcopy(self.relabeled_transmission_graph)

		for edge in seg_graph.edges(data=True):
			if segment not in edge[2]['segments']:
				seg_graph.remove_edge(edge[0], edge[1])

		return seg_graph

	def segment_transmission_paths(self, segment):
		"""
		This method will return the segment transmission paths for a specified 
		segment.

		INPUTS:
		-	INTEGER: segment 
				An integer that specifies the segment for which the 
				transmission paths are to be found.

		OUTPUTS:
		-	LIST of SETS: paths
				A list of disjoint sets that describes which nodes are 
				connected by paths.
		"""

		seg_graph = self.segment_transmission_graph(segment=segment)

		paths = identify_paths(seg_graph)

		return paths

	def segment_transmission_path_exists(self, node1, node2, segment):
		"""
		This method will return True if a full transmission path exists 
		between two nodes.

		INPUTS:
		-	STRING: node1, node2
				The nodes that we are using for path identification are 
				strings. Therefore, node1 and node2 are strings. They are the 
				labels of the nodes present in the graph.

		-	INTEGER: segment
				The segment number for which we want to know the segment 
				transmission path.

		OUTPUTS:
		-	BOOLEAN result that tells us if a path exists between the two 
			nodes.
		"""

		boolean = False
		paths = self.segment_transmission_paths(segment=segment)
		for path in paths:
			if node1 in path and node2 in path:
				boolean = True
				break

		return boolean

#################### BEGIN HELPER METHODS FOR PATH FINDING ####################
def identify_paths(graph):
	"""
	This method takes in a graph and returns a list of sets that identify 
	which nodes have a path between them.

	INPUTS:
	-	NETWORKX GRAPH: graph
			The graph on which paths are to be found.

	OUTPUTS:
	-	LIST of SETS: paths
			A list of sets in which nodes that are within the same set have a 
			path between them.
	"""
	
	paths = []

	# Step 1: Initialize singleton sets for each of the nodes in the 
	# transmission graph.
	for node in graph.nodes():
		paths.append(set([node]))


	# Step 2: Iterate over all the edges. Find the sets that contain the 
	# two nodes, and union them.
	for edge in graph.edges():
		path1 = find_set_with_element(paths, edge[0])
		path2 = find_set_with_element(paths, edge[1])

		if path1 != path2:
			new_path = path1.union(path2)

			paths.pop(paths.index(path1))
			paths.pop(paths.index(path2))
			paths.append(new_path)

	return paths

def find_set_with_element(paths, element):
	"""
	This method will return the set that contains the specified 
	element.

	INPUTS:
	-	LIST of SETS: paths 
			A list of sets in which nodes that are within the same set have a 
			path between them.

	-	NODE: element 
			A node within a NetworkX graph.

	OUTPUTS:
	-	SET: path
			The set of nodes that are connected with each other that contains 
			the query node 'element'.
	"""
	for path in paths:
		if element in path:
			return path
#################### END HELPER METHODS FOR PATH FINDING ######################









