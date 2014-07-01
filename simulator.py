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

		-	NETWORKX DIGRAPH: transmission_graph
				A Graph that stores the simulated transmission network. This 
				network is constructed using the parent-progeny relationship 
				stored in each virus.

		-	NETWORKX DIGRAPH: relabeled_transmission_graph
				Same as transmission_graph, but the nodes are now labeled with 
				strings.

		-	LIST of SETS: full_transmission_paths
				A disjoint set data structure that stores all of the full 
				transmission paths (involving all segments). 
		"""
		super(Simulator, self).__init__()
		
		self.current_time = 0

		self.pathogens = []

		self.transmission_graph = nx.DiGraph()

		self.relabeled_transmission_graph = nx.DiGraph()

		self.full_transmission_paths = []

		self.segment_paths = dict()

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

	def generate_transmission_graph(self):
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

		transmission_graph = nx.DiGraph()

		for pathogen in self.pathogens:
			transmission_graph.add_node(pathogen, \
				creation_time=pathogen.creation_time)

			# Pass if the parent is empty - this means that the pathogen was a 
			# seed pathogen 
			if len(pathogen.parent) == 0:
				pass

			else:
				for parent, segments in pathogen.parent.items():
					if len(segments) != 0:
						weight = edge_levenshtein_distance(parent, pathogen, \
							segments) 

						transmission_graph.add_edge(parent, pathogen, weight=weight, segments=segments)

			# Add an edge comprising of all segments if the progeny pathogen 
			# was replicated from a single parent pathogen.
			# if len(pathogen.parent) == 1:
			# 	weight = isolate_levdist(pathogen.parent, pathogen)
			# 	segments = pathogen.get_segment_numbers()
			# 	transmission_graph.add_edge(pathogen.parent[0], pathogen, \
			# 		weight=weight, segments=segments)

			# Add an edge comprising the segment that was transmitted for each 
			# parent pathogen.
			# if len(pathogen.parent) == 2:
			# 	transmission_graph.add_edge(pathogen.parent[0], pathogen, \
			# 		weight=0.5)
			# 	transmission_graph.add_edge(pathogen.parent[1], pathogen, \
			# 		weight=0.5)

		self.transmission_graph = transmission_graph

	def draw_transmission_graph(self, positions=False):
		"""
		This method draws the transmission graph to the screen using 
		matplotlib.

		INPUTS:
		-	BOOLEAN: positions
				If False: the circular layout will be drawn.

				If True: nodes will be restricted in the x-axis.

		"""
		transmission_graph = self.transmission_graph

		if positions == False:
			nx.draw_circular(transmission_graph)

		if positions == True:
			positions = dict()
			for pathogen in self.pathogens:
				positions[pathogen] = (pathogen.creation_time, randint(0, 20))

			nx.draw(transmission_graph, pos=positions)

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
		transmission_graph = self.generate_transmission_graph()

		if folder_name == None:
			output_handle = open('%s.gpickle' % outfile_name, 'w+')

		else:
			output_handle = open('%s/%s.gpickle' % (folder_name, outfile_name)\
			 , 'w+')

		nx.write_gpickle(transmission_graph, output_handle)

	def identify_reassortants(self):
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

	def relabel_transmission_graph(self):
		"""
		This method will assign the self.relabeled_transmission_graph 
		attribute with the nodes relabeled as strings rather than pathogen 
		objects. 
		"""
		# Create mapping from object to string
		mapping = dict()
		for node in self.transmission_graph.nodes():
			mapping[node] = str(node)

		# Relabel the transmission graph in a copy of the transmission graph
		relabeled_transmission_graph = nx.relabel_nodes(\
			self.transmission_graph, mapping)

		# Assign relabeled graph to self.relabeled_transmission_graph
		self.relabeled_transmission_graph = relabeled_transmission_graph
		

	def identify_full_transmission_paths(self):
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
		"""

		# Step 1: Initialize singleton sets for each of the nodes in the 
		# transmission graph.
		paths = []
		for node in self.relabeled_transmission_graph.nodes():
			paths.append(set([node]))

		

		# Step 2: Iterate over all the edges. Find the sets that contain the 
		# two nodes, and union them.
		for edge in self.relabeled_transmission_graph.edges(data=True):
			# This is the criteria that specifiees that the transmission paths 
			# are "full".
			if len(edge[2]['segments']) == 2:
				path1 = find_set_with_element(paths, edge[0])
				path2 = find_set_with_element(paths, edge[1])

				if path1 != path2:
					new_path = path1.union(path2)

					paths.pop(paths.index(path1))
					paths.pop(paths.index(path2))
					paths.append(new_path)

		# Step 3: Set the full_transmission_paths attribute to be the list of 
		# paths that are present in the graph.
		self.full_transmission_paths = paths

		return self.full_transmission_paths

	def transmission_path_exists(self, node1, node2):
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
		paths = self.identify_full_transmission_paths()
		for path in paths:
			if node1 in path and node2 in path:
				boolean = True
				break

		return boolean

	def identify_segment_paths(self):
		"""
		This method will identify the transmission paths for each segment. 
		"""

		segment_paths = dict()

		# Step 1: Initialize the dictionary with such that the keys are the 
		# segment numbers of the pathogens being considered.
		segment_numbers = set()
		for pathogen in self.transmission_graph.nodes():
			for segment in pathogen.segments:
				segment_numbers.add(segment.segment_number)

		for number in segment_numbers:
			segment_paths[number] = []

		# Step 2: For each segment, add nodes to the graph.
		for node in self.relabeled_transmission_graph.nodes():
			for number in segment_numbers:
				segment_numbers[number].append(set([node]))

		# Step 3: Union-Find all of the paths.
		for edge in self.relabeled_transmission_graph.edges(data=True):
			# path1 = 
			pass











def find_set_with_element(paths, element):
	"""
	This method will return the set that contains the specified 
	element.
	"""
	for path in paths:
		if element in path:
			return path












