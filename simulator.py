"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import networkx as nx
from random import randint

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
				initialized to an empty list, and  
		"""
		super(Simulator, self).__init__()
		
		self.current_time = 0

		self.pathogens = []

		self.transmission_graph = nx.DiGraph()

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

		NOTE: to draw to screen, call on draw_transmission_graph().
		"""

		transmission_graph = nx.DiGraph()

		for pathogen in self.pathogens:
			transmission_graph.add_node(pathogen, \
				creation_time=pathogen.creation_time)

			if len(pathogen.parent) == 0:
				pass
			if len(pathogen.parent) == 1:
				transmission_graph.add_edge(pathogen.parent[0], pathogen, \
					weight=1.0)

			if len(pathogen.parent) == 2:
				transmission_graph.add_edge(pathogen.parent[0], pathogen, \
					weight=0.5)
				transmission_graph.add_edge(pathogen.parent[1], pathogen, \
					weight=0.5)

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