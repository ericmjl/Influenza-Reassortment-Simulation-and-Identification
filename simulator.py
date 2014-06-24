"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

from Bio import Seq, SeqRecord, SeqIO
from Bio.Alphabet import IUPAC

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

	def increment_timestep(self):
		"""
		This is the customizable part of the simulator. Here, the actions that 
		take place each timestep are defined. When sub-classed, one can 
		customize different actions to be run at each time step.
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
		for i in range(timesteps):
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
			convenient_id = pathogen.convenient_id
			creation_date = pathogen.creation_date

			for segment in pathogen.segments:
				segment_number = segment.segment_number
				sequence_string = Seq(segment.compute_sequence())

				sequence = SeqRecord(sequence_string, IUPAC.nucleotide, \
					id="%s|%s|%s" % (convenient_id, segment_number, \
						creation_date))

				sequences.append(sequence)

		if folder_name == None:
			output_handle = open('%s' % outfile_name)

		else:
			output_handle = open('%s\%s' % (folder_name, outfile_name))

		SeqIO.write(output_handle, sequences, 'fasta')
		output_handle.close()

	def write_network(self, outfile_name, folder_name=None):
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
		transmission_graph = nx.DiGraph()

		for pathogen in pathogens:
			transmission_graph.add_node(pathogen)

			if virus.parent != Nonte:
				transmission_graph.add_edge(virus.parent, virus)

		if folder_name == None:
			output_handle = open('%s.gpickle' % outfile_name)

		else:
			output_handle = open('%s\%s.gpickle' % (folder_name, outfile_name))

		nx.write_gpickle(transmission_graph, output_handle)