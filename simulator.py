"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

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

	def add_pathogens(self, pathogens)


	def increment_timestep(self):
		"""
		This is the customizable part of the simulator. Here, the actions that 
		take place each timestep are defined. 
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
			self.current_time: 0
			self.pathogens: [] #empty list
		"""

		self.current_time = 0
		self.pathogens = []