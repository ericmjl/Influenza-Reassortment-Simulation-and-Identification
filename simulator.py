"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology
"""

class Simulator(object):
	"""The Simulator class defines how a simulation run happens."""
	def __init__(self, arg):
		"""
		Initialize the simulator.

		ATTRIBUTES:
		-	INT: current_time
				The current time in the simulator.
				
		-	LIST: viruses
				A list of viruses currently present in the simulator.
		"""
		super(Simulator, self).__init__()
		self.arg = arg
		