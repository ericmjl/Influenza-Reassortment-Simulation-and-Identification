"""
Author: Eric J. Ma
Affiliation: Massachusetts Institute of Technology

Purpose of this Python definitions file:
- To provide accuracy metric functions for use in simulations.
"""
	
# Prototype #2 for Accuracy Functions (created on 2 July 2014)

def fraction_accurate_reassortants(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in the simulator and reconstructor objects, and returns
	the fraction of reassortants that are correct.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
	
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(overlaps) / len(simulated)
			The fraction of simulation reassortants that were correctly 
			identified as such in the reconstruction.
	"""
	
	simulation = [str(item) for item in simulator.reassortants()]
	reconstruction = reconstructor.reassortants(reconstruction_type=reconstruction_type)
	
	overlaps = 0
	
	for isolate in reconstruction:
		if isolate in simulation:
			overlaps += 1
	
	if len(simulation) == 0:
		return 0
	else:
		return float(overlaps) / len(simulation)

def fraction_inaccurate_reassortants(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in the list of simulation reassortants and reconstruction
	reassortants, and returns the fraction of the reconstruction reassortants
	that were incorrect.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.

	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.

	OUTPUTS:
	-   FLOAT: float(incorrect) / len(reconstructed)
			The fraction reconstruction reassortants that were not present as 
			reassortants in the simulation.
	"""
	simulation = [str(item) for item in simulator.reassortants()]
	reconstruction = reconstructor.reassortants(reconstruction_type=reconstruction_type)
		
	incorrect = 0
	for isolate in reconstruction:
		if isolate not in simulation:
			incorrect += 1
	
	if len(reconstruction) == 0:
		return 0
	else:
		return float(incorrect) / len(reconstruction)


def fraction_accurate_edges(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in a simulator and reconstructor object, and returns the
	fraction of accurate edges that were identified in the Reconstructor's 
	reconstruction network.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
			
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(overlaps) / len(simulation)
			The fraction of simulation edges that were correctly identified in 
			the specified reconstruction.
	"""
	
	simulation = simulator.relabeled_transmission_graph.edges(data=True)
	
	if reconstruction_type == 'reconstruction':
		reconstruction = reconstructor.pruned_condensed_graph.edges(data=True)
	if reconstruction_type == 'reassigned_source':
		reconstruction = reconstructor.reassigned_source_graph.edges(data=True)
		
	overlaps = 0
	for edge in reconstruction:
		if edge in simulation:
			overlaps += 1
	if len(simulation) == 0:
		return 0
	else:
		return float(overlaps) / len(simulation)

def fraction_inaccurate_edges(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in a simulator and reconstructor object, and returns the
	fraction of edges in the reconstruction that were not present in the 
	simulation.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
			
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(overlaps) / len(simulation)
			The fraction of simulation edges that were correctly identified in the specified reconstruction.
	"""
	
	simulation = simulator.relabeled_transmission_graph.edges(data=True)
	
	if reconstruction_type == 'reconstruction':
		reconstruction = reconstructor.pruned_condensed_graph.edges(data=True)
	if reconstruction_type == 'reassigned_source':
		reconstruction = reconstructor.reassigned_source_graph.edges(data=True)
		
	incorrect = 0
	
	for edge in reconstruction:
		if edge not in simulation:
			incorrect += 1

	if len(reconstruction) == 0:
		return 0
	else:
		return float(incorrect) / len(reconstruction)

def path_accuracy(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in a simulator and reconstructor object, and returns the
	fraction of edges in the reconstruction that represented a path in the 
	original simulation. This becomes especially pertinent for the case where 
	sampling occurs.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
			
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(num_correct) / float(num_considered)
			The fraction of edges in the reconstruction that represented an 
			accurate path in the simulation.
	"""
	
	# simulation = simulator.relabeled_transmission_graph.edges(data=True)
	# simulation_full = simulator.full_transmission_graph.edges()
	# simulation_reas = simulator.reassortant_edges
	
	if reconstruction_type == 'reconstruction':
		reconstruction = reconstructor.pruned_condensed_graph.edges(data=True)
	if reconstruction_type == 'reassigned_source':
		reconstruction = reconstructor.reassigned_source_graph.edges(data=True)

	num_considered = 0
	num_correct = 0
	
	for edge in reconstruction:
		num_considered += 1
		if reconstructor.is_full_transmission_edge(edge):
			if simulator.full_transmission_path_exists(edge[0], edge[1]):
				num_correct += 1
			else:
				num_correct_segments = 0
				for segment in edge[2]['segments']:
					if simulator.segment_transmission_path_exists(edge[0], edge[1], int(segment)):
						num_correct_segments += 1
						
				if num_correct_segments == len(edge[2]['segments']):
					num_correct += 1
		if not reconstructor.is_full_transmission_edge(edge):
			num_correct_segments = 0
			for segment in edge[2]['segments']:
				if simulator.segment_transmission_path_exists(edge[0], edge[1], int(segment)):
					num_correct_segments += 1
					
			if len(edge[2]['segments']) == num_correct_segments:
				num_correct += 1
				
			else:
				pass
			
	return float(num_correct) / float(len(reconstruction))

def fraction_accurate_reassortant_edges(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in the simulator and reconstructor objects, and returns
	the fraction of ground truth reassortant edges that were found in the 
	reconstruction.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
	
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(overlaps) / len(simulated)
			The fraction of simulation reassortant edges that were correctly 
			identified
			as such in the reconstruction.
	"""
	
	simulation = simulator.reassortant_edges
	reconstruction = reconstructor.reassortant_edges(reconstruction_type)

	overlaps = 0
	for edge in reconstruction:
		if edge in simulation:
			overlaps += 1

	if len(simulation) == 0:
		return 0
	else:
		return float(overlaps)/len(simulation)
	
def fraction_inaccurate_reassortant_edges(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in the simulator and reconstructor objects, and returns
	the fraction of reconstruction reassortant edges that are incorrect.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
	
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(incorrect) / len(reconstruction)
			The fraction of reconstruction reassortant edges that were 
			incorrectly identified as reassortant edges.
	"""
	
	simulation = simulator.reassortant_edges
	reconstruction = reconstructor.reassortant_edges(reconstruction_type)
	
	incorrect = 0
	
	for edge in reconstruction:
		if edge not in simulation:
			incorrect += 1

	if len(reconstruction) == 0:
		return 0
	else:
		return float(incorrect) / len(reconstruction)


def reassortant_path_accuracy(simulator, reconstructor, reconstruction_type='reconstruction'):
	"""
	This method takes in a simulator and reconstructor object, and returns the
	fraction of reassortant edges in the reconstruction that represented the 
	segment transmission path in the simulation.
	
	INPUTS:
	-   OBJECTS: simulator, reconstructor
			The simulator and reconstructor objects that hold the graphs.
			
	-   STRING: reconstruction_type
			A string specifying the type of reconstruction that we want to 
			evaluate accuracy for.
			Currently, we allow:
			-   'reconstruction': the best reconstruction possible.
			-   'reassigned_source': a shuffled version of the reconstruction, 
				in which the edges are shuffled by ignoring genetic similarity.
			
	OUTPUTS:
	-   FLOAT: float(num_correct) / float(num_considered)
			The fraction of segment transmission edges in the reconstruction 
			that represented an accurate path in the simulation.
	"""

	simulation = simulator.relabeled_transmission_graph.edges(data=True)

	reconstruction = reconstructor.reassortant_edges(reconstruction_type)

	num_considered = 0
	num_correct = 0

	for edge in reconstruction:
		for segment in edge[2]['segments']:
			num_considered += 1
			if simulator.segment_transmission_path_exists(edge[0], edge[1], int(segment)):
				num_correct += 1

	if len(reconstruction) == 0:
		return 0

	else:
		return float(num_correct) / float(num_considered)
