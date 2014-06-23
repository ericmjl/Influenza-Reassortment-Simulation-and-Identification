from random import random
from datetime import datetime

import hashlib

def generate_id():
	"""
	This method sets the ID of the virus to be a unique string based on
	the string representation of the current time and a randomly chosen
	number.
	"""
	random_number = str(random())
	current_time = str(datetime.now())

	unique_string = random_number + current_time

	unique_id = hashlib.new('sha512')
	unique_id.update(unique_string)
	
	id = unique_id.hexdigest()

	return id