import math
import numpy as np
from HMM_transitions import HMM_transitions


class pHMM_Viterbi:

	def __init__(self):
		self.dictionary_of_transition_paths = {'MM': 0, 'MI': 1, 'MD': 2, 'IM': 3, 'II': 4, 'ID': 5, 'DM': 6, 'DI': 7,
											   'DD': 8}

		return

	# The transision matrix is assumed to be log10 with pseudocounts
	# The format of the transition matrix is assumed to be the 9 transition possibilities
	# i.e. MM, MI, MD, IM, II, etc as 9 rows, and the matched states M0 to the maximum M
	# as columns

	def calculate_Viterbi(self, sequence, number_of_matched_states, transition_matrix, emission_matrix,
						  symbol_dictionary):

		length_of_sequence = len(sequence)

		# set up the Viterbi matrix as a [sequence_length+1] number of rows x [number_of_matched_states+2] number of columns matrix
		# each cell contain a tupple x,y,z (initialised to (0,0,0)) representing the M, I and D transitions
		#
		# row 0 [column 0[(VM,VI,VD)], column 1[(VM,VI,VD)], column 2[(VM,VI,VD)], etc]
		# row 1 etc
		#
		#

		viterbi_matrix = []
		direction_of_movement_matrix = []
		for row in range(length_of_sequence + 1):
			row = []
			direction_row = []
			for column in range(number_of_matched_states + 1):
				entry = []
				direction_entry = []
				for state in range(3):
					entry.append('*')
					direction_entry.append('*')
				row.append(entry)
				direction_row.append(direction_entry)
			viterbi_matrix.append(row)
			direction_of_movement_matrix.append(direction_row)
		viterbi_matrix[0][0][0] = 0

		self.fill_first_row(viterbi_matrix, direction_of_movement_matrix, transition_matrix, number_of_matched_states)
		self.fill_first_column(viterbi_matrix, direction_of_movement_matrix, transition_matrix)
		self.complete_viterbi_matrix(sequence, viterbi_matrix, direction_of_movement_matrix, transition_matrix, emission_matrix, symbol_dictionary)

		return viterbi_matrix, direction_of_movement_matrix

	def fill_first_row(self, viterbi_matrix, direction_of_movement_matrix, transition_matrix, number_of_columns):
		for column in range(1,
							number_of_columns + 1):  # we start at 1 because state 0 has no delete state, i.e., it remains undefined
			viterbi_matrix[0][column][2] = self.calculate_VD(viterbi_matrix, direction_of_movement_matrix, 0, column, transition_matrix)
		return

	def fill_first_column(self, viterbi_matrix, direction_of_movement_matrix, transition_matrix):
		number_of_rows = len(viterbi_matrix)
		for row in range(1,
						 number_of_rows):  # we start at 1 because state 0 has no delete state, i.e., it remains undefined
			viterbi_matrix[row][0][1] = self.calculate_VI(viterbi_matrix, direction_of_movement_matrix, row, 0, transition_matrix)
		return

	def complete_viterbi_matrix(self, sequence, viterbi_matrix, direction_of_movement_matrix, transition_matrix, emission_matrix, symbol_dictionary):
		number_of_columns = len(viterbi_matrix[0])
		number_of_rows = len(viterbi_matrix)

		for row in range(1, number_of_rows):
			for column in range(1, number_of_columns):
				viterbi_matrix[row][column][0] = self.calculate_VM(sequence, viterbi_matrix, direction_of_movement_matrix, row, column, transition_matrix, emission_matrix,symbol_dictionary)

				# DO the D state before the I state, because in state 1 you can calculate state D1 using I0 data
				# and need D1 to calculate I1.
				viterbi_matrix[row][column][2] = self.calculate_VD(viterbi_matrix, direction_of_movement_matrix, row, column, transition_matrix)
				viterbi_matrix[row][column][1] = self.calculate_VI(viterbi_matrix, direction_of_movement_matrix, row, column, transition_matrix)
		return

	def calculate_VD(self, viterbi_matrix, direction_of_movement_matrix, row, column, transition_matrix):

		# If no previous probability is defined, we want to make sure that
		# this parameter will not be the maximum value in the new probability,
		# so set it to a big negative value

		big_negative_number = -1000

		# M_I_D_parameters is [Vmd,Vid,Vdd] as M_I_D_parameters at index 0, 1 and 2
		M_I_D_parameters = [big_negative_number,big_negative_number,big_negative_number]
		key = ['MD', 'ID', 'DD']
		for parameter_index in range(3):
			if viterbi_matrix[row][column - 1][parameter_index] != '*':
				M_I_D_parameters[parameter_index] = transition_matrix[self.dictionary_of_transition_paths[key[parameter_index]]][column - 1] + viterbi_matrix[row][column - 1][parameter_index]
		# Make a tupel of the 3 parameter values

		result = (M_I_D_parameters[0], M_I_D_parameters[1], M_I_D_parameters[2])

		# determine the direction D, V or H depending on whether the M, I or D derived
		# parameter is the maximum
		direction_of_movement_matrix[row][column][2] = ['D','V','H',][result.index(max(result))]

		# return the maximum of the result tuple as VD
		max_result = max(result)

		return max_result

	def calculate_VI(self, viterbi_matrix, direction_of_movement_matrix, row, column, transition_matrix):

		big_negative_number = -1000
		key = ['MI', 'II', 'DI']
		M_I_D_parameters = [big_negative_number, big_negative_number, big_negative_number]
		#print('calculate_VI', 'row=', row, 'column=', column)
		for parameter_index in range(3):
			if viterbi_matrix[row - 1][column][parameter_index] != '*':
				if parameter_index != 1:
					p = viterbi_matrix[row-1][column][parameter_index]
				else:
					p = 0
				M_I_D_parameters[parameter_index] = transition_matrix[self.dictionary_of_transition_paths[key[parameter_index]]][column] + p
		result = (M_I_D_parameters[0], M_I_D_parameters[1], M_I_D_parameters[2])
		direction_of_movement_matrix[row][column][1] = ['D','V','H',][result.index(max(result))]
		max_result = max(result)

		return max_result

	def calculate_VM(self, sequence, viterbi_matrix, direction_of_movement_matrix, row, column, transition_matrix, emission_matrix,
					 symbol_dictionary):

		# The entry in the transition_matrix in row x gives the transition probability to get to row x+1.  In the viterbi_matrix,
		# the entry at row x is the probability to have gotten to row x.  Thus, the probability to get to row x+1 is
		# transition_matrix[row x]*viterbi_matrix[row x]

		emit = emission_matrix[symbol_dictionary[sequence[row-1]]][column]
		big_negative_number = -1000
		key = ['MM', 'IM', 'DM']
		M_I_D_parameters = [big_negative_number, big_negative_number, big_negative_number]
		for parameter_index in range(3):
			if viterbi_matrix[row - 1][column - 1][parameter_index] != '*':
				M_I_D_parameters[parameter_index] = emit + viterbi_matrix[row-1][column-1][parameter_index] + transition_matrix[self.dictionary_of_transition_paths[key[parameter_index]]][column - 1]
		result = (M_I_D_parameters[0], M_I_D_parameters[1], M_I_D_parameters[2])
		direction_of_movement_matrix[row][column][0] = ['D','V','H',][result.index(max(result))]

		return max(result)

	def get_viterbi_path(self, viterbi_matrix, direction_matrix):
		number_of_columns = len(viterbi_matrix[0])-1
		number_of_rows = len(viterbi_matrix)-1
		column = number_of_columns
		row = number_of_rows
		list_of_directions = []

		while (row != 0) and (column != 0):


			value = viterbi_matrix[row][column]
			max_value = max(value)
			max_index = viterbi_matrix[row][column].index(max_value)
			if row == number_of_rows and column == number_of_columns:
				list_of_directions.append(['M','I','D'][max_index]+str(column))
			direction = direction_matrix[row][column][max_index]
			#print('row=',row,'column=',column,'max_index=',max_index,'direction=',direction )
			previous_direction = 'D'
			current_column = 0
			if direction == 'D':
				list_of_directions.append('M'+str(column-1))
				row -= 1
				column -= 1
			elif direction == 'V':
				list_of_directions.append('I'+str(column-1))
				row -= 1
			elif direction == 'H':
				list_of_directions.append('D'+str(column-1))
				column -= 1
		reversed_list_of_directions = list_of_directions[::-1]
		return(reversed_list_of_directions)

