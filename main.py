from pHMM_Viterbi import pHMM_Viterbi
from HMM_transitions import HMM_transitions
import prettytable
import math

def transpose_list_of_lists(start_list):
    transposed_list_of_lists = []
    number_of_rows = len(start_list)
    number_of_columns = len(start_list[0])
    for row in range(number_of_columns):
        new_row = []
        for column in range(number_of_rows):
            new_row.append(start_list[column][row])
        transposed_list_of_lists.append(new_row)
    return(transposed_list_of_lists)


hmm_transitions = HMM_transitions('sequence14.txt',True, True, 'D')
query_sequence = 'CCTGTC'
print('Sequence alignment')
print(hmm_transitions.sequence_array,'\n')
emmision_blocks = hmm_transitions.get_emmision_blocks()
molecule_type = emmision_blocks[1]
transition_matrix = hmm_transitions.get_transition_probabilities()
number_of_matched_columns = hmm_transitions.number_of_matched_columns

#print('Transition count matrix:','\n')
#print(transition_matrix)
transposed_transition_probabilities = transpose_list_of_lists(transition_matrix)
number_of_rows = len(transposed_transition_probabilities)
number_of_columns = len(transposed_transition_probabilities[0])

for row in range(number_of_rows):
    for column in range(number_of_columns):
        if transposed_transition_probabilities[row][column] != 0:
            transposed_transition_probabilities[row][column] = math.log10(transposed_transition_probabilities[row][column])
        else:
            transposed_transition_probabilities[row][column] = 0
transition_matrix = transposed_transition_probabilities

transition_matrix_table = prettytable.PrettyTable()
field_row = []
for column in range(number_of_columns):
    field_row.append(str(column))
transition_matrix_table.field_names = field_row
for row in range(number_of_rows):
    transition_matrix_table.add_row(transposed_transition_probabilities[row])
print('Transition matrix:')
print(transition_matrix_table,'\n')

list_of_lists = []
for entry in range(len(emmision_blocks[0])):
    if len(emmision_blocks[0][entry][0]) == 0:
        if(emmision_blocks[1] == 'D'):
            emmision_blocks[0][entry][0] = [0]*4
        else:
            emmision_blocks[0][entry][0] = [0] * 20

    list_of_lists.append(emmision_blocks[0][entry][0])

emission_matrix = transpose_list_of_lists(list_of_lists)
for row in range(len(emission_matrix)):
    for column in range(len(emission_matrix[0])):
        if emission_matrix[row][column] != 0:
            emission_matrix[row][column] = math.log10(emission_matrix[row][column])
        else:
            emission_matrix[row][column] = 0

if molecule_type == 'D':
    symbol_dictionary = hmm_transitions.dictionary_of_nucleotides
else:
    symbol_dictionary = hmm_transitions.dictionary_of_amino_acids


viterbi = pHMM_Viterbi()
viterbi_matrix, direction_of_movement_matrix = viterbi.calculate_Viterbi(query_sequence, number_of_matched_columns, transition_matrix, emission_matrix, symbol_dictionary)

print('Viterbi matrix:')
number_of_rows = len(viterbi_matrix)
number_of_columns = len(viterbi_matrix[0])
viterbi_table = prettytable.PrettyTable()
field_row = []
for state in range(number_of_columns):
    field_row.append(str(state))
viterbi_table.field_names = field_row
for row in range(number_of_rows):
    for parameter in range(3):
        table_row = []
        for column in range(number_of_columns):
            if viterbi_matrix[row][column][parameter] != '*':
                entry = ['M: ','I: ','D: '][parameter] + str('%.3f' % viterbi_matrix[row][column][parameter])
            else:
                entry = ['M: ', 'I: ', 'D: '][parameter] + '*'
            table_row.append(entry)
        viterbi_table.add_row(table_row)
        if parameter == 2 and row != number_of_rows-1:
            table_row=[]
            for c in range(number_of_columns):
                table_row.append('--------')
            viterbi_table.add_row(table_row)
print(viterbi_table,'\n')


                #+ '\n' + 'I: ' + str(viterbi_matrix[row][column][1]) + '\n' + 'D: ' + str(viterbi_matrix[row][column][2])

        # for i in range(3):
        #     entry = ['M: ','I: ','D: '][i]
        #     if type(viterbi_matrix[row][column][0]) == float:
        #         entry += '{:.3f}'.format(viterbi_matrix[row][column][i])
        #     elif type(viterbi_matrix[row][column][i]) == str and viterbi_matrix[row][column][i] == '*':
        #         entry += '*'
        #     if i < 2:
        #         entry += '\n'
#         table_row.append(entry)
#     viterbi_table.add_row(table_row)
# print(viterbi_table)

direction_of_movement_table = prettytable.PrettyTable()
field_row = []
for column in range(number_of_columns):
    field_row.append('M'+str(column))
direction_of_movement_table.field_names = field_row
for row in range(number_of_rows):
    for parameter in range(3):
        table_row = []
        for column in range(number_of_columns):
            entry = ['M: ','I: ','D: '][parameter] + str(direction_of_movement_matrix[row][column][parameter])
            table_row.append(entry)
        direction_of_movement_table.add_row(table_row)
        if parameter == 2 and row != number_of_rows-1:
            table_row=[]
            for c in range(number_of_columns):
                table_row.append('----')
            direction_of_movement_table.add_row(table_row)
print('Direction matrix:')
print(direction_of_movement_table,'\n')


# Start in the last block and work your way to the 0,0 block

directions = viterbi.get_viterbi_path(viterbi_matrix, direction_of_movement_matrix)
print('Direction:','\n',directions)












