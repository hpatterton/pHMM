import numpy as np
import copy
from collections import Counter
from HMM_Sequence import HMM_Sequence

class HMM_transitions:

    def __init__(self, filepath, emission_pseudocounts = True, transition_pseudocounts = True, molecule = 'D', normalise_emissions_to_background = True, sequence_conservation = False):  # should we calculate pseudocounts or raw ratios?
        self.sequence_array = []
        self.read_multiple_alignment(filepath)
        #   Should a block be IDed as conserved in terms of aa composition, or % '-'s?
        self.number_of_sequences = len(self.sequence_array)
        self.number_of_columns = len(self.sequence_array[0])
        self.delete_list = ['']*self.number_of_sequences
        self.dictionary_of_transition_paths = {'MM': 0, 'MI': 1, 'MD': 2, 'IM': 3, 'II': 4, 'ID': 5, 'DM': 6, 'DI': 7, 'DD': 8}
        self.number_of_paths = len(self.dictionary_of_transition_paths)
        self.transition_paths = [] # remember to initialize to 0
        self.conservation_cutoff = 0.5
        self.column_categories = []
        self.molecule_type = molecule
        self.emission_pseudocounts = emission_pseudocounts
        self.transition_pseudocounts = transition_pseudocounts
        self.normalise_emissions_to_background = normalise_emissions_to_background
        if self.molecule_type == 'D':
            self.background_emission_frequency = 0.25
        else:
            self.background_emission_frequency = 0.05
        self.sequence_conservation = sequence_conservation
        self.number_of_matched_columns = self.make_column_category_list()
        self.transition_probabilities = [[0]*len(self.dictionary_of_transition_paths) for i in range(self.number_of_matched_columns+1)]
        self.dictionary_of_amino_acids = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
                                                  'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
                                                  'V': 17, 'W': 18, 'Y': 19}
        self.dictionary_of_nucleotides = {'G': 0, 'A': 1, 'T': 2, 'C': 3}

        return

    def get_sequence(self):
        mystring = []
        for i in range(len(self.sequence_array)):
            mytemp = ''
            for j in range(len(self.sequence_array[0])):
                mytemp = mytemp + self.sequence_array[i][j]
            mystring.append(mytemp)
        return mystring

    def read_multiple_alignment(self, filepath):
        read_sequence = HMM_Sequence()
        if read_sequence.is_fasta_format(filepath):
            sequences = read_sequence.read_fasta_file(filepath)
        else:
            sequences = read_sequence.read_text_file(filepath)
        if sequences:
            number_of_sequences = len(sequences)
            self.sequence_array = np.array([list(sequences[i]) for i in range(number_of_sequences)])
            #print(self.sequence_array)
            return True
        else:
            return False

    def is_conserved(self, row):
        conserved = False
        row_as_list = row.tolist()
        if self.sequence_conservation == False:
            number_of_hyphens = row_as_list.count('-')
            if (self.number_of_sequences-number_of_hyphens)/self.number_of_sequences >= self.conservation_cutoff:
                conserved = True
        return conserved

    def make_column_category_list(self):
        number_of_matched_columns = 0
        for i in range(self.number_of_columns):
            if self.is_conserved(self.sequence_array[:,i]):
                self.column_categories.append(0) # 0 is for a conserved columns
                number_of_matched_columns += 1
            else:
                self.column_categories.append(1) # 1 is for an insert column
        return number_of_matched_columns

    def interpret_insert_block(self, array):
        number_of_pass_throughs = 0
        number_of_inserts = 0
        number_of_return_loops = 0
        for i in range(self.number_of_sequences): # add all the rows composed entirely of '-'
            length_of_row = len(array[i])
            if array[i].count('-') == length_of_row:
                number_of_pass_throughs += 1
            if array[i].count('-') < length_of_row:
                number_of_inserts += 1
                number_of_return_loops += (length_of_row - array[i].count('-')- 1)
        return (number_of_pass_throughs, number_of_inserts, number_of_return_loops)

    # make iD list where any row with insertion becomes the number of insertions, all else 'P' for pass-through
    def collapse_insert_array_to_list(self, insert_array):
        insert_list = []
        number_of_columns = len(insert_array[0])
        number_of_rows = len(insert_array)
        for i in range(number_of_rows):
            if insert_array[i].count('-') < number_of_columns:
                insert_list.append(str(number_of_columns-insert_array[i].count('-')))
            else:
                insert_list.append('P')
        return insert_list

    def get_next_block(self, index): # the index here is to the column categories list

        if index >= (self.number_of_columns): # we have a start and stop to add
            return ([],-1, False)
        if self.column_categories[index] == 0:  # conserved, return only column
            return (self.sequence_array[:,index].tolist(), index+1, True)
        else:
            # how many non_conserved columns do we have
            end_of_block = 1
            while (index + end_of_block < self.number_of_columns) and (self.column_categories[index + end_of_block] == 1):
                end_of_block += 1
            return (self.sequence_array[:,index:index+end_of_block].tolist(), index + end_of_block, False)

    def get_transition_probabilities(self):
        profile_index = 0
        current_index = 0

        #START
        current_block, next_index, next_block_conserved = self.get_next_block(current_index)
        current_index = next_index
        if next_block_conserved:
            self.calculate_start_M(profile_index, current_block)
        else:
            next_block, next_index, next_block_conserved = self.get_next_block(current_index)
            self.calculate_start_EM(profile_index, current_block, next_block)
            current_block = copy.deepcopy(next_block) # to ensure a deep copy

        current_block_conserved = next_block_conserved
        next_block, next_index, next_block_conserved = self.get_next_block(next_index)
        profile_index += 1

        while next_block != []: # test if next block is END

            if current_block_conserved:
                if next_block_conserved:
                    self.calculate_body_MM(profile_index, current_block, next_block)
                else:
                    self.calculate_body_ME(profile_index, next_index, current_block, next_block)
            else: # E
                next_block, next_index, next_block_conserved = self.get_next_block(next_index)
                if next_block == []:
                    self.calculate_E_end(profile_index, current_block)
                    return
                else:
                    current_block = copy.deepcopy(next_block)
                    current_block_conserved = next_block_conserved
                    next_block, next_index, next_block_conserved = self.get_next_block(next_index)
                    if next_block == []:
                        self.calculate_M_end(profile_index, current_block)
                        return
                    else:
                        if next_block_conserved:
                            self.calculate_body_MM(profile_index, current_block, next_block)
                        else:
                            self.calculate_body_ME(profile_index, next_index, current_block, next_block)

            profile_index += 1

            # If the current_block-next_block pair is conserved-non_conserved, you have already acalculated the contribution
            # of the non_conserved block.  The current block must now be the one after non_conserved.  You know that this will
            # be a conserved block, because that always follows a non_conserved, excdpt when non_conserved was the last block.
            # In that case the next_block call will return [], and the function will terminate

            if next_block_conserved == True:
                current_block = copy.deepcopy(next_block)
                current_block_conserved = next_block_conserved
                next_block, next_index, next_block_conserved = self.get_next_block(next_index)
            else:
                current_block, next_index, next_block_conserved = self.get_next_block(next_index)
                current_block_conserved = next_block_conserved
                next_block, next_index, next_block_conserved = self.get_next_block(next_index)

        if current_block_conserved:
            self.calculate_M_end(profile_index, current_block)
        else:
            profile_index -= 1
            self.calculate_E_end(profile_index, current_block)
        #print(self.transition_probabilities)
        return self.transition_probabilities

    def get_zeroed_transition_dictionary(self):
        transition_dictionary = {'MD':0,'MI':0,'MM':0,'ID':0,'II':0,'IM':0,'DD':0,'DI':0,'DM':0}
        return transition_dictionary

    def calculate_transition_probabilities(self, profile_index, transition_numbers, transition_pseudocounts=True):
        #print(profile_index,transition_numbers)
        if transition_pseudocounts == False:
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['MM']] = (self.number_of_sequences - transition_numbers['MI'] - transition_numbers['MD']) / self.number_of_sequences
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['MD']] = transition_numbers['MD'] / self.number_of_sequences
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['MI']] = transition_numbers['MI'] / self.number_of_sequences
            if (transition_numbers['II'] + transition_numbers['IM']) > 0:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['ID']] = transition_numbers['ID'] / (transition_numbers['II'] + transition_numbers['IM'])
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['II']] = transition_numbers['II'] / (transition_numbers['II'] + transition_numbers['IM'])
            else:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['ID']] = 0
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['II']] = 0
            if (transition_numbers['II'] + transition_numbers['IM']) > 0:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['IM']] = (transition_numbers['IM'] - transition_numbers['ID']) / (transition_numbers['II'] + transition_numbers['IM'])
            else:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['IM']] = 0
            if (self.delete_list.count('D') + transition_numbers['DD'] + transition_numbers['DM'] + transition_numbers['DI']) > 0:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DD']] = transition_numbers['DD'] / (self.delete_list.count('D') + transition_numbers['DM'] + transition_numbers['DI'])
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DM']] = transition_numbers['DM'] / (self.delete_list.count('D') + transition_numbers['DM'] + transition_numbers['DI'])
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DI']] = transition_numbers['DI'] / (self.delete_list.count('D') + transition_numbers['DM'] + transition_numbers['DI'])
            else:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DD']] = 0
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DM']] = 0
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DI']] = 0
        else:
        #
        # For transition pseudo counts, add one to the transition number from state under consideration to the next state
        # for example count of 4 transitions from M0 to M1 becomes 4+1=5
        # Divide by the sum of all the transitions from this state to the other 3 possible states
        # For example, if there are 4 M1 to M2, 2 M1 to I1 and 1 M1 to D2 transitions, the denominator becomes
        # (4+1)+(2+1)+(1+1)
        # So, with transition pseudocounts, the fractions becomes (4+1)/[(4+1)+(2+1)+(1+1)]
        #
            MM = self.number_of_sequences - transition_numbers['MI'] - transition_numbers['MD']
            From_M_transition_denominator = (MM+1)+(transition_numbers['MI']+1)+(transition_numbers['MD']+1)
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['MM']] = (MM+1) / From_M_transition_denominator
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['MD']] = (transition_numbers['MD']+1) / From_M_transition_denominator
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['MI']] = (transition_numbers['MI']+1) / From_M_transition_denominator

            From_I_transition_denominator = (transition_numbers['IM']+1)+(transition_numbers['II']+1)+(transition_numbers['ID']+1)
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['ID']] = (transition_numbers['ID']+1) / From_I_transition_denominator
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['II']] = (transition_numbers['II']+1) / From_I_transition_denominator
            self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['IM']] = (transition_numbers['IM'] - transition_numbers['ID']+1) / From_I_transition_denominator

            # There is no D0 state, so simply set the values to 0
            if profile_index == 0:
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DD']] = 0
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DM']] = 0
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DI']] = 0
            else:
                From_D_transition_denominator = (transition_numbers['DM']+1)+(transition_numbers['DI']+1)+(transition_numbers['DD']+1)
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DD']] = (transition_numbers['DD']+1) / From_D_transition_denominator
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DM']] = (transition_numbers['DM']+1) / From_D_transition_denominator
                self.transition_probabilities[profile_index][self.dictionary_of_transition_paths['DI']] = (transition_numbers['DI']+1) / From_D_transition_denominator
        return self.transition_probabilities[profile_index]

    def calculate_start_M(self, profile_index, current_block):

        transitions = self.get_zeroed_transition_dictionary()
        number_of_deletions = 0
        for i in range(len(self.delete_list)):
            if current_block[i] == '-':
                self.delete_list[i] = 'D'
                number_of_deletions += 1
        transitions['MD'] = number_of_deletions
        transitions['MM'] = self.number_of_sequences
        self.calculate_transition_probabilities(profile_index, transitions)
        return

    def calculate_start_EM(self, profile_index, current_block, next_block):
        transitions = self.get_zeroed_transition_dictionary()
        number_of_pass_throughs, number_of_inserts, number_of_return_loops = self.interpret_insert_block(current_block)
        insert_list = self.collapse_insert_array_to_list(current_block)
        for i in range(self.number_of_sequences):
            if insert_list[i] == 'P' and next_block[i] == '-':
                transitions['MD'] += 1
                self.delete_list[i] = 'D'
            if insert_list[i] != 'P' and next_block[i] == '-':
                transitions['ID'] += 1
                self.delete_list[i] = 'D'
        transitions['II'] = number_of_return_loops
        transitions['IM'] = number_of_inserts
        transitions['MI'] = number_of_inserts
        transitions['MM'] = self.number_of_sequences
        self.calculate_transition_probabilities(profile_index, transitions)
        return

    def calculate_body_MM(self, profile_index, current_block, next_block):

        transitions = self.get_zeroed_transition_dictionary()
        for i in range(len(self.delete_list)):
            if current_block[i] == '-' and next_block[i] != '-':
                transitions['DM'] += 1
                self.delete_list[i] = ''
            if current_block[i] != '-' and next_block[i] == '-':
                transitions['MD'] += 1
                self.delete_list[i] = 'D'
            if current_block[i] == '-' and next_block[i] == '-':
                transitions['DD'] += 1
            transitions['MM'] = self.number_of_sequences

        self.calculate_transition_probabilities(profile_index, transitions)
        return

    def calculate_body_ME(self, profile_index, next_index, current_block, next_block):
        transitions = self.get_zeroed_transition_dictionary()
        number_of_pass_throughs, number_of_inserts, number_of_return_loops = self.interpret_insert_block(next_block)
        insert_list = self.collapse_insert_array_to_list(next_block)
        next_next_block, next_next_index, next_conserved = self.get_next_block(next_index)
        if next_next_block == []:
            return
        else:
            for i in range(len(self.delete_list)):
                if current_block[i] == '-' and insert_list[i] == 'P' and next_next_block[i] == '-':
                    transitions['DD'] += 1
                if current_block[i] == '-' and insert_list[i] == 'P' and next_next_block[i] != '-':
                    transitions['DM'] += 1
                    self.delete_list[i] = ''
                if insert_list[i] != 'P' and next_next_block[i] == '-':
                    transitions['ID'] += 1
                    self.delete_list[i] = 'D'
                if self.delete_list[i] == 'D' and insert_list[i] != 'P':
                    transitions['DI'] += 1
                    self.delete_list[i] = ''
                if current_block[i] != '-' and insert_list[i] == 'P' and next_next_block[i] == '-':
                    transitions['MD'] += 1
                    self.delete_list[i] = 'D'
            transitions['MI'] = number_of_inserts
            transitions['II'] = number_of_return_loops
            transitions['IM'] = number_of_inserts
            transitions['MM'] = self.number_of_sequences

            self.calculate_transition_probabilities(profile_index, transitions)

        return

    def calculate_E_end(self, profile_index, current_block):

        transitions = self.get_zeroed_transition_dictionary()
        number_of_pass_throughs, number_of_inserts, number_of_return_loops = self.interpret_insert_block(current_block)
        insert_list = self.collapse_insert_array_to_list(current_block)

        for i in range(len(self.delete_list)):
            if self.delete_list[i] == 'D' and insert_list[i] != 'P':
                transitions['DI'] += 1
                self.delete_list[i] = ''
            if self.delete_list[i] == 'D' and insert_list[i] == 'P':
                transitions['DM'] += 1
                self.delete_list[i] = ''
        transitions['MI'] = number_of_inserts
        transitions['II'] = number_of_return_loops
        transitions['IM'] = number_of_inserts
        transitions['MM'] = self.number_of_sequences

        self.calculate_transition_probabilities(profile_index, transitions)
        return

    def calculate_M_end(self, profile_index, current_block):
        transitions = self.get_zeroed_transition_dictionary()
        transitions['DM'] = self.delete_list.count('D')
        for i in range(len(self.delete_list)):
            self.delete_list[i] = ''
        transitions['MM'] = self.number_of_sequences
        self.calculate_transition_probabilities(profile_index, transitions)
        return

    def get_pseudocounts(self, numerator, denominator, molecule_type):
        if molecule_type == 'P':
            return((numerator+1)/(denominator+20))
        else:
            return((numerator+1)/(denominator+4))

# list of symbols can be a 1D or 2D list array
    def get_emmission_list(self,list_of_symbols, molecule_type):
        if molecule_type == 'P':
            dictionary = self.dictionary_of_amino_acids
        elif molecule_type == 'D':
            dictionary = self.dictionary_of_nucleotides
        else:
            return []
        emmision_array = [0 for i in range(len(dictionary))]
        #print('emission array',emmision_array)
        symbol_count = Counter(self.get_flattened_filtered_list(list_of_symbols))
        total_symbols = sum(symbol_count.values())
        denominator = total_symbols
        for character in dictionary:
            numerator = symbol_count[character]
            if self.emission_pseudocounts:
                if self.normalise_emissions_to_background:
                    emmision_array[dictionary[character]] = self.get_pseudocounts(numerator, denominator, molecule_type)/self.background_emission_frequency
                else:
                    emmision_array[dictionary[character]] = self.get_pseudocounts(numerator, denominator,
                                                                                  molecule_type)
            else:
                if self.normalise_emissions_to_background:
                    emmision_array[dictionary[character]] = (numerator/denominator)/self.background_emission_frequency
                else:
                    emmision_array[dictionary[character]] = (numerator / denominator)

        return emmision_array

    def get_flattened_filtered_list(self, sequence_list): # not a numpy array
        flattened_list = []
        for i in sequence_list:
            flattened_list += i
        flattened_filtered_list = [symbol for symbol in flattened_list if symbol != '-']
        return flattened_filtered_list

    def guess_molecule_type(self, flattened_filtered_list):
        return self.molecule_type
        # test_set = set(self.dictionary_of_nucleotides)
        # length = len(flattened_filtered_list)
        # index = 0
        # type = '?'
        # while index < length and flattened_filtered_list[index] in test_set:
        #     type = 'D'
        #     index += 1
        # if index < length:
        #     type = '?'
        # if type == 'D':
        #     return (type)
        # dictionary_of_amino_acids = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
        #                                   'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
        #                                   'V': 17, 'W': 18, 'Y': 19}
        # test_set = set(self.dictionary_of_amino_acids)
        # length = len(flattened_filtered_list)
        # index = 0
        # while index < length and flattened_filtered_list[index] in test_set:
        #     type = 'P'
        #     index += 1
        # if index < length:
        #     type = '?'
        # return type

    def get_emmision_blocks(self):
        index = 0
        profile_index = 0
        number_of_emmision_blocks = self.number_of_matched_columns+1
        if self.column_categories[0] == 1:
            number_of_emmision_blocks += 1
        emmision_blocks = [[[],[]] for i in range(number_of_emmision_blocks)]
        next_block, next_index, is_conserved = self.get_next_block(index)
        flattened_filtered_list = self.get_flattened_filtered_list(next_block)
        molecule_type = self.guess_molecule_type(flattened_filtered_list)
        emmision_list = self.get_emmission_list(flattened_filtered_list, molecule_type)

        # first block
        if index == 0 and is_conserved == False:
            emmision_blocks[profile_index][1] = emmision_list
        else:
            profile_index += 1
            emmision_blocks[profile_index][0] = emmision_list
            index = next_index
        # rest of the blocks
        while index < len(self.column_categories):
            next_block, next_index, is_conserved = self.get_next_block(index)
            flattened_filtered_list = self.get_flattened_filtered_list(next_block)
            molecule_type = self.guess_molecule_type(flattened_filtered_list)
            emmision_list = self.get_emmission_list(flattened_filtered_list, molecule_type)
            if is_conserved:
                profile_index += 1
                emmision_blocks[profile_index][0] = emmision_list
            else:
                emmision_blocks[profile_index][1] = emmision_list
            index = next_index
        return emmision_blocks, molecule_type




























