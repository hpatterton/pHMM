class HMM_Sequence:

    def __init__(self):
        self.dictionary_of_amino_acids = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
                                                  'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
                                                  'V': 17, 'W': 18, 'Y': 19}
        self.dictionary_of_nucleotides = {'G': 0, 'A': 1, 'T': 2, 'C': 3}
        return

    def read_fasta_file(self, filepath):
        f = open(filepath, 'r')
        sequence = []
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                if len(current_sequence) > 0:
                    sequence.append(current_sequence)
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        sequence.append(current_sequence)

        sequence_length = len(sequence[0])
        number_of_sequences = len(sequence)
        for i in range(1, len(sequence)):
            if len(sequence[i]) != sequence_length:
                return []
        return sequence

    def is_fasta_format(self, filepath):
        f = open(filepath, 'r')
        name = []
        sequence = []
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                name.append(line[1:].rstrip('\n'))
                if len(current_sequence) > 0:
                    sequence.append(current_sequence)
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        sequence.append(current_sequence)
        if len(name) == len(sequence):
            return True
        else:
            return False

    def read_text_file(self, filepath):
        f = open(filepath, 'r')
        sequence = []
        for line in f:
            sequence.append(line.rstrip('\n'))

        sequence_length = len(sequence[0])
        number_of_sequences = len(sequence)
        for i in range(1, number_of_sequences):
            if len(sequence[i]) != sequence_length:
                return []
        return sequence

