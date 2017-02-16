# coding: utf-8
from __future__ import division
import numpy as np
import pandas as pd
import seaborn as sns
import random
import matplotlib.pyplot as plt

from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from sklearn.model_selection import train_test_split

amino_acids = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
non_amino_letters = ['B', 'J', 'O', 'U', 'X', 'Z']

def count_hydrophobic(sequence):
	"""Returns count of hydrophobic amino acids"""
	count = sequence.count('F') + sequence.count('Y') + sequence.count('W')
	count += sequence.count('H') + sequence.count('K') + sequence.count('T')
	count += sequence.count('C') + sequence.count('G') + sequence.count('A')
	count += sequence.count('V') + sequence.count('I') + sequence.count('L')
	count += sequence.count('M')
	return(count)


def fasta_to_pandas(file):
	"""Takes a fasta sequence file, creates some features and converts to pandas dataframe"""
	fasta_sequences = SeqIO.parse(open(file),'fasta')
	
	# Create ordered dict (include new features)
	fasta_data = OrderedDict()
	fasta_data['name'] = []
	fasta_data['sequence'] = []
	fasta_data['sequence_length'] = []
	fasta_data['molecular_weight'] = []
	fasta_data['isolectric_point'] = []
	fasta_data['aromaticity'] = []
	fasta_data['hydrophobicity'] = []
	fasta_data['second_helix'] = []
	fasta_data['second_turn'] = []
	fasta_data['second_sheet'] = []
	fasta_data['pct_pos_charged'] = []
	fasta_data['pct_neg_charged'] = []
	fasta_data['pct_hydrophobic'] = []

	
	for amino in amino_acids:
		key_string = amino + '_count'
		fasta_data[key_string] = []
		key_string = amino + '_first50_count'
		fasta_data[key_string] = []
		key_string = amino + '_last50_count'
		fasta_data[key_string] = []
		
	# Parse fasta files and compute features
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		fasta_data['name'].append(name)
		fasta_data['sequence'].append(sequence)
		fasta_data['sequence_length'].append(len(sequence))
		
		# Remove letters not corresponding to amino acids
		for letter in non_amino_letters:
			sequence = sequence.replace(letter, '')
		
		# Compute other features
		proto_param = ProtParam.ProteinAnalysis(sequence)
		fasta_data['molecular_weight'].append(proto_param.molecular_weight())
		fasta_data['isolectric_point'].append(proto_param.isoelectric_point())
		fasta_data['aromaticity'].append(proto_param.aromaticity())
		fasta_data['hydrophobicity'].append(proto_param.gravy())

		secondary_struct_list = proto_param.secondary_structure_fraction()
		fasta_data['second_helix'].append(secondary_struct_list[0])
		fasta_data['second_turn'].append(secondary_struct_list[1])
		fasta_data['second_sheet'].append(secondary_struct_list[2])

		fasta_data['pct_pos_charged'].append((sequence.count('H')+sequence.count('K')+sequence.count('R'))/len(sequence))
		fasta_data['pct_neg_charged'].append((sequence.count('D')+sequence.count('E'))/len(sequence))
		fasta_data['pct_hydrophobic'].append(count_hydrophobic(sequence)/len(sequence))

		for amino in amino_acids:
			key_string = amino + '_count'
			fasta_data[key_string].append(sequence.count(amino)/len(sequence))
			key_string = amino + '_first50_count'
			fasta_data[key_string].append(sequence[:51].count(amino)/50)
			key_string = amino + '_last50_count'
			fasta_data[key_string].append(sequence[:51].count(amino)/50)
		
	return(pd.DataFrame.from_dict(fasta_data))


def main():


	print("#-------- Loading and parsing sequence files")
	# Loading train files
	train_file_names = ['cyto', 'mito', 'nuclear', 'secreted']
	train_dfs = []

	for file_name in train_file_names:
		print('Loading file %s...' % file_name)
		fasta_df = fasta_to_pandas('../data/'+file_name+'.fasta')
		fasta_df.insert(loc=2, column='class', value=file_name)
		train_dfs.append(fasta_df)    
	train_full = pd.concat(train_dfs)

	# Loading blind test file
	print('Loading file blind test...')
	blind_test = fasta_to_pandas('../data/blind_test.fasta')
	print("#-------- Done.")

	# Counts per class
	for file_name in train_file_names:
		class_count = sum(train_full['class'] == file_name)
		print('Count for class %s: %d' % (file_name, class_count))

	# # Starting sequence frequency
	# unique_seq = set([seq[-3:] for seq in train_full.sequence])

	# Split into train and test sets and export as CSV files
	y_train_full = train_full['class']
	X_train_full = train_full.drop('class', axis=1)
	print("Full training data size: %d" % y_train_full.shape[0])
	X_train, X_test, y_train, y_test = train_test_split(X_train_full, y_train_full, test_size=0.25, random_state=0)

	train = X_train
	print("Training data size: %d" % y_train.shape[0])
	train.insert(loc=2, column='class', value=y_train)
	test = X_test
	print("Test data size: %d" % y_test.shape[0])
	test.insert(loc=2, column='class', value=y_test)

	train.to_csv('../data/train.csv', index=False)
	test.to_csv('../data/test.csv', index=False)
	print("Blind test data size: %d" % blind_test.shape[0])
	blind_test.to_csv('../data/blind_test.csv', index=False)


if __name__ == "__main__":
	main()