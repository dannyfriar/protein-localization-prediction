# coding: utf-8
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from sklearn.model_selection import train_test_split

amino_acids = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
non_amino_letters = ['B', 'J', 'O', 'U', 'X', 'Z']


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
        
        for amino in amino_acids:
            key_string = amino + '_count'
            fasta_data[key_string].append(sequence.count(amino))
            key_string = amino + '_first50_count'
            fasta_data[key_string].append(sequence[:51].count(amino))
            key_string = amino + '_last50_count'
            fasta_data[key_string].append(sequence[:51].count(amino))
        
    return(pd.DataFrame.from_dict(fasta_data))


def main():


    print("#-------- Loading and parsing sequence files")
    # Loading train files
    train_file_names = ['cyto', 'mito', 'nuclear', 'secreted']
    train_dfs = []

    for file_name in train_file_names:
        print('Loading file %s...' % file_name)
        fasta_df = fasta_to_pandas('data/'+file_name+'.fasta')
        fasta_df.insert(loc=2, column='class', value=file_name)
        train_dfs.append(fasta_df)    
    train_full = pd.concat(train_dfs)

    # Loading blind test file
    print('Loading file blind test...')
    blind_test = fasta_to_pandas('data/blind_test.fasta')
    print("#-------- Done.")

    # Counts per class
    for file_name in train_file_names:
        class_count = sum(train_full['class'] == file_name)
        print('Count for class %s: %d' % (file_name, class_count))

    # Starting sequence frequency
    unique_seq = set([seq[-3:] for seq in train_full.sequence])

    # Split into train and test sets and export as CSV files
    y_train_full = train_full['class']
    X_train_full = train_full.drop('class', axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X_train_full, y_train_full, test_size=0.25, random_state=0)

    train = X_train
    train.insert(loc=2, column='class', value=y_train)
    test = X_test
    test.insert(loc=2, column='class', value=y_test)

    train.to_csv('data/train.csv', index=False)
    test.to_csv('data/test.csv', index=False)
    blind_test.to_csv('data/blind_test.csv', index=False)


if __name__ == "__main__":
    main()