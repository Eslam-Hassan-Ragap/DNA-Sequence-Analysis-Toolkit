# DNA-Sequence-Analysis-Toolkit

## Overview
DNA-Sequence-Analysis-Toolkit is a Python library/module designed to process and analyze DNA sequences stored in FASTA files. It provides a suite of functions to read FASTA files, count sequences, calculate sequence lengths, identify open reading frames (ORFs), and detect sequence repeats. This tool is intended for bioinformaticians, researchers, and developers working with genetic data.

## Features
* ### ***Reading FASTA files***: Parses a FASTA file and stores sequences in a dictionary.<br/>
+ ***Sequence Counting***: Counts the number of sequences in a FASTA file or dictionary.<br/>
- ***Sequence Lengths***: Calculates the length of each sequence.<br/>
- ***Longest and Shortest Sequences***: Identifies sequences with the maximum and minimum lengths.<br/>
+ ***Forward Reading Frames***: Generates three forward reading frames for each sequence.<br/>
* ***Open Reading Frames <ins>(ORFs)</ind>***: Identifies ORFs in forward reading frames.<br/>
* ***Longest ORF***: Finds the longest ORF in a given sequence or dataset.<br/>
- ***Repeat Analysis***: Counts repeated sequences of specified lengths.<br/>

## Functions
### ***File Handling***<br/>
read_fasta(file_path): Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values

### ***Basic Statistics***<br/>
`num_seq(seq)`: Counts the number of sequences in a FASTA file or dictionary
`lengths(genes)`: Calculates the length of each sequence
`max_seq(genes)`: Finds the sequence with the maximum length
`min_seq(genes)`: Finds the sequence with the minimum length

### ***Advanced Analysis***<br/>
`get_forward_reading_frames(genes)`: Generates forward reading frames for DNA sequences
`get_forward_orf(rf)`: Finds open reading frames (ORFs) in each reading frame
`max_len_ORF(genes, rf=0, seq_id=None, starting_position=False)`: Locates the longest ORF in gene sequences
`count_repeates(genes, length=1, repeat=None)`: Counts repeats of specific lengths or sequences
