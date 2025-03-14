# DNA Sequence Analysis Toolkit

## Overview
The **DNA Sequence Analysis Toolkit** is a Python library/module crafted for processing and analyzing DNA sequences stored in FASTA files. This toolkit offers a robust set of tools for bioinformaticians, researchers, and developers working with genetic data, enabling them to perform tasks such as parsing sequences, calculating statistics, identifying open reading frames (ORFs), and detecting sequence repeats efficiently.

## Key Features
- **FASTA File Parsing**: Reads and organizes sequences from FASTA files into a dictionary for easy access.
- **Sequence Counting**: Determines the total number of sequences in a file or dataset.
- **Length Analysis**: Computes the length of each DNA sequence.
- **Extremes Identification**: Identifies the longest and shortest sequences in a dataset.
- **Forward Reading Frames**: Generates three forward reading frames for each DNA sequence.
- **Open Reading Frame <ins>(ORF)</ins> Detection**: Identifies ORFs within forward reading frames.
- **Longest ORF Search**: Locates the longest ORF across a sequence or dataset.
- **Repeat Detection**: Analyzes and counts repeated subsequences of specified lengths.

## Functions
### File Handling
- **`read_fasta(file_path)`**  
  - **Description**: Parses a FASTA file and returns a dictionary where keys are sequence IDs and values are the corresponding DNA sequences.
  - **Parameters**: 
    - `file_path` (str): Path to the FASTA file.
  - **Returns**: Dictionary of sequence IDs and sequences.

### Basic Statistics
- **`num_seq(seq)`**  
  - **Description**: Counts the total number of sequences in a FASTA file or dictionary.
  - **Parameters**: 
    - `seq` (dict or str): Dictionary of sequences or path to a FASTA file.
  - **Returns**: Integer representing the number of sequences.

- **`lengths(genes)`**  
  - **Description**: Calculates the length of each DNA sequence in the dataset.
  - **Parameters**: 
    - `genes` (dict): Dictionary of sequence IDs and sequences.
  - **Returns**: Dictionary mapping sequence IDs to their lengths.

- **`max_seq(genes)`**  
  - **Description**: Identifies the sequence with the maximum length.
  - **Parameters**: 
    - `genes` (dict): Dictionary of sequences.
  - **Returns**: Tuple of the sequence ID and its sequence.

- **`min_seq(genes)`**  
  - **Description**: Identifies the sequence with the minimum length.
  - **Parameters**: 
    - `genes` (dict): Dictionary of sequences.
  - **Returns**: Tuple of the sequence ID and its sequence.

### Advanced Analysis
- **`get_forward_reading_frames(genes)`**  
  - **Description**: Generates the three forward reading frames for each DNA sequence.
  - **Parameters**: 
    - `genes` (dict): Dictionary of sequences.
  - **Returns**: Dictionary with sequence IDs mapped to a list of three reading frames.

- **`get_forward_orf(rf)`**  
  - **Description**: Identifies open reading frames (ORFs) in each reading frame.
  - **Parameters**: 
    - `rf` (list): List of reading frames for a sequence.
  - **Returns**: List of ORFs found in the reading frames.

- **`max_len_ORF(reading_frames, num_rf=0, seq_id=None, starting_position=False)`**  
  - **Description**: Finds the longest ORF in a specified sequence or across all sequences in the dataset.
  - **Parameters**: 
    - `reading_frames` (dict): Dictionary of sequences.
    - `num_rf` (int, optional): Reading frame to analyze (1, 2, or 3). Defaults to 0 (any frame).
    - `seq_id` (str, optional): Specific sequence ID to analyze (default is None, analyzes all).
    - `starting_position` (bool, optional): If True, returns the starting position of the longest ORF (default is False).
  - **Returns**: Tuple containing the longest ORF sequence and its length (and starting position if requested).

- **`count_repeats(genes, length=1, repeat=None)`**  
  - **Description**: Counts occurrences of repeated subsequences of a specified length or a specific sequence.
  - **Parameters**: 
    - `genes` (dict): Dictionary of sequences.
    - `length` (int, optional): Length of repeats to count (default is 1).
    - `repeat` (str, optional): Specific sequence to count (default is None, counts all repeats of given length).
  - **Returns**: Dictionary of repeats and their counts (only includes repeats occurring more than once).

## Usage Notes
>[!NOTE]  
>- This toolkit is optimized for DNA sequence analysis and assumes input sequences are valid DNA (A, T, C, G).  
>- Functions are designed with flexible parameters to support diverse analysis workflows.  
>- The `count_repeats` function filters out repeats that occur only once, focusing on statistically significant patterns.

## Installation
To use the toolkit, ensure you have Python 3.x installed. Clone the repository or download the module, then import it into your project:
```python
import dna_sequence_analysis_toolkit as dna_toolkit
```

## Example
```python
# Example usage
file_path = "example.fasta"
sequences = dna_toolkit.read_fasta(file_path)
print(f"Number of sequences: {dna_toolkit.num_seq(sequences)}")
print(f"Longest sequence: {dna_toolkit.max_seq(sequences)}")
orfs = dna_toolkit.get_forward_orf(dna_toolkit.get_forward_reading_frames(sequences))
print(f"Longest ORF: {dna_toolkit.max_len_ORF(sequences)}")
```

## Author
**Eslam Hassan**  
- Contact: [heslam607@gmail.com](e-mail) 