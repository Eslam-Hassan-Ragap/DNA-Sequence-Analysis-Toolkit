# DNA Sequence Analysis Toolkit

## Overview
The **DNA Sequence Analysis Toolkit** is a Python-based module designed to help researchers, bioinformaticians, and students perform efficient analysis of DNA sequences stored in FASTA files. It provides functions for parsing, exploring, and analyzing genetic data—covering sequence statistics, reading frames, open reading frames (ORFs), and repeat pattern discovery.

## Key Features
- **FASTA Parsing** – Import and organize DNA sequences into Python dictionaries.
- **Sequence Statistics** – Count sequences, compute lengths, and identify extremes.
- **Reading Frame Generation** – Automatically produce all three forward reading frames.
- **ORF Detection** – Identify and extract open reading frames from DNA sequences.
- **Longest ORF Search** – Find the longest ORF and optionally its genomic start position.
- **Repeat Analysis** – Detect repeated subsequences of any length or specific sequence.
- **Top Repeat Finder** – Identify the most frequent repeat patterns.
- **Max Repeat Finder** – Retrieve the subsequences with the absolute maximum repeat counts.

---

## Functions

### File Handling
#### `read_fasta(file_path)`
Reads a FASTA file and returns a dictionary where keys are sequence IDs and values are the corresponding sequences.  
**Parameters:**  
- `file_path` *(str)* – Path to the FASTA file.  
**Returns:** `dict`

---

### Basic Statistics
#### `num_seq(seq)`
Counts the total number of sequences in a FASTA file or a sequence dictionary.  
**Parameters:**  
- `seq` *(dict or str)* – Dictionary or FASTA file path.  
**Returns:** `int`

#### `lengths(genes)`
Calculates the length of each DNA sequence.  
**Returns:** `dict` mapping IDs → lengths.

#### `max_seq(genes)`
Finds the sequence(s) with the maximum length.  
**Returns:** `dict` of sequence ID(s) and length.

#### `min_seq(genes)`
Finds the sequence(s) with the minimum length.  
**Returns:** `dict` of sequence ID(s) and length.

---

### Reading Frames and ORFs
#### `get_forward_reading_frames(genes)`
Generates all three forward reading frames for each DNA sequence.  
**Returns:** `dict` mapping IDs → reading frames.

#### `get_forward_orf(rf)`
Detects open reading frames (ORFs) in each forward reading frame using standard start (`ATG`) and stop codons (`TAA`, `TGA`, `TAG`).  
**Returns:** `dict` mapping sequence IDs → ORFs found.

#### `max_len_ORF(reading_frames, num_rf=0, seq_id=None, starting_position=False)`
Finds the longest ORF across all sequences or within a specific reading frame/sequence.  
**Parameters:**  
- `num_rf` *(int, optional)* – Reading frame (1–3).  
- `seq_id` *(str, optional)* – Specific sequence ID.  
- `starting_position` *(bool, optional)* – If True, returns 1-based start position.  
**Returns:**  
- If `starting_position=False`: `dict` with ID → ORF length.  
- If `starting_position=True`: tuple of (`dict`, `start_position`).

---

### Repeat Analysis
#### `count_repeates(genes, length=1, repeat=None)`
Counts repeated subsequences in DNA sequences.  
**Parameters:**  
- `genes` *(dict)* – DNA sequences.  
- `length` *(int)* – Length of repeat to search (default `1`, meaning all).  
- `repeat` *(str)* – Specific repeat sequence (optional).  
**Returns:** `dict` of repeats and their counts per sequence.

#### `most_frequent_repeats(repeates, top_n=5)`
Finds the top N most frequent repeat sequences across all input sequences.  
**Returns:** `dict` of repeats → total count.

#### `max_repet(repeates)`
Finds repeats with the **maximum total occurrence** across all sequences.  
**Returns:** `dict` of repeat(s) with the maximum frequency.

---

## Example Usage
```python
import sequence_analysis_tools as seq_tools

file_path = "example.fasta"
genes = seq_tools.read_fasta(file_path)

print("Total sequences:", seq_tools.num_seq(genes))
print("Longest sequence:", seq_tools.max_seq(genes))
print("Shortest sequence:", seq_tools.min_seq(genes))

rf = seq_tools.get_forward_reading_frames(genes)
orfs = seq_tools.get_forward_orf(rf)
longest_orf = seq_tools.max_len_ORF(rf, starting_position=True)
print("Longest ORF and position:", longest_orf)

repeats = seq_tools.count_repeates(genes, length=6)
top_repeats = seq_tools.most_frequent_repeats(repeats)
print("Top repeated motifs:", top_repeats)
```

---

## Notes
- The toolkit assumes input sequences are valid DNA (A, T, C, G).
- All functions are independent, allowing flexible workflow integration.
- Repeat functions ignore single-occurrence motifs by default.

---

## Author
**Eslam Hassan**  
📧 [eslam.hassan.ragap@gmail.com](mailto:eslam.hassan.ragap@gmail.com)
