def read_fasta(file_path):
    """Function:
    Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values.

    Params:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: Dictionary with sequence IDs as keys and sequences as values.
    """
    
    with open(file_path, "r") as file_fa:
        genes={}
        id_gene=None
        for line in file_fa:
            line=line.strip()
            if(line.startswith('>')):
                id_gene=line[1:line.index(" ")]
                if id_gene not in genes:
                    genes[id_gene]=""
            elif id_gene:    
                genes[id_gene]+=line
    return genes

def num_seq(seq):
    """Function:
    Count the number of sequences in a FASTA file or a dictionary of sequences.

    Args:
        seq (dict or str): FASTA file or dictionary of sequences.

    Returns:
        int: Number of sequences in the file or dictionary.
    """
    if type(seq)==dict: # Check if seq is a dictionary
        return len(seq)
    elif type(seq)==str: # Check if seq is a string
        return len(read_fasta(seq))
    
def lengths(genes): 
    """Fuction:
    Calculate the length of each sequence in a dictionary of sequences or a FASTA file.

    Args:
        genes (dict or str): Dictionary of sequences or path to a FASTA file.

    Returns:
        dict: Dictionary with sequence IDs as keys and their lengths as values.
    """
    if type(genes) == str: # Check if genes is a FASTA file path
        genes=read_fasta(genes) # Read the FASTA file and store the sequences in a dictionary
    len_genes={}
    for id_gene,seq_gene in genes.items(): # Iterate over the sequences in the dictionary 
        len_genes[id_gene]=len(seq_gene)
    
    return len_genes

def max_seq(genes):
    """Function:
    Find the sequence with the maximum length in a dictionary of sequences.

    Args:
        genes (dict or str): Dictionary of sequences or path to a FASTA file.

    Returns:
        dict: Dictionary with the sequence ID of maximum length and its length.
    """
    if type(genes)==str:  # Check if genes is a FASTA file path
        genes=read_fasta(genes) # Read the FASTA file and store the sequences in a dictionary
        
    max_len_gene={}
    
    max_len_seq=len(max(genes.values(),key=len)) # Find the length of the longest sequence
    
    for id_gene , seq_gene in genes.items():     # Iterate over the sequences in the dictionary
        
        if len(seq_gene)==max_len_seq:          # Check if the length of the sequence is equal to the maximum length
            max_len_gene[id_gene]=max_len_seq
            
    return max_len_gene

def min_seq(genes):
    """Funciont:
    Find the sequence with the minimum length in a dictionary of sequences.

    Args:
        genes (dict or str): Dictionary of sequences or path to a FASTA file.

    Returns:
        dict: Dictionary with the sequence ID of minimum length and its length.
    """
    if type(genes)==str:  # Check if genes is a FASTA file path
        genes=read_fasta(genes) # Read the FASTA file and store the sequences in a dictionary
        
    min_len_gene={}
    
    min_len_seq=len(min(genes.values(),key=len)) # Find the length of the shortest sequence
    
    for id_gene , seq_gene in genes.items(): # Iterate over the sequences in the dictionary
        if len(seq_gene)==min_len_seq:  # Check if the length of the sequence is equal to the minimum length
            min_len_gene[id_gene]=min_len_seq
            
    return min_len_gene