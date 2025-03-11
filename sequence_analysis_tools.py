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
        file (dict or str): FASTA file or dictionary of sequences.

    Returns:
        int: Number of sequences in the file or dictionary.
    """
    if type(seq)==dict: # Check if seq is a dictionary
        return len(seq)
    elif type(seq)==str: # Check if seq is a string
        return len(read_fasta(seq))