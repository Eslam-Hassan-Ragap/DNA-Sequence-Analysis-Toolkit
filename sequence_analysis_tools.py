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
    Calculate the length of each sequence in a dictionary of sequences.

    Args:
        genes (dict): Dictionary of sequences.

    Returns:
        dict: Dictionary with sequence IDs as keys and their lengths as values.
    """
    len_genes={}
    
    for id_gene,seq_gene in genes.items(): # Iterate over the sequences in the dictionary 
        len_genes[id_gene]=len(seq_gene)
    
    return len_genes

def max_seq(genes):
    """Function:
    Find the sequence with the maximum length in a dictionary of sequences.

    Args:
        genes (dict): Dictionary of sequences.

    Returns:
        dict: Dictionary with the sequence ID of maximum length and its length.
    """
        
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
        genes (dict): Dictionary of sequences.

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

def get_forward_reading_frames(genes):
    """Function:
    Generate forward reading frames for each sequence in a dictionary of sequences.

    Args:
        genes (dict): Dictionary where keys are sequence IDs and values are DNA sequences.

    Returns:
        dict: Dictionary where each sequence ID maps to its forward reading frames.
    """
    f_readings={} 
    
    for id_gene,seq_gene in genes.items(): # Iterate over the sequences in the dictionary
        subreading={} 
        for i in range(3): # Iterate over a sequence to generate forward reading frames
            num_frame=f"f_reading_frame_{i+1}"
            subreading[num_frame]=seq_gene[i:] # Generate forward reading frames
            
        f_readings[id_gene]=subreading # Store the forward reading frames for the sequence
        
    return f_readings

def get_forward_orf(rf):
    """Function:
    Find the open reading frames (ORFs) in each forward reading frame of a sequence.

    Args:
        rf (dict): Dictionary where each sequence ID maps to its forward reading frames.

    Returns:
        dict: Dictionary where each sequence ID maps to its ORFs in the forward reading frames.
    """
        
    orfs={} # Initialize an empty dictionary to store the ORFs
    stop_codons = {"TAA", "TGA", "TAG"}
    for gene_id, frames in rf.items(): # Iterate over the sequences in the dictionary
        sub_orfs = {} # Initialize an empty dictionary to store ORFs for the sequence
        for frame_id, frame in frames.items(): # Iterate over the forward reading frames
            orf_sequence = ""
            
            for i in range(0, len(frame), 3): # Iterate over the reading frame to find ORFs by codon
                codon = frame[i:i+3]
                
                if codon == "ATG" or orf_sequence: # Check if the codon is a start codon or if an ORF has started
                    orf_sequence += codon
                    if codon in stop_codons: # Check if the codon is a stop codon
                        sub_orfs[f"O_{frame_id}"] = orf_sequence # Store the ORF in the dictionary
                        break
        
        
                sub_orfs.setdefault(f"O_{frame_id}", "") # Ensure an entry exists for frames without valid ORFs
        orfs[gene_id] = sub_orfs
    return orfs