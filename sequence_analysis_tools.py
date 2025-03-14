def read_fasta(file_path):
    """Function:
    Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values.

    Params:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: Dictionary with sequence IDs as keys and sequences as values.
    """
    
    with open(file_path, "r") as file_fa:
        genes = {}
        id_gene = None
        for line in file_fa:
            line = line.strip()
            if line.startswith('>'): # Check if the line is a sequence ID
                id_gene = line[1:line.index(" ")] if " " in line else line[1:]
                if id_gene not in genes:
                    genes[id_gene] = ""
            elif id_gene:
                genes[id_gene] += line
    return genes

def num_seq(seq):
    """Function:
    Count the number of sequences in a FASTA file or a dictionary of sequences.

    Args:
        seq (dict or str): Dictionary of sequences or path to a FASTA file.

    Returns:
        int: Number of sequences in the file or dictionary.
    """
    if isinstance(seq, dict): # Check if seq is a dictionary
        return len(seq)
    elif isinstance(seq, str): # Check if seq is a FASTA file path
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
    if isinstance(genes,str):  # Check if genes is a FASTA file path
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
            
            for i in range(0, len(frame) - (len(frame) % 3), 3): # Iterate over the reading frame to find ORFs by codon
                codon = frame[i:i+3]
                
                if codon == "ATG" or orf_sequence: # Check if the codon is a start codon or if an ORF has started
                    orf_sequence += codon
                    if codon in stop_codons: # Check if the codon is a stop codon
                        sub_orfs[f"O_{frame_id}"] = orf_sequence # Store the ORF in the dictionary
                        break
        
        
            sub_orfs.setdefault(f"O_{frame_id}", "") # Ensure an entry exists for frames without valid ORFs
        orfs[gene_id] = sub_orfs
    return orfs

def max_len_ORF(reading_frames, num_rf=0, seq_id=None, starting_position=False):
    """
    Finds the longest Open Reading Frame (ORF) in the given set of reading frames.
    
    Args:
        reading_frames (dict): Dictionary where each sequence ID maps to its forward reading frames.
        num_rf (int, optional): Specifies the reading frame (1, 2, or 3). Defaults to 0 (any frame).
        seq_id (str, optional): Specific gene ID to search for the longest ORF. Defaults to None.
        starting_position (bool, optional): If True, returns the starting position of the ORF. Defaults to False.
    
    Returns:
        dict or tuple: Dictionary with the gene ID and its longest ORF length.
                      If `starting_position=True`, returns a tuple (max_len_orf, start_position).
    """
  
    orfs = get_forward_orf(reading_frames) # Get ORFs for the gene sequences

    max_orf_length = 0
    max_orf_info = {}
    longest_orf_name = None
    gene_with_longest_orf = None
    
    if not seq_id:  # Find longest ORF across all ORFs
        for gene_id, orf_dict in orfs.items():
            if num_rf:  # If a specific reading frame is requested
                orf_name = f"O_f_reading_frame_{num_rf}"
                orf_seq = orf_dict.get(orf_name, "")
                if len(orf_seq) > max_orf_length:
                    max_orf_length = len(orf_seq)
                    gene_with_longest_orf = gene_id
                    longest_orf_name = orf_name
            else:  # Check all frames
                longest_orf_seq = max(orf_dict.values(), key=len, default="") # Find longest ORF in the gene
                if len(longest_orf_seq) > max_orf_length:
                    max_orf_length = len(longest_orf_seq)
                    gene_with_longest_orf = gene_id
                    longest_orf_name = [k for k, v in orf_dict.items() if v == longest_orf_seq][0] # Get ORF name
    
    else:  # If a specific gene ID is provided
        if seq_id not in orfs:
            return {}  # Return empty if gene ID not found
        
        longest_orf_seq = max(orfs[seq_id].values(), key=len, default="") # Find longest ORF in the gene
        max_orf_length = len(longest_orf_seq)
        longest_orf_name = [k for k, v in orfs[seq_id].items() if v == longest_orf_seq][0] # Get ORF name
        gene_with_longest_orf = seq_id

    if not gene_with_longest_orf:
        return {}  # No valid ORF found

    # Prepare the result dictionary
    max_orf_info = {gene_with_longest_orf: {longest_orf_name: max_orf_length}}

    # Find starting position if requested
    if starting_position:
        reading_frame_seq = reading_frames[gene_with_longest_orf].get(f"f_reading_frame_{longest_orf_name[-1]}", "")
        start_pos = reading_frame_seq.find(orfs[gene_with_longest_orf][longest_orf_name]) + 1  # 1-based index
        return max_orf_info, start_pos

    return max_orf_info

def count_repeates(genes,length=1,repeat=None):
    """Function:
    Count the number of repeats of a specific length or a specific sequence in each sequence of a dictionary

    Args:
        genes (dict): Dictionary where keys are sequence IDs and values are DNA sequences.
        length (int, optional): length of the repeat. Defaults to 1.  
        repeat (str, optional):  sequence of the repeat. Defaults to None.

    Returns:
        dict: Dictionary where each sequence ID maps to a dictionary of repeats and their counts.
    """
    repeates={}
    if repeat:  # Check if a specific repeat sequence is provided
        length=len(repeat)
        for id_gene,seq_gene in genes.items(): # Iterate over the sequences in the dictionary
            sub_repeates={}
            for rep in range(0,len(genes)-(length-1)): # Iterate over the sequence to find repeats
                if seq_gene[rep:rep+length]==repeat:
                    sub_repeates[seq_gene[rep:rep+length]]=sub_repeates.get(seq_gene[rep:rep+length],0)+1 # Store the repeat and its count
            repeates[id_gene]=sub_repeates
    else: # Find all repeats of a specific length or all lengths
        if length==1: # Find all repeats of all lengths
            for id_gene,seq_gene in genes.items():
                sub_repeates={}
                for Len in range(2,len(seq_gene)+1):      
                    for rep in range(0,len(seq_gene)-(Len-1)):
                        sub_repeates[seq_gene[rep:rep+Len]]=sub_repeates.get(seq_gene[rep:rep+Len],0)+1
                    sub_repeates={k:v for k,v in sub_repeates.items() if v>1} # Store the repeat and its count which is greater than 1
                repeates[id_gene]=sub_repeates
        else:   # Find repeats of a specific length
            for id_gene,seq_gene in genes.items():
                sub_repeates={}
                for rep in range(0,len(seq_gene)-(length-1)):
                    sub_repeates[seq_gene[rep:rep+length]]=sub_repeates.get(seq_gene[rep:rep+length],0)+1
                sub_repeates={k:v for k,v in sub_repeates.items() if v>1}
                repeates[id_gene]=sub_repeates
                
    return repeates