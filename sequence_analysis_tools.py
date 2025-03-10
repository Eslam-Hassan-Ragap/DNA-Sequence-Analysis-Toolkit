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