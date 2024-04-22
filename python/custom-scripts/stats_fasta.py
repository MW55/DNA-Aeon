# This function should check fastaFiles properties

# in the checks, it should also check the GC content of the sequence(between 40-60%)

def stats_fasta(fasta_file):
    with open(fasta_file, 'r') as file :
        content = file.read()
        entries = content.split('>') 
        headers = []
        sequences = []
        test = False
        for entry in entries:
            lines = entry.split('\n',1)
            if len(lines) < 2:
                continue
            headers.append(lines[0])

            sequences.append(lines[1])
    return (headers, sequences)
                               
def main():
    total_elements = 0
    gc_content = 0 
    hd_count = 0
    sq_count = 0    
    length = 0
    fasta_file = './data/encoded.fasta'
    headers, sequences = stats_fasta(fasta_file)
    print("done")
    for header in headers:
        print(header)
        hd_count += 1
    for sequence in sequences:
        print(sequence)
        length += len(sequence)
        sq_count += 1
        for base in sequence:
            if base == 'G' or base == 'C':
                gc_content += 1
            total_elements += 1
    
    print(f"GC content: {gc_content/total_elements*100:.2f}% [40-60%]")
    print("length: ", length/sq_count)
    print("headers: ", hd_count)
    print("sequences: ", sq_count)
    print("total bases: ", total_elements)
    print("done")    
    
if __name__ == "__main__":
    main()