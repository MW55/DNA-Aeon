# This function should check fastaFiles properties

# we read char by char to find the header and the sequence
# we have to accumulate the sequence in array of array
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
    print("length: ", length/sq_count)
    print("headers: ", hd_count)
    print("sequences: ", sq_count)
    print("done")    

if __name__ == "__main__":
    main()