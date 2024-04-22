# This function check that the input file and the encoded/decoded version are the same (no error)

def compare_files(file1, file2, chunk_size=1024):
    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        while True:
            chunk1 = f1.read(chunk_size)
            chunk2 = f2.read(chunk_size)
            if chunk1 != chunk2 :
                return False
            if not chunk1 and not chunk2 :
                return True
            
def main():
    file1 = './data/results/D'
    file2 = './data/D'
    if compare_files(file1,file2):
        print("The files are identical.")
    else: 
        print("The file are different.")

if __name__ == "__main__":
    main()
