import os

IN_DIR = "RU10_Dorn"
NUM_PACKETS_PER_LINE = 1

if __name__ == "__main__":
    i = 0
    text = ""
    fname = ""
    with open("../.OUTFILES/RU10_output_102k.fasta", "w") as f:
        for file in sorted(os.listdir(IN_DIR)):
            inputfile = open(IN_DIR + "/" + file, "r")
            text += inputfile.read()
            fname += "_" + str(file.split(".")[0])
            if i % NUM_PACKETS_PER_LINE == 1 or NUM_PACKETS_PER_LINE == 1:
                f.write(">" + fname[1:] + "\n")
                f.write(text + "\n")
                text = ""
                fname = ""
            i += 1
            if i > 101999:
                break
            inputfile.close()
