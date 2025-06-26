
def read_fasta(file_path):
    sequences = {}  
    ID = ""  
    Seq = ""  
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()  
                
                if line.startswith('>'): 
                    if ID != "" and Seq != "":
                        sequences[ID] = Seq  
                    ID = line  
                    Seq = ""   
                else:
                    Seq += line  
        if ID != "" and Seq != "":
            sequences[ID] = Seq
    except:
        print("Could not find or read the input file!")
        return {}  
    return sequences


def filter_sequences(sequences, min_length):
    long_sequences = {}
    for name, seq in sequences.items():
        if len(seq) >= min_length:
            long_sequences[name] = seq
    return long_sequences


def write_fasta(sequences, output_file):
    try:
        with open(output_file, 'w') as file:
            for name, seq in sequences.items():
                file.write(name + "\n")  # Write sequence name
                # Write sequence in chunks of 80 characters
                for i in range(0, len(seq), 80):
                    file.write(seq[i:i + 80] + "\n")
    except:
        print("Could not save the output file!")


def min_length():
    while True:
        try:
            length = int(input("Enter minimum sequence length (a positive number): "))
            if length > 0:
                return length
            else:
                print("""Please enter a positive number.""")
        except:
            print("Please enter a valid number.")


if __name__ == "__main__":

 input_file = input("Enter the name of your FASTA file: ")
 output_file = "filtered_sequences.fasta"
 sequences = read_fasta(input_file)
 if sequences:
     min_length = min_length()
     long_sequences = filter_sequences(sequences, min_length)
    
     write_fasta(long_sequences, output_file)
    
     print("\n")
     print("Summary:")
     print("Total sequences read:", len(sequences))
     print("Sequences kept (longer than", min_length, "):", len(long_sequences))
     print("Results saved to:", output_file)
 else:
     print("No sequences were processed.")