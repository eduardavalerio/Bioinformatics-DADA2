input_file = "barcodes_lib1.txt"
output_file = "barcodes_lib1.fasta"

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        line = line.strip()
        if line.startswith(">"):  # Header line (barcode)
            f_out.write(line + "\n")
        else:  # Sequence line (aaaa:aaaa)
            sequence_part = line.split(":")[0]  # Take the part before the colon
            f_out.write(sequence_part.upper() + "\n")  # Convert to uppercase

# Output 
# >barcode01
# AAAAAAAA
# >barcode02
# AAAAACCC
# >barcode03
# GCAATTTT
