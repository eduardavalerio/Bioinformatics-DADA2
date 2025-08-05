# Input 
# aaaaaaaa:aaaaaaaa
# aaaaaccc:aaaaaccc
# gcaatttt:gcaatttt

input_file = "tags.txt"    # Your input file
output_file = "tags_lib1.fasta"  # Output FASTA file
start_number = 1            # Set the number of barcode (e.g., 1 for >barcode1, 2 for >barcode2, etc.)

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for i, line in enumerate(f_in, start=start_number):
        line = line.strip()
        if ":" in line:  # Skip lines without colons (if any)
            seq_part = line.split(":")[0].upper()  # Take left of ":", uppercase
            f_out.write(f">tag{i}\n{seq_part}\n")  # Write FASTA entry

# Output 
# >tag01
# AAAAAAAA
# >tag02
# AAAAACCC
# >tag03
# GCAATTTT
