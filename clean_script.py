# Define the input and output file paths
input_file = 'formatted_filtered_combined_protein_models.fasta.gz'
output_file = 'formatted_filtered_combined_protein_models_cleaned.fasta.gz'

# Function to clean the header
def clean_header(header):
    # Split the header by '|'
    parts = header.split('|')
    
    # Check if the header contains at least 3 parts
    if len(parts) >= 3:
        # Combine the first three parts with '|'
        cleaned_header = '|'.join(parts[:3])
        return cleaned_header
    else:
        # If the header doesn't have enough parts, return it unchanged
        return header

# Open the input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Initialize variables
    header = None

    # Iterate through the lines in the input file
    for line in infile:
        if line.startswith('>'):
            # This is a header line
            header = line.strip()
            # Clean the header
            cleaned_header = clean_header(header)
            # Write the cleaned header to the output file
            outfile.write(cleaned_header + '\n')
        else:
            # This is a sequence line, so copy it to the output file
            outfile.write(line)

# Print a message indicating that the cleaning is complete
print("Headers have been cleaned and saved to", output_file)
