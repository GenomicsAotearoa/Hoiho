# Read genome names from the file
with open("seabird_genomes", "r") as f:
    content = f.read()

# Split the content into individual genome names
genome_list = content.split()

# Filter out the entries starting with "Anc"
filtered_genomes = [g for g in genome_list if not g.startswith("Anc")]

# Print the result
print(",".join(filtered_genomes))
