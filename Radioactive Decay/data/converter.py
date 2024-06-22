# Open the input file
with open('datatest.txt', 'r') as file:
    lines = file.readlines()

# Extract second and third columns
data = [(line.split()[1], line.split()[2]) for line in lines]

# Write the extracted data to a new text file
with open('output.txt', 'w') as file:
    for item in data:
        file.write(f"{item[0]}\t{item[1]}\n")
