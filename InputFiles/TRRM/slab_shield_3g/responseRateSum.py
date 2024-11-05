import json

# Load the JSON data from a file
file_name = '/Users/ru/SCONE/InputFiles/TRRM/slab_shield_3g/slab_shield_3g.json'  # Replace 'XXX' with your actual file name

with open(file_name, 'r') as f:
    data = json.load(f)

# Access the 'responseRate' field
response_rates = data['responseRate']['responseRate']

# Initialize a variable to store the sum
total_sum = 0.0

# Iterate over the nested lists to extract first elements and sum them
for outer_list in response_rates:
    for pair in outer_list:
        first_element = pair[0]
        total_sum += first_element

# Print the total sum
print("Total sum of first elements:", total_sum)
