import os
import pandas as pd
import sys

# Directory to load data
FunTB_dir = os.getcwd()

# Input arguments - python Sample_Grouping_Creation.py Metadata File.csv
metadata_file = sys.argv[1]  # Metadata file (.csv)
file_path = os.path.join(FunTB_dir, 'Metadata_files', metadata_file)  # Metadata file path
metadata_df = pd.read_csv(file_path, encoding='utf-8')  # Metadata file load

Selected_columns = sys.argv[2:]

# General functions

def filter_dataframe(df, filter_dict):
    """
    Filter a DataFrame based on specified filter conditions.
    Parameters:
        - df: DataFrame to be filtered
        - filter_dict: Dictionary containing column names as keys and filter conditions as values
    Returns:
        - Filtered DataFrame
    """
    for column, value in filter_dict.items():
        if isinstance(value, str):
            df = df[df[column] == value]
        else:
            df = df[df[column] == int(value)]
    return df

def get_filters_conditions(df, group, column_name, filter_conditions):
    """
    Prompts the user to enter a value for a column and updates the filter conditions.
    """
    valid_values = set(df[column_name].values)
    while True:
        value = input(f"{column_name} value for group '{group}' ({valid_values}): ")
        try:
            if value.isdigit() and int(value) in valid_values:
                filter_conditions[column_name] = int(value)
                break
            elif value in valid_values:
                filter_conditions[column_name] = value
                break
            else:
                print(f"Invalid value. Choose from {valid_values}.")
        except ValueError:
            print(f"Invalid input, please enter a valid value.")
    return filter_conditions

def save_group_Ids(file_name, samples_Ids):
    with open(file_name + ".txt", "w") as output:
        # Write as a Python list with quoted strings
        output.write("[\n" + ",\n".join([f"    '{id}'" for id in samples_Ids]) + "\n]")

def get_groups_names(number_of_groups):
    """
    Prompts the user to enter names for a specified number of groups.
    """
    return [input(f"Set name of group {i+1}: ") for i in range(number_of_groups)]

# Main program
number_of_groups = int(input('Set number of groups: '))
Groups_names = get_groups_names(number_of_groups)

for group in Groups_names:
    filter_conditions = {}
    for column_name in Selected_columns:
        filter_conditions = get_filters_conditions(metadata_df, group, column_name, filter_conditions)
    
    filtered_data = filter_dataframe(metadata_df, filter_conditions)

    # Save group IDs
    save_file = os.path.join(FunTB_dir, 'Samples_lists_files', group)
    save_group_Ids(save_file, filtered_data.iloc[:, 0].tolist())  # Assuming IDs are in the first column