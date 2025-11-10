import pandas as pd

# --- Configuration ---
FILE_1 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/NT3-AB-NA_R1_001_loops_LVK_LVY_loop_summary.csv'
#FILE_2 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/NT3-AB-NA_R2_001_loops_LVK_LVY_loop_summary.csv'
FILE_2 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/NT3-AB-NB_R1_001_loops_LVK_LVY_loop_summary.csv'
#FILE_2 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/NT3-AB-NB_R2_001_loops_LVK_LVY_loop_summary.csv'

#FILE_1 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/TIMP3-AB-POS_R1_001_loops_LVK_LVY_loop_summary.csv'
#FILE_2 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/TIMP3-AB-POS_R2_001_loops_LVK_LVY_loop_summary.csv'
#FILE_2 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/TIMP3-AB-NEG_R1_001_loops_LVK_LVY_loop_summary.csv'
#FILE_2 = 'Local/Azenta/Code and Output/FASTQ_files/outputs/TIMP3-AB-NEG_R2_001_loops_LVK_LVY_loop_summary.csv'

# Option 1: By column name (string)
COLUMN_TO_COMPARE = 'Loop_Seq'

# Option 2: By column index (integer)
# COLUMN_TO_COMPARE = 2  # (This would be the 3RD column, as indexing is 0-based)


# Set this to None if your CSV files do NOT have a header row
HEADER_SETTING = 0  # 0 means the first row is a header. None if no header

# --- Header Check ---
# If you use a column *NAME* (like 'User ID'), 
# you MUST have a header row.
if isinstance(COLUMN_TO_COMPARE, str) and HEADER_SETTING is None:
    print(f"Error: You specified a column name ('{COLUMN_TO_COMPARE}')")
    print("but set HEADER_SETTING to None. Please set HEADER_SETTING = 0.")
    raise SystemExit

# --- Logic ---
try:
    # Read *only* the specified column from each file
    # .squeeze() converts the single-column DataFrame into a Series
    col1 = pd.read_csv(
        FILE_1, 
        header=HEADER_SETTING, 
        usecols=[COLUMN_TO_COMPARE]
    ).squeeze()
    
    col2 = pd.read_csv(
        FILE_2, 
        header=HEADER_SETTING, 
        usecols=[COLUMN_TO_COMPARE]
    ).squeeze()

    # Convert the columns to sets to find unique items
    set1 = set(col1)
    set2 = set(col2)

    # --- Calculations ---
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    
    # Calculate Jaccard Similarity
    if len(union) == 0:
        jaccard_similarity = 1.0 if len(intersection) == 0 else 0.0
    else:
        jaccard_similarity = len(intersection) / len(union)

    # Find items unique to each file
    only_in_file1 = set1.difference(set2)
    only_in_file2 = set2.difference(set1)

    # --- Results ---
    print(f"--- Column Similarity Report for '{COLUMN_TO_COMPARE}' ---")
    print(f"\nJaccard Similarity Score: {jaccard_similarity:.4f}")
    print(f" (1.0 = identical sets, 0.0 = no overlap)")

    print(f"\n--- Counts ---")
    print(f"Total unique items in File 1: {len(set1)}")
    print(f"Total unique items in File 2: {len(set2)}")
    print(f"Items in *both* files (Intersection): {len(intersection)}")
    print(f"Total unique items in *all* files (Union): {len(union)}")
    print(f"Items *only* in File 1: {len(only_in_file1)}")
    print(f"Items *only* in File 2: {len(only_in_file2)}")


except FileNotFoundError as e:
    print(f"Error: {e}")
except pd.errors.EmptyDataError:
    print("Error: One of the files is empty.")
except ValueError as e:
    # This often happens if the COLUMN_TO_COMPARE doesn't exist
    print(f"ValueError: {e}")
    print(f"Check if '{COLUMN_TO_COMPARE}' is the correct column name/index.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")