import pandas as pd
import pickle

# === Settings ===
input_file = "out.csv"       # change to your actual CSV filename
output_file = "output.csv"
search_sequence = "YYCAR" #"DTAVYYCAR"
columns_to_check = ['Frame fwd 1', 'Frame rev 1', 
                        'Frame fwd 2', 'Frame rev 2',
                        'Frame fwd 3', 'Frame rev 3']

# === Step 1: Load reference strings and labels ===
file_motifs = "motifs.pkl"
file_labels = "labels.pkl"
reference_motifs = None
with open(file_motifs, 'rb') as file:
    reference_motifs = pickle.load(file)
reference_labels = None
with open(file_labels, 'rb') as file:
    reference_labels = pickle.load(file)


# === Step 2: Load input data, data starts on row 2 ===
header_row = pd.read_csv(input_file, nrows=1).columns.tolist()
df = pd.read_csv(input_file, skiprows=1, header=None)
df.columns = header_row

results = []

# === Step 3: Scan ===
for idx, row in df.iterrows():
    for col_frame in columns_to_check:
        try:
            cell_value = str(row[col_frame])
        except KeyError:
            continue  # skip missing column
        if pd.notna(cell_value) and search_sequence in cell_value:
            upstream = cell_value.split(search_sequence)[0]
            best_label = ''
            for motifs, label in zip(reference_motifs, reference_labels):
                if motifs:
                    count_matches = sum(1 for motif in motifs if motif in upstream)
                    if count_matches / len(motifs) >= 0.6:
                        best_label = label
                        break  # stop at first match
            new_row = {
                'Sequence': row['Sequence'],
                'Count': row['Count'],
                'Matched_AA_Seq': cell_value,
                'Reading_Frame': col_frame,
                'Reference_Label': best_label
            }
            results.append(new_row)

# === Step 4: Save output ===
if results:
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_file, index=False)
    print(f"Done! Found {len(results)} matches and saved to {output_file}")
else:
    print("No matches found.")
