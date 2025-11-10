import pandas as pd
import pickle

def find_unique_motifs(strings):
    n = len(strings)
    motifs = [[] for _ in range(n)]
    
    # Build a map from substring to set of string indices where it appears
    substring_to_strings = dict()
    
    for idx, s in enumerate(strings):
        seen = set()
        for start in range(len(s)):
            for end in range(start + 1, len(s) + 1):
                sub = s[start:end]
                if sub not in seen:
                    seen.add(sub)
                    if sub not in substring_to_strings:
                        substring_to_strings[sub] = set()
                    substring_to_strings[sub].add(idx)
    
    # For each string, find shortest unique substrings
    for idx, s in enumerate(strings):
        unique_substrings = []
        shortest_length = None
        for start in range(len(s)):
            for end in range(start + 1, len(s) + 1):
                sub = s[start:end]
                if substring_to_strings[sub] == {idx}:
                    if shortest_length is None or len(sub) < shortest_length:
                        shortest_length = len(sub)
                        unique_substrings = [sub]
                    elif len(sub) == shortest_length:
                        unique_substrings.append(sub)
        motifs[idx] = unique_substrings
    
    return motifs


# === Load Excel file ===
excel_path = "../Data/Original_Naive_library_Backbone.xlsx" 
df = pd.read_excel(excel_path, header=None)

# Data starts from row 3 to index 2 (0-based)
df = df.iloc[2:7]  # skip first two rows

names = df.iloc[:, 0].astype(str).tolist()     # Column A
strings = df.iloc[:, 2].astype(str).tolist()   # Column C

# === Find motifs ===
motifs = find_unique_motifs(strings)

file_motifs = "motifs.pkl"
with open(file_motifs, 'wb') as file:
    pickle.dump(motifs, file)
file_labels = "labels.pkl"
with open(file_labels, 'wb') as file:
    pickle.dump(names, file)

# === Print nicely ===
for name, string, motif_list in zip(names, strings, motifs):
    print(f"Name: {name}")
    print(f"AA Seq: {string}")
    print(f"Shortest unique motifs: {motif_list}")
    print("-" * 40)