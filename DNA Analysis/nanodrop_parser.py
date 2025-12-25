import os
import sys
import pandas as pd
import argparse
import glob

def find_tsv_files(start_dir):
    """Finds all files ending with _table.tsv in the given directory recursively."""
    files = []
    for root, dirs, filenames in os.walk(start_dir):
        for filename in filenames:
            if filename.endswith("_table.tsv"):
                files.append(os.path.join(root, filename))
    return files

def parse_nanodrop_file(file_path):
    """Parses a NanoDrop _table.tsv file."""
    print(f"Parsing: {file_path}")
    try:
        # Pre-scan file to find the data start line (line with many columns)
        data_start_line = 0
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                if not line.strip(): continue
                cols = line.split('\t')
                if len(cols) > 50: # Arbitrary threshold, spectra has > 2000
                    data_start_line = i
                    break
        
        # Read the file skipping lines
        df = pd.read_csv(file_path, sep='\t', header=None, skiprows=data_start_line)
        
        # Determine number of columns
        num_cols = df.shape[1]
        
        # Expected metadata columns: 6
        # 0: DateTime
        # 1: Sample ID
        # 2: A260/A280
        # 3: A260/A230
        # 4: A260
        # 5: A280
        # 6+: Spectra starting at 220.0 nm with 0.5 nm step
        
        meta_cols = ['DateTime', 'Sample ID', 'A260/A280', 'A260/A230', 'A260', 'A280']
        
        if num_cols < 6:
            print(f"Warning: File {file_path} has fewer columns ({num_cols}) than expected.")
            return None

        # Rename the known columns
        rename_map = {i: name for i, name in enumerate(meta_cols)}
        df = df.rename(columns=rename_map)
        
        # Generate wavelength headers for the rest
        # Start 220.0, step 0.5
        wavelength_cols = {}
        current_wl = 220.0
        for i in range(6, num_cols):
            wavelength_cols[i] = f"{current_wl:.1f}"
            current_wl += 0.5
            
        df = df.rename(columns=wavelength_cols)
        
        # Ensure numeric columns for calculations
        cols_to_numeric = ['A260', 'A280', 'A260/A280', 'A260/A230']
        for col in cols_to_numeric:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        # Drop rows where A260 or A280 is NaN (removes potential header rows parsed as data)
        df = df.dropna(subset=['A260', 'A280'])

        # Add a calculated Concentration column (assuming dsDNA factor 50)
        df['Approx_Conc_dsDNA_ng_uL'] = df['A260'] * 50.0
        
        return df
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Parse and display NanoDrop data.")
    parser.add_argument("path", help="Path to file or directory to search.")
    # Removed manual flags, now automatic
    args = parser.parse_args()
    
    path = args.path
    if os.path.isfile(path):
        files = [path]
        output_dir = os.path.dirname(path)
    elif os.path.isdir(path):
        files = find_tsv_files(path)
        output_dir = path
    else:
        print(f"Error: {path} is not a valid file or directory.")
        return

    if not files:
        print("No _table.tsv files found.")
        return

    all_data = []
    
    for f in files:
        df = parse_nanodrop_file(f)
        if df is not None:
            df['Source_File'] = os.path.basename(f)
            all_data.append(df)
    
    if not all_data:
        print("No valid data parsed.")
        return
        
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Define columns to display in summary AND save to CSV
    # Using the "prior terminal summary" columns (richer set)
    # excluding raw spectra columns from this list
    display_cols = ['Sample ID', 'DateTime', 'A260', 'A280', 'A260/A280', 'A260/A230', 'Approx_Conc_dsDNA_ng_uL', 'Source_File']
    
    # Format for terminal display
    pd.options.display.float_format = '{:.2f}'.format
    print("\nSummary Data:")
    # Check if columns exist before printing (paranoia check)
    valid_cols = [c for c in display_cols if c in combined_df.columns]
    print(combined_df[valid_cols].to_string(index=False))
    
    # Save Summary CSV (Only the summary columns, not the whole dataframe)
    output_csv = os.path.join(output_dir, "nanodrop_summary.csv")
    combined_df[valid_cols].to_csv(output_csv, index=False)
    print(f"\nSaved summary to {output_csv}")
    
    # Plotting
    try:
        import matplotlib.pyplot as plt
        
        # Identify wavelength columns (numeric column names)
        # Robustly exclude metadata columns
        exclude_cols = set(display_cols)
        
        wl_cols = []
        for col in combined_df.columns:
            if col in exclude_cols:
                continue
            try:
                # Wavelengths data stored as column names like "220.0"
                # Float conversion test
                float(col)
                wl_cols.append(col)
            except ValueError:
                continue
        
        if not wl_cols:
            print("No spectral data found to plot.")
            return

        # Sort wavelengths numerically
        wl_cols.sort(key=float)
        x_vals = [float(c) for c in wl_cols]
        
        plt.figure(figsize=(12, 8))
        
        # Plot all samples
        for index, row in combined_df.iterrows():
            # Ensure y_vals are float
            y_vals = pd.to_numeric(row[wl_cols], errors='coerce')
            
            label = f"{row['Sample ID']}"
            if len(files) > 1:
                label += f" ({row['Source_File']})"
            
            plt.plot(x_vals, y_vals, label=label, linewidth=1, alpha=0.8)
            
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Absorbance (10mm)")
        plt.title(f"Spectra Summary")
        
        # Legend logic
        if len(combined_df) < 30: 
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            print("Legend hidden due to high sample count (>30).")
            
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        output_plot = os.path.join(output_dir, "nanodrop_spectra.png")
        plt.savefig(output_plot, dpi=300)
        print(f"Saved spectra plot to {output_plot}")
        
    except ImportError:
        print("Error: matplotlib is required for plotting.")
    except Exception as e:
        print(f"Error occurred during plotting: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
