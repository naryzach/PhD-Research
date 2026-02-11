import os
import shutil
import glob

# Paths relative to this script
SOURCE_DIR = r"../Local/BD6_Flow"
DEST_DIR = r"../Local/BD6_Flow_Renamed"

def main():
    # Determine absolute paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    src_path = os.path.join(base_dir, SOURCE_DIR)
    dest_path = os.path.normpath(os.path.join(base_dir, DEST_DIR))

    # Clean destination if it exists
    if os.path.exists(dest_path):
        print(f"Cleaning existing directory: {dest_path}")
        shutil.rmtree(dest_path)
    
    os.makedirs(dest_path)
    print(f"Created directory: {dest_path}")

    # Find FCS files
    fcs_files = glob.glob(os.path.join(src_path, "*.fcs"))
    
    if not fcs_files:
        print(f"No FCS files found in {src_path}")
        return

    print(f"Found {len(fcs_files)} files. Analyzing names...")
    
    # Map: lowercase_target_name -> list of (source_path, original_target_name_without_ext, extension)
    from collections import defaultdict
    name_map = defaultdict(list)
    
    for file_path in fcs_files:
        filename = os.path.basename(file_path)
        
        # 1. Start with stripping well number
        parts = filename.split(' ', 1)
        if len(parts) < 2:
            print(f"Skipping {filename}: Could not parse well number (no space found).")
            continue
            
        new_name_with_ext = parts[1]
        
        # 2. Rename 'NC' to 'Negative Control'
        if new_name_with_ext.startswith("NC"):
            new_name_with_ext = new_name_with_ext.replace("NC", "Negative Control", 1)
        
        base_name, ext = os.path.splitext(new_name_with_ext)
        
        # Key for collision detection (lowercase)
        key = base_name.lower().strip()
        name_map[key].append((file_path, base_name, ext))

    # 3. Process and Copy
    count = 0
    for key, items in name_map.items():
        # Check for collision
        is_collision = len(items) > 1
        
        for i, (src_path, base_name, ext) in enumerate(items):
            final_name = ""
            
            if is_collision:
                # User request: "number them sequentially if they don't already have numbers"
                # If we have a collision, we MUST differentiator them to avoid overwrite.
                # We simply append " 1", " 2", etc.
                # Use the original base_name to preserve casing (e.g. EtOh vs EtOH)
                final_name = f"{base_name} {i+1}{ext}"
                print(f"Collision resolved: {os.path.basename(src_path)} -> {final_name}")
            else:
                final_name = f"{base_name}{ext}"
                
            dest_file_path = os.path.join(dest_path, final_name)
            
            try:
                shutil.copy2(src_path, dest_file_path)
                # print(f"Copied: {final_name}")
                count += 1
            except Exception as e:
                print(f"Error copying {final_name}: {e}")

    print(f"\nProcessing complete. {count} files copied and renamed.")

if __name__ == "__main__":
    main()
