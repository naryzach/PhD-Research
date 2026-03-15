import os
import re

def find_streamlit_secrets(directory):
    """
    Scans Python files for st.secrets and st.connection usage.
    Handles:
    - st.secrets["a"]["b"]
    - st.secrets.a.b
    - st.secrets.get("a")
    - email_secrets = st.secrets["email"]; email_secrets.get("sender")
    - "key" in st.secrets["section"] (Membership test as a key source)
    - st.connection("name") -> [connections.name]
    """
    start_pattern = re.compile(r'(st\.secrets|st\.connection)')
    # Alias pattern: var = st.secrets["section"] or var = st.secrets.section
    alias_pattern = re.compile(r'([a-zA-Z0-9_]+)\s*=\s*st\.secrets(?:\[["\']([^"\']+)["\']\]|\.([a-zA-Z0-9_]+))')
    
    results = {}

    for root, dirs, files in os.walk(directory):
        if '.git' in dirs: dirs.remove('.git')
        if '.venv' in dirs: dirs.remove('.venv')
        
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        content = f.read()
                        
                        folder_results = []
                        
                        # Pass 1: Find aliases
                        aliases = {}
                        for match in alias_pattern.finditer(content):
                            var_name = match.group(1)
                            section = match.group(2) or match.group(3)
                            if section: aliases[var_name] = [section]
                        
                        # Pass 2: Detect st.secrets/st.connection chains
                        for match in start_pattern.finditer(content):
                            start_type = match.group(0)
                            pos = match.end()
                            
                            if start_type == 'st.secrets':
                                chain = []
                                while pos < len(content):
                                    bracket_match = re.match(r'\s*\[["\']([^"\']+)["\']\]', content[pos:])
                                    if bracket_match:
                                        chain.append(bracket_match.group(1))
                                        pos += bracket_match.end()
                                        continue
                                    
                                    dot_match = re.match(r'\s*\.([a-zA-Z0-9_]+)', content[pos:])
                                    if dot_match:
                                        attr = dot_match.group(1)
                                        if attr == 'get':
                                            get_match = re.match(r'\s*\(\s*["\']([^"\']+)["\']', content[pos+dot_match.end():])
                                            if get_match:
                                                chain.append(get_match.group(1))
                                                pos += dot_match.end() + get_match.end()
                                                continue
                                        else:
                                            chain.append(attr)
                                            pos += dot_match.end()
                                            continue
                                    break
                                if chain: folder_results.append(chain)
                            
                            elif start_type == 'st.connection':
                                conn_match = re.match(r'\s*\(\s*["\']([^"\']+)["\']', content[pos:])
                                if conn_match:
                                    conn_name = conn_match.group(1)
                                    folder_results.append(['connections', conn_name, 'url'])
                                    pos += conn_match.end()
                        
                        # Pass 3: Look for aliases usage
                        for var_name, base_chain in aliases.items():
                            # Find all var_name.get("key") or var_name["key"]
                            # var_name.get(...)
                            get_pat = re.compile(re.escape(var_name) + r'\.get\(["\']([^"\']+)["\']')
                            for k in get_pat.findall(content):
                                folder_results.append(base_chain + [k])
                            
                            # var_name[...]
                            bracket_pat = re.compile(re.escape(var_name) + r'\[["\']([^"\']+)["\']\]')
                            for k in bracket_pat.findall(content):
                                folder_results.append(base_chain + [k])

                        # Pass 4: Membership tests like "key" in st.secrets["section"]
                        # regex for: ["']key["'] in st.secrets[...] or st.secrets.section
                        member_pat = re.compile(r'["\']([^"\']+)["\']\s+in\s+st\.secrets((?:\[["\'](?:[^"\']+)["\']\]|\.[a-zA-Z0-9_]+)+)')
                        for key, chain_str in member_pat.findall(content):
                            # Extract keys from chain_str
                            base_chain = re.findall(r'\[["\']([^"\']+)["\']\]|\.([a-zA-Z0-9_]+)', chain_str)
                            # base_chain is list of tuples (bracket_val, dot_val)
                            flat_chain = [b or d for b, d in base_chain]
                            folder_results.append(flat_chain + [key])

                        if folder_results:
                            if root not in results:
                                results[root] = []
                            results[root].extend(folder_results)
                except Exception as e:
                    print(f"Error reading {path}: {e}")
    
    return results

def merge_trees(existing, new):
    """Deep merges two trees."""
    for key, value in new.items():
        if key in existing:
            if isinstance(existing[key], dict) and isinstance(value, dict):
                merge_trees(existing[key], value)
            elif isinstance(value, dict):
                # Upgrade flat value to dict if it looks like a placeholder
                if isinstance(existing[key], str) and "[YOUR" in existing[key].upper():
                    existing[key] = value
        else:
            existing[key] = value

def parse_simple_toml(content):
    tree = {}
    current_node = tree
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('[') and line.endswith(']'):
            section_name = line[1:-1]
            parts = section_name.split('.')
            current_node = tree
            for part in parts:
                if part not in current_node:
                    current_node[part] = {}
                elif not isinstance(current_node[part], dict):
                    # Conflict: flat becomes dict
                    current_node[part] = {"_value": current_node[part]}
                current_node = current_node[part]
        elif '=' in line:
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip().strip('"').strip("'")
            current_node[key] = val
    return tree

def write_simple_toml(tree, prefix=''):
    lines = []
    
    # 1. Flat keys at this level
    for key in sorted(tree.keys()):
        val = tree[key]
        if not isinstance(val, dict):
            lines.append(f'{key} = "{val}"')
    
    # 2. Sub-sections
    for section in sorted(tree.keys()):
        sub_tree = tree[section]
        if isinstance(sub_tree, dict):
            full_section = f"{prefix}.{section}" if prefix else section
            if lines: lines.append("")
            lines.append(f"[{full_section}]")
            # Recursive call for children
            # But we want to keep them in the same block if it's a sub-section of TOML
            # Actually, standard TOML writer for deep nesting is better done flatly for simplicity
            lines.extend(format_toml_node(sub_tree, full_section))
    return "\n".join(lines) + "\n"

def format_toml_node(node, current_prefix):
    lines = []
    # Dump flat keys
    for key in sorted(node.keys()):
        val = node[key]
        if not isinstance(val, dict):
            lines.append(f'{key} = "{val}"')
    
    # Dump sub-sections
    for section in sorted(node.keys()):
        sub_tree = node[section]
        if isinstance(sub_tree, dict):
            full_prefix = f"{current_prefix}.{section}"
            lines.append("")
            lines.append(f"[{full_prefix}]")
            lines.extend(format_toml_node(sub_tree, full_prefix))
    return lines

def setup_secrets(results):
    """
    Creates/Updates .streamlit/secrets.toml files.
    """
    for folder, secret_chains in results.items():
        st_dir = os.path.join(folder, '.streamlit')
        if not os.path.exists(st_dir):
            os.makedirs(st_dir)
            print(f"Created directory: {st_dir}")
            
        secrets_file = os.path.join(st_dir, 'secrets.toml')
        gitignore_file = os.path.join(st_dir, '.gitignore')
        
        if not os.path.exists(gitignore_file):
            with open(gitignore_file, 'w') as f:
                f.write('*\n')
            print(f"Created .gitignore in: {st_dir}")

        # Build new tree
        new_tree = {}
        # Sort chains by length so shorter ones are processed first, but can be overridden by longer ones
        for chain in sorted(secret_chains, key=len):
            current = new_tree
            for i, key in enumerate(chain):
                if i == len(chain) - 1:
                    # Last element: set value if not already a dict
                    if key not in current or not isinstance(current[key], dict):
                        current[key] = f"[YOUR {key.upper()}]"
                else:
                    # Intermediate element: must be a dict
                    if key not in current or not isinstance(current[key], dict):
                        current[key] = {}
                    current = current[key]

        # Read existing
        existing_tree = {}
        if os.path.exists(secrets_file):
            with open(secrets_file, 'r') as f:
                existing_tree = parse_simple_toml(f.read())
        
        # Merge
        merge_trees(existing_tree, new_tree)
        
        # Write back
        new_content = write_simple_toml(existing_tree)
        
        # Strip extra leading/trailing newlines for cleanliness
        new_content = new_content.strip() + "\n"
        
        with open(secrets_file, 'w') as f:
            f.write(new_content)
        
        print(f"Updated/Verified secrets: {secrets_file}")

if __name__ == "__main__":
    import os
    import re
    repo_root = os.path.abspath(os.path.dirname(__file__))
    print(f"Scanning repository: {repo_root}")
    secret_map = find_streamlit_secrets(repo_root)
    if not secret_map:
        print("No st.secrets usage found.")
    else:
        setup_secrets(secret_map)
        print("\nSetup complete. Please populate the secrets.toml files with actual credentials.")
