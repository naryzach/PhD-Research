import os
from playwright.sync_api import sync_playwright

cif_path = '/home/ryangustafson/Documents/GitHubProj/PhD-Research/Demonstrations/SharedAssets/models/WT/mmp2_ab_wt.cif'
html_path = 'temp_3dmol.html'
output_path = 'mmp2_active_site.png'

# Read CIF content
with open(cif_path, 'r') as f:
    cif_data = f.read()

# Escape for JS template literal
cif_data_js = cif_data.replace('`', '\\`').replace('$', '\\$')

# Create HTML with 3Dmol.js
html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; background: transparent; overflow: hidden; }}
        #viewer {{ width: 800px; height: 600px; position: relative; }}
    </style>
</head>
<body>
    <div id="viewer"></div>
    <script>
        let viewer = $3Dmol.createViewer("viewer", {{
            defaultcolors: $3Dmol.elementColors.rasmol,
            backgroundColor: 'white' // we will make it transparent in playwright, but 3dmol uses it for fog/bg. Actually 3dmol supports alpha.
        }});
        viewer.setBackgroundColor(0x000000, 0.0); // Transparent
        let cifData = `{cif_data_js}`;
        viewer.addModel(cifData, "cif");
        
        // Style: Cartoon colored by secondary structure, maybe highlight some residues?
        viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'ssJmol', radius: 1.2}}}});
        
        // Since we don't have explicit ZN, let's just show a beautiful view of the protein
        viewer.zoomTo();
        
        // Spin it a bit to a nice angle (optional)
        // viewer.rotate(90, 'x');
        
        viewer.render();
    </script>
</body>
</html>
"""

with open(html_path, 'w') as f:
    f.write(html_content)

def render():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.set_viewport_size({"width": 800, "height": 600})
        
        absolute_path = 'file://' + os.path.abspath(html_path)
        page.goto(absolute_path, wait_until='networkidle')
        
        # Wait a moment for 3Dmol to render WebGL
        page.wait_for_timeout(2000)
        
        # Screenshot with omit_background=True for transparency
        page.locator("#viewer").screenshot(path=output_path, omit_background=True)
        
        browser.close()

if __name__ == '__main__':
    render()
    print(f"Successfully generated {output_path}")
    if os.path.exists(html_path):
        os.remove(html_path)
