# =============================================================
# AMA ABSTRACT POSTER BUILD SCRIPT
# ⚠️  FIXED DIMENSIONS: 48" W × 36" H (landscape, 4:3 ratio)
#     Do NOT change width/height in this file or in style.css
#     without AMA submission approval.
#     At 96 DPI: 4608 × 3456 px
# =============================================================
import yaml
from jinja2 import Environment, FileSystemLoader
from playwright.sync_api import sync_playwright
import os

def render_html():
    # Load data
    with open('poster_data.yaml', 'r') as f:
        data = yaml.safe_load(f)

    # Setup Jinja2 environment
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template('template.html')

    # Render HTML
    output_html = template.render(data=data)
    with open('output.html', 'w') as f:
        f.write(output_html)
    return 'output.html'

def generate_pdf_png(html_file):
    with sync_playwright() as p:
        # Launch browser
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        
        # Open HTML file
        absolute_path = 'file://' + os.path.abspath(html_file)
        page.goto(absolute_path, wait_until='networkidle')
        
        # Set viewport to 4608x3456 (48" x 36" at 96 DPI)
        # ⚠️  AMA POSTER FIXED DIMENSIONS — DO NOT CHANGE
        page.set_viewport_size({"width": 4608, "height": 3456})
        
        # Export PDF — 48in × 36in landscape
        # ⚠️  AMA POSTER FIXED DIMENSIONS — DO NOT CHANGE
        page.pdf(path='poster_output.pdf', width='48in', height='36in', print_background=True)
        
        # Export PNG
        page.screenshot(path='poster_output.png', full_page=True)
        
        browser.close()

if __name__ == '__main__':
    html_file = render_html()
    generate_pdf_png(html_file)
    print("Successfully generated output.html, poster_output.pdf, and poster_output.png")
