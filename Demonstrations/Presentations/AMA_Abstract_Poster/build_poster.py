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
        
        # Set viewport to 1920x1440 (4:3 poster)
        page.set_viewport_size({"width": 1920, "height": 1440})
        
        # Export PDF
        page.pdf(path='poster_output.pdf', width='48in', height='36in', print_background=True)
        
        # Export PNG
        page.screenshot(path='poster_output.png', full_page=True)
        
        browser.close()

if __name__ == '__main__':
    html_file = render_html()
    generate_pdf_png(html_file)
    print("Successfully generated poster_output.pdf and poster_output.png")
