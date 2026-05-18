# =============================================================
# AMA ABSTRACT POSTER BUILD SCRIPT
# ⚠️  FIXED DIMENSIONS: 48" W × 36" H (landscape, 4:3 ratio)
#     Do NOT change width/height in this file or in style.css
#     without AMA submission approval.
#     At 96 DPI: 4608 × 3456 px
#
# Outputs:
#   poster_output.pdf      — vector PDF, best for printing (48" × 36")
#   poster_output.png      — 96 DPI preview (4608 × 3456 px)
#   poster_hires.png       — 288 DPI print-ready (13824 × 10368 px)  [--hires]
# =============================================================
import argparse
import yaml
from jinja2 import Environment, FileSystemLoader
from playwright.sync_api import sync_playwright
import os

def render_html():
    with open('poster_data.yaml', 'r') as f:
        data = yaml.safe_load(f)
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template('template.html')
    output_html = template.render(data=data)
    with open('output.html', 'w', encoding='utf-8') as f:
        f.write(output_html)
    return 'output.html'

def generate_pdf_png(html_file, hires=False):
    absolute_path = 'file://' + os.path.abspath(html_file)

    with sync_playwright() as p:
        # ── Standard render: PDF + 96 DPI preview (original approach) ────────
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(absolute_path, wait_until='networkidle')
        # ⚠️  AMA POSTER FIXED DIMENSIONS — DO NOT CHANGE
        page.set_viewport_size({"width": 4608, "height": 3456})
        page.pdf(path='poster_output.pdf', width='48in', height='36in', print_background=True)
        page.screenshot(path='poster_output.png', full_page=True)
        browser.close()

        # ── High-res PNG: same page load with device_scale_factor=2 ──────────
        # Kept as a separate session so it doesn't disturb the standard render.
        # deviceScaleFactor=2 → 192 DPI (9216×6912 px); sufficient for large-format print.
        if hires:
            print("Capturing 192 DPI PNG (9216 × 6912 px) …")
            browser = p.chromium.launch(headless=True)
            ctx = browser.new_context(
                viewport={"width": 4608, "height": 3456},
                device_scale_factor=2,
            )
            page = ctx.new_page()
            page.goto(absolute_path, wait_until='networkidle')
            page.wait_for_timeout(4000)
            page.screenshot(path='poster_hires.png')
            browser.close()
            print("  ✓  poster_hires.png  (192 DPI, 9216 × 6912 px)")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--hires', action='store_true',
        help='Also export poster_hires.png at 288 DPI (13824×10368 px). Takes ~60 s extra.'
    )
    args = parser.parse_args()

    html_file = render_html()
    generate_pdf_png(html_file, hires=args.hires)

    outputs = "output.html, poster_output.pdf, and poster_output.png"
    if args.hires:
        outputs += ", and poster_hires.png"
    print(f"Successfully generated {outputs}")
