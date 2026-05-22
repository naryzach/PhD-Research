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
#   poster_hires.png       — 192 DPI print-ready (9216 × 6912 px)   [--hires]
#   poster_digital.pdf     — 16:9 PDF (64" × 36") for digital dir.  [--digital]
#   poster_digital.png     — 16:9 preview                           [--digital]
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

def generate_digital_pdf(html_file):
    """
    Generates poster_digital.pdf: 16:9 canvas (64" × 36" = 6144 × 3456 px).
    The 48"×36" poster is centered with 8" side margins on the same gradient
    background. CSS overrides are injected at runtime — no separate template needed.
    """
    DIGITAL_W_PX = 6144     # 64" × 96 dpi
    POSTER_W_PX  = 4608     # 48" × 96 dpi
    MARGIN_PX    = (DIGITAL_W_PX - POSTER_W_PX) // 2   # 768 px each side

    absolute_path = 'file://' + os.path.abspath(html_file)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(absolute_path, wait_until='networkidle')
        page.set_viewport_size({"width": DIGITAL_W_PX, "height": 3456})

        # Inject print-mode overrides appended last so these !important rules win cascade.
        # html becomes the 64"-wide canvas; body is the 48" poster centered within it.
        page.evaluate(f"""
            const s = document.createElement('style');
            s.textContent = `
                @media print {{
                    html {{
                        width: {DIGITAL_W_PX}px !important;
                        overflow: visible !important;
                        background: linear-gradient(135deg, #f0f9ff 0%, #dbeafe 100%) !important;
                    }}
                    body {{
                        width: {POSTER_W_PX}px !important;
                        height: 3456px !important;
                        margin-left: {MARGIN_PX}px !important;
                        margin-top: 0px !important;
                        transform: none !important;
                    }}
                }}
            `;
            document.head.appendChild(s);
        """)

        page.pdf(
            path='poster_digital.pdf',
            width='64in',
            height='36in',
            print_background=True,
        )
        page.screenshot(path='poster_digital.png', full_page=True)
        browser.close()

    print("  ✓  poster_digital.pdf  (16:9, 64\" × 36\", 8\" side margins)")
    print("  ✓  poster_digital.png  (16:9 preview)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--hires', action='store_true',
        help='Also export poster_hires.png at 192 DPI (9216×6912 px). Takes ~60 s extra.'
    )
    parser.add_argument(
        '--digital', action='store_true',
        help='Also export poster_digital.pdf/.png at 16:9 (64"×36") for digital directory submission.'
    )
    args = parser.parse_args()

    html_file = render_html()
    generate_pdf_png(html_file, hires=args.hires)

    if args.digital:
        print("\nGenerating 16:9 digital version …")
        generate_digital_pdf(html_file)

    outputs = "output.html, poster_output.pdf, and poster_output.png"
    if args.hires:
        outputs += ", poster_hires.png"
    if args.digital:
        outputs += ", poster_digital.pdf, and poster_digital.png"
    print(f"\nSuccessfully generated {outputs}")
