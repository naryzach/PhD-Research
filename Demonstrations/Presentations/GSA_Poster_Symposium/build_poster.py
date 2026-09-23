"""Build the Graduate Poster Symposium poster from poster_data.yaml.

Outputs poster_output.pdf (48 x 36 in, vector) and poster_output.png (96 dpi preview).
Add --hires for a 192 dpi PNG. Size is an assumption; change WIDTH_IN/HEIGHT_IN here and in
style.css together if the symposium specifies a different size.
"""
import argparse
import os
import yaml
from jinja2 import Environment, FileSystemLoader
from playwright.sync_api import sync_playwright

WIDTH_IN, HEIGHT_IN, DPI = 48, 36, 96
W, H = WIDTH_IN * DPI, HEIGHT_IN * DPI
HERE = os.path.dirname(os.path.abspath(__file__))


def render_html():
    with open(os.path.join(HERE, "poster_data.yaml"), encoding="utf-8") as f:
        data = yaml.safe_load(f)
    env = Environment(loader=FileSystemLoader(HERE))
    html = env.get_template("template.html").render(data=data)
    out = os.path.join(HERE, "output.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    return out


def export(html_file, hires=False):
    url = "file:///" + html_file.replace("\\", "/")
    with sync_playwright() as p:
        b = p.chromium.launch()
        page = b.new_page(viewport={"width": W, "height": H})
        page.goto(url, wait_until="load")
        page.wait_for_timeout(1500)
        page.pdf(path=os.path.join(HERE, "poster_output.pdf"), width=f"{WIDTH_IN}in",
                 height=f"{HEIGHT_IN}in", print_background=True)
        page.screenshot(path=os.path.join(HERE, "poster_output.png"))
        overflow = page.evaluate("""() => {
            const out = [];
            document.querySelectorAll('.col').forEach((c, i) => {
                out.push({col: i + 1, scrollH: c.scrollHeight, clientH: c.clientHeight});
            });
            return out;
        }""")
        b.close()
        if hires:
            b = p.chromium.launch()
            ctx = b.new_context(viewport={"width": W, "height": H}, device_scale_factor=2)
            pg = ctx.new_page()
            pg.goto(url, wait_until="load")
            pg.wait_for_timeout(2500)
            pg.screenshot(path=os.path.join(HERE, "poster_hires.png"))
            b.close()
    for o in overflow:
        state = "OVERFLOWS" if o["scrollH"] > o["clientH"] + 2 else "fits"
        print(f"column {o['col']}: content {o['scrollH']} px vs available {o['clientH']} px ({state})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--hires", action="store_true")
    args = ap.parse_args()
    export(render_html(), hires=args.hires)
    print("built poster_output.pdf and poster_output.png")
