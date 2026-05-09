# AMA Abstract Poster

This directory contains the files needed to generate a 48" W x 36" H academic poster from a YAML configuration file. 
The poster is built using HTML/CSS (with Jinja2 templating) and exported to PDF/PNG using Playwright.

## Required Dependencies

To build the poster, you need the following Python packages installed in your environment (e.g., the project's `.venv`):

- `jinja2`
- `pyyaml`
- `playwright`

You can install them via pip:
```bash
pip install -r requirements.txt
```

### System Dependencies
Playwright also requires the Chromium browser binary and its underlying Linux system C-libraries. To install these, run:

```bash
# Install the browser binary (can be done without sudo)
playwright install chromium

# Install system dependencies (REQUIRES SUDO)
playwright install-deps
```

## How to Build
Once dependencies are satisfied, generate the `poster_output.pdf` and `poster_output.png` by running:

```bash
python build_poster.py
```
