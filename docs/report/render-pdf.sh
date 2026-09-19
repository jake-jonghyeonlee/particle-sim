#!/bin/zsh
set -eu
cd "${0:A:h}"
# Set REPORT_PYTHON to a Python with playwright and pymupdf installed.
report_python="${REPORT_PYTHON:-work/pdf-tools/bin/python}"
"$report_python" render_pdf.py
