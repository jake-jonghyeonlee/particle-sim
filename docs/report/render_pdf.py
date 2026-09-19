"""Render locally with Chrome; reject overflow, missing assets or lost PDF text.

Install report-only dependencies with: uv pip install --python <python> playwright pymupdf
Run with that Python. Chrome must be installed; no browser download is needed.
"""
from collections import Counter
import json
from pathlib import Path
import re
import unicodedata

import pymupdf
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent

def normalized(text):
    return re.sub(r'\s+', '', unicodedata.normalize('NFKC', text)).replace('\u200b', '')

with sync_playwright() as driver:
    browser = driver.chromium.launch(channel='chrome')
    page = browser.new_page()
    page.emulate_media(media='print')
    page.goto((ROOT / 'report.html').as_uri(), wait_until='networkidle')
    page.evaluate('document.fonts.ready')
    assets = page.evaluate('''() => [...document.images].map(i => ({
        src:i.getAttribute('src'), loaded:i.complete && i.naturalWidth > 0}))''')
    assert all(i['loaded'] for i in assets), assets
    layout = page.evaluate('''() => [...document.querySelectorAll('.page')].map((s,i) => {
        const footer=s.querySelector('footer').getBoundingClientRect();
        const body=[...s.children].filter(e=>e.tagName!=='FOOTER');
        const rect=s.getBoundingClientRect();
        return {page:i+1, text:s.innerText,
            footer_clearance_px:footer.top-Math.max(...body.map(e=>e.getBoundingClientRect().bottom)),
            horizontal_overflow:[...s.querySelectorAll('*')].filter(e=>{
                const r=e.getBoundingClientRect();return r.left<rect.left-1 || r.right>rect.right+1;
            }).map(e=>e.tagName)};
    })''')
    assert all(s['footer_clearance_px'] >= 4 and not s['horizontal_overflow'] for s in layout), layout
    temporary = ROOT / 'work/report-candidate.pdf'
    temporary.parent.mkdir(exist_ok=True)
    page.pdf(path=str(temporary), prefer_css_page_size=True, print_background=True,
             display_header_footer=False)
    browser.close()

with pymupdf.open(temporary) as pdf:
    assert len(pdf) == len(layout), (len(pdf), len(layout))
    for expected, actual in zip(layout, pdf, strict=True):
        missing = Counter(normalized(expected.pop('text'))) - Counter(normalized(actual.get_text()))
        assert not missing, (expected['page'], missing)
        expected['missing_text_characters'] = sum(missing.values())
        expected['pdf_images'] = len(actual.get_images())
        for word in actual.get_text('words'):
            assert actual.rect.contains(pymupdf.Rect(word[:4])), (expected['page'], word)
        actual.get_pixmap(matrix=pymupdf.Matrix(1.25, 1.25)).save(
            ROOT / f"work/page-{expected['page']:02}.png")
result = {'pages':len(layout), 'images':assets, 'layout':layout}
(ROOT / 'evidence/pdf-validation.json').write_text(json.dumps(result, ensure_ascii=False, indent=2))
temporary.replace(ROOT / 'particle-sim-report.pdf')
print(f"Verified {len(layout)} pages, {len(assets)} images; no text loss or overflow.")
