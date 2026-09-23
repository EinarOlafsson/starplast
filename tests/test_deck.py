"""The published guide must remain navigable and agree across portable formats."""
import json
from pathlib import Path
import re
import zipfile
from xml.etree import ElementTree as ET

from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
DECK=ROOT/'docs/deck'


def test_every_slide_has_images_transcript_and_wrapping_github_navigation():
    deck=json.loads((DECK/'slides.json').read_text())
    count=deck['count']
    assert 1<=count<=50 and len(deck['slides'])==count
    for i,slide in enumerate(deck['slides'],1):
        assert (DECK/slide['image']).is_file()
        assert (DECK/slide['thumb']).is_file()
        assert len(slide['text'])>80
        page=(DECK/'pages'/f'{i:02d}.md').read_text()
        previous=count if i==1 else i-1
        next_slide=1 if i==count else i+1
        assert f'({previous:02d}.md)' in page and f'({next_slide:02d}.md)' in page
        assert f'/deck/#{i}' in page
        assert slide['source'].startswith('https://einarolafsson.github.io/starplast/')


def test_pdf_and_editable_powerpoint_preserve_slide_count_and_titles():
    deck=json.loads((DECK/'slides.json').read_text())
    pdf=PdfReader(DECK/'starplast_deck.pdf')
    assert len(pdf.pages)==deck['count']
    assert pdf.metadata.author=='Einar Olafsson'
    with zipfile.ZipFile(DECK/'starplast_deck.pptx') as archive:
        slides=[p for p in archive.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml',p)]
        assert len(slides)==deck['count']
        for i,(page,slide) in enumerate(zip(pdf.pages,deck['slides']),1):
            assert slide['title'] in ' '.join(page.extract_text().split())
            tree=ET.fromstring(archive.read(f'ppt/slides/slide{i}.xml'))
            text=' '.join(n.text or '' for n in tree.iter('{http://schemas.openxmlformats.org/drawingml/2006/main}t'))
            assert slide['title'] in text  # Native editable text, not a flattened slide image.


def test_viewer_is_offline_capable_and_readme_uses_portable_links():
    viewer=(DECK/'index.html').read_text()
    assert '__SLIDES__' not in viewer
    assert 'const DECK = {' in viewer and 'fetch(' not in viewer
    assert 'starplast_deck.pdf' in viewer and 'starplast_deck.pptx' in viewer
    assert 'ArrowRight' in viewer and 'touchend' in viewer and 'showModal' in viewer
    readme=(ROOT/'README.md').read_text()
    assert 'https://einarolafsson.github.io/starplast/deck/' in readme
    assert 'https://raw.githubusercontent.com/EinarOlafsson/starplast/main/docs/deck/slides/slide_01.jpg' in readme
