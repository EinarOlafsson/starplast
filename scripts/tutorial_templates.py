"""HTML for the tutorials: a GUI walkthrough page, a notebook page, the index and the long guide.

The layout follows spaCR's tutorials on purpose -- a settings panel whose every field explains
itself on hover, a figure panel stepped with a slider, a console whose every line explains itself,
and a button bar -- so someone who has used one reads the other without learning a new format.
Everything shown is produced by `build_tutorials.py` from the running application; nothing here is
content, only the frame it is put in.
"""
from __future__ import annotations

import html
import json

STYLE = """
@font-face { font-family: 'Open Sans'; src: local('Open Sans'); }
:root { --bg:#1e1e1e; --panel:#2a2a2a; --card:#333; --ink:#eaeaea; --muted:#a8a8a8;
        --accent:#4aa3ff; --pass:#3fb37f; --fail:#e0625a; --line:#444; }
html, body { margin:0; padding:0; background:var(--bg); color:var(--ink);
             font-family:'Open Sans', 'Segoe UI', Roboto, sans-serif; }
a { color: var(--accent); }
.topbar { display:flex; align-items:center; gap:16px; padding:10px 20px; background:#161616;
          border-bottom:1px solid var(--line); }
.topbar h1 { font-size:18px; margin:0; font-weight:600; }
.topbar .nav { margin-left:auto; display:flex; gap:14px; font-size:14px; }
.container { display:flex; height:calc(100vh - 50px); }
.left-panel { width:30%; min-width:300px; background:var(--panel); padding:18px; box-sizing:border-box;
              border-right:1px solid var(--line); overflow-y:auto; }
.left-panel h2, .right-panel h2 { font-size:16px; margin:0 0 10px 0; }
.intro { font-size:13px; color:var(--muted); line-height:1.5; margin-bottom:14px; }
.field { display:flex; justify-content:space-between; gap:10px; padding:7px 9px; margin:4px 0;
         background:#232323; border-radius:6px; font-size:13px; cursor:help; }
.field:hover { outline:1px solid var(--accent); }
.field .label { color:var(--muted); }
.field .value { font-family:monospace; text-align:right; word-break:break-all; max-width:55%; }
.right-panel { flex:1; display:flex; flex-direction:column; padding:18px; box-sizing:border-box;
               min-width:0; }
.figure { flex:1.4; background:var(--card); border-radius:12px; position:relative; display:flex;
          flex-direction:column; align-items:center; justify-content:center; min-height:0;
          padding:10px; box-sizing:border-box; }
.figure img { max-width:100%; max-height:calc(100% - 50px); object-fit:contain; border-radius:6px; }
.caption { font-size:13px; color:var(--muted); margin-top:8px; text-align:center; max-width:90%; }
.slider-row { display:flex; align-items:center; gap:12px; justify-content:center; height:44px; }
.slider-row input { width:60%; }
.bottom { flex:1; display:flex; gap:14px; min-height:0; margin-top:6px; }
.console { flex:1; background:#222; border-radius:10px; display:flex; flex-direction:column;
           min-height:0; }
.console .label { padding:6px 10px; font-weight:600; font-size:14px; background:#1b1b1b;
                  border-radius:10px 10px 0 0; }
.console .out { flex:1; overflow-y:auto; padding:8px 10px; font-family:monospace; font-size:12.5px; }
.console .line { padding:2px 4px; border-radius:4px; cursor:help; white-space:pre-wrap; }
.console .line:hover { background:#2f3b48; }
.buttons { width:30%; min-width:240px; display:flex; flex-direction:column; gap:8px; }
.buttons button { background:#3a3a3a; color:var(--ink); border:1px solid #555; border-radius:6px;
                  padding:8px; font-size:13px; cursor:pointer; text-align:left; }
.buttons button:hover { border-color:var(--accent); }
.buttons .desc { font-size:12.5px; color:var(--muted); min-height:60px; line-height:1.45; }
.help-box { position:fixed; max-width:420px; background:#101820; border:1px solid var(--accent);
            color:var(--ink); padding:9px 11px; border-radius:6px; font-size:12.5px; line-height:1.45;
            display:none; z-index:10; pointer-events:none; }
.verdict-pass { color:var(--pass); font-weight:700; } .verdict-fail { color:var(--fail); font-weight:700; }
/* notebook */
.nb { max-width:1100px; margin:0 auto; padding:20px 24px 80px 24px; }
.nb h1 { font-size:26px; } .nb h2 { font-size:20px; margin-top:34px; border-bottom:1px solid var(--line);
         padding-bottom:4px; }
.nb p, .nb li { line-height:1.6; font-size:15px; color:#dcdcdc; }
.cell { margin:14px 0; }
.cell .in, .cell .outp { display:flex; gap:10px; }
.cell .prompt { width:62px; flex:none; text-align:right; font-family:monospace; font-size:12.5px;
                color:#7aa7d9; padding-top:9px; }
.cell .outp .prompt { color:#d98b7a; }
pre.code { flex:1; margin:0; background:#262a30; border:1px solid #3a3f47; border-radius:6px;
           padding:9px 11px; overflow-x:auto; font-size:13px; line-height:1.5; }
pre.code .ln { display:block; cursor:help; border-radius:3px; }
pre.code .ln:hover { background:#324257; }
.outbox { flex:1; background:#1a1a1a; border-left:3px solid #555; padding:8px 11px; overflow-x:auto;
          font-size:13px; }
.outbox pre { margin:0; white-space:pre-wrap; font-size:12.5px; }
.outbox table { border-collapse:collapse; font-size:12.5px; }
.outbox th, .outbox td { border:1px solid #444; padding:3px 7px; text-align:left; }
.outbox th { background:#2c2c2c; }
.outbox img { max-width:100%; border-radius:6px; }
.note { background:#23303d; border-left:3px solid var(--accent); padding:8px 12px; border-radius:4px;
        font-size:14px; line-height:1.55; }
/* index and guide */
.page { max-width:1100px; margin:0 auto; padding:24px 24px 90px 24px; }
.page h1 { font-size:28px; } .page h2 { font-size:21px; margin-top:38px;
           border-bottom:1px solid var(--line); padding-bottom:5px; }
.page h3 { font-size:17px; margin-top:24px; }
.page p, .page li { line-height:1.65; font-size:15px; color:#dcdcdc; }
.page img.shot { max-width:100%; border:1px solid #444; border-radius:8px; margin:10px 0; }
.page table { border-collapse:collapse; font-size:13.5px; margin:10px 0; width:100%; }
.page th, .page td { border:1px solid #444; padding:5px 8px; text-align:left; vertical-align:top; }
.page th { background:#2a2a2a; }
.cards { display:grid; grid-template-columns:repeat(auto-fill, minmax(300px, 1fr)); gap:14px; }
.card { background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:14px 16px; }
.card h3 { margin:0 0 6px 0; font-size:16px; }
.card p { font-size:13.5px; color:var(--muted); margin:0 0 10px 0; }
.card a { margin-right:14px; font-size:14px; }
.toc { columns:2; font-size:14px; }
"""

SCRIPT = """
const help = document.getElementById('help-overlay');
function showHelp(ev, text) { help.innerHTML = text; help.style.display = 'block';
  const x = Math.min(ev.clientX + 16, window.innerWidth - 440);
  const y = Math.min(ev.clientY + 16, window.innerHeight - help.offsetHeight - 10);
  help.style.left = x + 'px'; help.style.top = y + 'px'; }
function hideHelp() { help.style.display = 'none'; }
document.querySelectorAll('[data-tip]').forEach(el => {
  el.addEventListener('mousemove', ev => showHelp(ev, el.getAttribute('data-tip')));
  el.addEventListener('mouseleave', hideHelp); });
"""


def _e(text) -> str:
    return html.escape(str(text))


def page(title: str, body: str, nav: str = "", script: str = "") -> str:
    """One complete HTML page in the tutorial style."""
    return (f"<!DOCTYPE html>\n<html lang='en'><head><meta charset='UTF-8'>"
            f"<meta name='viewport' content='width=device-width, initial-scale=1'>"
            f"<title>{_e(title)}</title><style>{STYLE}</style></head><body>"
            f"<div class='topbar'><h1>{_e(title)}</h1><div class='nav'>{nav}</div></div>"
            f"{body}<div id='help-overlay' class='help-box'></div>"
            f"<script>{SCRIPT}{script}</script></body></html>")


def nav_links(index="index.html", guide="complete_guide.html") -> str:
    return (f"<a href='{index}'>All tutorials</a><a href='{guide}'>Complete guide</a>"
            f"<a href='../strategies.md'>Strategy catalogue</a>")


def gui_page(title: str, intro: str, fields: list, slides: list, console: list,
             buttons: list) -> str:
    """spaCR-style GUI walkthrough.

    fields  -- [(label, value, tip)] shown in the settings panel
    slides  -- [(image path, caption)] stepped with the slider
    console -- [(line, tip)] shown in the console
    buttons -- [(label, description)]
    """
    f_html = "".join(f"<div class='field' data-tip='{_e(tip)}'><span class='label'>{_e(l)}"
                     f"</span><span class='value'>{_e(v)}</span></div>" for l, v, tip in fields)
    c_html = "".join(f"<div class='line' data-tip='{_e(tip)}'>{_e(line)}</div>"
                     for line, tip in console)
    b_html = "".join(f"<button onmouseenter='describe({i})' onclick='describe({i})'>{_e(l)}</button>"
                     for i, (l, _d) in enumerate(buttons))
    body = (f"<div class='container'><div class='left-panel'><h2>Settings</h2>"
            f"<div class='intro'>{intro}</div>{f_html}</div>"
            f"<div class='right-panel'><div class='figure'><img id='slide' src='{_e(slides[0][0])}'>"
            f"<div class='caption' id='caption'>{_e(slides[0][1])}</div></div>"
            f"<div class='slider-row'><span>1</span><input id='slider' type='range' min='0' "
            f"max='{len(slides) - 1}' value='0' oninput='go(this.value)'><span>{len(slides)}</span>"
            f"</div><div class='bottom'><div class='console'><div class='label'>Console</div>"
            f"<div class='out'>{c_html}</div></div><div class='buttons'>{b_html}"
            f"<div class='desc' id='desc'>Hover or click a button to see what it does. Hover any "
            f"setting or console line for an explanation; move the slider to step through the "
            f"run.</div></div></div></div></div>")
    script = (f"const slides = {json.dumps(slides)}; const descs = {json.dumps([d for _l, d in buttons])};"
              "function go(i) { document.getElementById('slide').src = slides[i][0];"
              " document.getElementById('caption').textContent = slides[i][1]; }"
              "function describe(i) { document.getElementById('desc').textContent = descs[i]; }")
    return page(title, body, nav_links(), script)


def code_block(lines: list) -> str:
    """[(code line, explanation)] -> a code cell whose lines explain themselves on hover."""
    return "<pre class='code'>" + "".join(
        f"<span class='ln' data-tip='{_e(tip)}'>{_e(code) if code else '&#8203;'}</span>"
        for code, tip in lines) + "</pre>"


def notebook_page(title: str, cells: list) -> str:
    """cells: dicts of kind 'md' (html) or 'code' (lines, output html, count)."""
    parts = [f"<div class='nb'>"]
    for c in cells:
        if c["kind"] == "md":
            parts.append(c["html"])
            continue
        parts.append(f"<div class='cell'><div class='in'><div class='prompt'>In [{c['count']}]:"
                     f"</div>{code_block(c['lines'])}</div>")
        if c.get("output"):
            parts.append(f"<div class='outp'><div class='prompt'>Out[{c['count']}]:</div>"
                         f"<div class='outbox'>{c['output']}</div></div>")
        parts.append("</div>")
    parts.append("</div>")
    return page(title, "".join(parts), nav_links())
