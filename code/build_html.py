#!/usr/bin/env python3
"""Build survey.html (a self-contained web page) from survey.tex via pandoc.

Post-processing of pandoc's output:
  * theorem-like environments are renumbered section-wise (2.3, 2.4, ...) as in
    the PDF, and \\ref links are updated accordingly;
  * \\cite spans (left empty by pandoc) become [n] links to the bibliography,
    numbered in the order of the \\bibitem entries;
  * the page is wrapped in our own head/style (no <html>/<body> wrapper:
    the artifact host supplies those).
"""
import re, subprocess, sys, os, html

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, "survey.tex")
OUT = os.path.join(ROOT, "survey.html")

body = subprocess.run(
    ["pandoc", TEX, "-f", "latex", "-t", "html5", "--mathjax", "--toc",
     "--toc-depth=2", "--shift-heading-level-by=1", "--wrap=none"],
    capture_output=True, text=True, check=True, cwd=ROOT).stdout

tex = open(TEX, encoding="utf-8").read()

# ---------------------------------------------------------------- equations
# Pandoc wraps LaTeX equation environments in \[...\], which leaves invalid
# nested display delimiters and prints \eqref targets as raw label names.
equation_labels = re.findall(
    r"\\begin\{equation\}\s*\\label\{([^}]*)\}", tex, flags=re.S
)
equation_num = {label: i + 1 for i, label in enumerate(equation_labels)}

def equation_repl(m):
    label, content = m.group(1), m.group(2).strip()
    number = equation_num[label]
    return (
        f'<span id="{label}" class="math display">'
        f'\\[{content}\\tag{{{number}}}\\]</span>'
    )

body = re.sub(
    r'<span class="math display">\\\[\\begin\{equation\}\s*'
    r'\\label\{([^}]*)\}\s*(.*?)\\end\{equation\}\\\]</span>',
    equation_repl,
    body,
    flags=re.S,
)

def eqref_repl(m):
    label = m.group(1)
    number = equation_num.get(label)
    if number is None:
        return m.group(0)
    return f'<a href="#{label}" class="ref">({number})</a>'

body = re.sub(
    r'<a href="#([^"]*)" data-reference-type="eqref" '
    r'data-reference="[^"]*">.*?</a>',
    eqref_repl,
    body,
)

# ---------------------------------------------------------------- citations
bibkeys = re.findall(r"\\bibitem\{([^}]*)\}", tex)
bibnum = {k: i + 1 for i, k in enumerate(bibkeys)}

def cite_repl(m):
    keys = m.group(1).split()
    parts = []
    for k in keys:
        n = bibnum.get(k)
        if n is None:
            parts.append(html.escape(k))
        else:
            parts.append(f'<a class="cite" href="#bib-{k}">{n}</a>')
    return "[" + ", ".join(parts) + "]"

body = re.sub(r'<span class="citation" data-cites="([^"]*)">\s*</span>', cite_repl, body)

# bibliography: pandoc emits <div class="thebibliography"><p><span>99</span></p><p>...</p>...
def bib_repl(m):
    inner = m.group(1)
    paras = re.findall(r"<p>(.*?)</p>", inner, flags=re.S)
    paras = [p for p in paras if not re.fullmatch(r"\s*<span>99</span>\s*", p)]
    assert len(paras) == len(bibkeys), (len(paras), len(bibkeys))
    items = "".join(f'<li id="bib-{k}"><span class="bibnum">[{i+1}]</span><span class="bibtext">{p}</span></li>'
                    for i, (k, p) in enumerate(zip(bibkeys, paras)))
    return f'<h2 id="references">References</h2><ol class="bibliography">{items}</ol>'

body, n = re.subn(r'<div class="thebibliography">(.*?)</div>\s*$', bib_repl, body, flags=re.S)
assert n == 1

# ---------------------------------------------------------------- theorem numbering
KINDS = ["theorem", "proposition", "lemma", "corollary", "conjecture", "question",
         "definition", "example", "remark", "observation"]
NAMES = {k: k.capitalize() for k in KINDS}

# walk the document in order: section headings (h2 after shift) and theorem divs
pattern = re.compile(r'<h2 id="([^"]*)"(?: class="([^"]*)")?>|<div id="([^"]*)" class="(' + "|".join(KINDS) + r')">\s*<p><strong>([A-Za-z]+) (\d+)</strong>')
sec = 0
cnt = 0
newnum = {}          # old pandoc number (as string) -> new label
labelnum = {}        # div id -> new label
def walk(m):
    global sec, cnt
    if m.group(1) is not None:                      # h2
        if m.group(2) and "unnumbered" in m.group(2):
            return m.group(0)
        sec += 1; cnt = 0
        return m.group(0)
    div_id, kind, name, old = m.group(3), m.group(4), m.group(5), m.group(6)
    cnt += 1
    new = f"{sec}.{cnt}"
    newnum[old] = new
    labelnum[div_id] = new
    return f'<div id="{div_id}" class="thm {kind}">\n<p><span class="thm-head">{name} {new}</span>'
body = pattern.sub(walk, body)

# references to theorems: <a href="#label" data-reference-type="ref" data-reference="label">N</a>
def ref_repl(m):
    label, txt = m.group(1), m.group(2)
    if label in labelnum:
        return f'<a href="#{label}" class="ref">{labelnum[label]}</a>'
    return f'<a href="#{label}" class="ref">{txt}</a>'
body = re.sub(r'<a href="#([^"]*)" data-reference-type="ref" data-reference="[^"]*">([^<]*)</a>', ref_repl, body)

# theorem heads with optional names "(Jaeger [5])" -> keep, then the statement is <em>...</em>
body = body.replace('<div class="proof">\n<p><em>Proof.</em>', '<div class="proof">\n<p><span class="proof-head">Proof.</span>')
body = body.replace("◻", '<span class="qed">&#9633;</span>')

# section numbering for h2/h3 (mirrors LaTeX): h2 numbered unless class unnumbered; h3 numbered within
sec = 0; sub = 0
def head_repl(m):
    global sec, sub
    level, hid, cls, txt = m.group(1), m.group(2), m.group(3) or "", m.group(4)
    if "unnumbered" in cls or hid in ("references",):
        return m.group(0)
    if level == "2":
        sec += 1; sub = 0
        return f'<h2 id="{hid}"><span class="secnum">{sec}</span>{txt}</h2>'
    sub += 1
    return f'<h3 id="{hid}"><span class="secnum">{sec}.{sub}</span>{txt}</h3>'
body = re.sub(r'<h([23]) id="([^"]*)"(?: class="([^"]*)")?>(.*?)</h\1>', head_repl, body, flags=re.S)

# TOC: pandoc puts <nav id="TOC" role="doc-toc"> ... ; add numbers via CSS counters instead.
# title block: pandoc (without -s) does not emit the title; we add our own.

head = r'''<title>The Alspach–Zhang Conjecture</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;0,8..60,700;1,8..60,400;1,8..60,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400&display=swap">
<style>
:root{
  --bg:#f6f7f5; --bg-raised:#ffffff; --ink:#1c232b; --ink-2:#4b5760; --ink-3:#7a8790;
  --rule:#d5dbd8; --rule-strong:#9aa6a1; --accent:#0e6b67; --accent-ink:#0b4d4a; --accent-soft:#e3efed;
  --thm-bg:#f0f4f2; --code-bg:#eef1ef;
  --serif:"Source Serif 4","Iowan Old Style","Palatino Linotype",Georgia,serif;
  --sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
  --mono:"IBM Plex Mono",Menlo,Consolas,monospace;
  color-scheme:light;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#151a1e; --bg-raised:#1c2227; --ink:#e6eaec; --ink-2:#b4bec4; --ink-3:#87939b;
    --rule:#2e373d; --rule-strong:#55626a; --accent:#63c2ba; --accent-ink:#8fd6d0; --accent-soft:#1d2f2e;
    --thm-bg:#1b2327; --code-bg:#1f272c; color-scheme:dark;
  }
}
:root[data-theme="dark"]{
  --bg:#151a1e; --bg-raised:#1c2227; --ink:#e6eaec; --ink-2:#b4bec4; --ink-3:#87939b;
  --rule:#2e373d; --rule-strong:#55626a; --accent:#63c2ba; --accent-ink:#8fd6d0; --accent-soft:#1d2f2e;
  --thm-bg:#1b2327; --code-bg:#1f272c; color-scheme:dark;
}
html{background:var(--bg)}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--serif);font-size:17px;line-height:1.55;
  -webkit-font-smoothing:antialiased}
.page{max-width:46rem;margin:0 auto;padding:3.5rem 1.25rem 5rem}
a{color:var(--accent);text-decoration:none}
a:hover{text-decoration:underline}
a:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:2px}
header.masthead{border-bottom:1px solid var(--rule-strong);padding-bottom:1.5rem;margin-bottom:2rem}
.eyebrow{font-family:var(--sans);font-size:.72rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);margin:0 0 .9rem}
h1.doc-title{font-size:2.1rem;line-height:1.15;font-weight:700;margin:0 0 .6rem;text-wrap:balance;letter-spacing:-.01em}
.subtitle{font-size:1.15rem;color:var(--ink-2);margin:0 0 1rem;font-style:italic}
.meta{font-family:var(--sans);font-size:.85rem;color:var(--ink-3);margin:0}
.abstract{background:var(--bg-raised);border:1px solid var(--rule);padding:1.1rem 1.4rem;margin:0 0 2rem;font-size:.97rem;line-height:1.55}
.abstract .label{font-family:var(--sans);font-size:.72rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--ink-3);display:block;margin-bottom:.4rem}
.abstract p{margin:0}
nav#TOC{font-family:var(--sans);font-size:.9rem;margin:0 0 2.5rem;padding:1rem 0 1.25rem;border-top:1px solid var(--rule);border-bottom:1px solid var(--rule)}
nav#TOC::before{content:"Contents";display:block;font-size:.72rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--ink-3);margin-bottom:.6rem}
nav#TOC ul{list-style:none;margin:0;padding:0}
nav#TOC > ul{columns:2;column-gap:2.5rem}
nav#TOC > ul > li{break-inside:avoid;margin:0 0 .5rem;font-weight:500}
nav#TOC ul ul{margin:.2rem 0 0 0;font-weight:400}
nav#TOC ul ul li{margin:0 0 .15rem;color:var(--ink-2)}
nav#TOC a{color:inherit;display:grid;grid-template-columns:2.4em 1fr;gap:.2rem}
.tocnum{color:var(--ink-3);font-variant-numeric:tabular-nums}
h2,h3{font-weight:600;line-height:1.25;text-wrap:balance;margin:2.6rem 0 .9rem}
h2{font-size:1.5rem;padding-top:1.2rem;border-top:1px solid var(--rule)}
h3{font-size:1.15rem;margin-top:2rem}
.secnum{display:inline-block;min-width:2.2em;color:var(--accent);font-family:var(--sans);font-weight:500;font-variant-numeric:tabular-nums}
h3 .secnum{min-width:2.6em}
p{margin:0 0 1rem}
ul,ol{padding-left:1.5rem;margin:0 0 1rem}
li{margin:0 0 .35rem}
strong{font-weight:600}
.thm{background:var(--thm-bg);border-left:3px solid var(--accent);padding:.85rem 1.1rem .5rem;margin:1.2rem 0 1.2rem;border-radius:0 3px 3px 0}
.thm p{margin-bottom:.6rem}
.thm-head{font-family:var(--sans);font-weight:600;font-size:.92rem;color:var(--accent-ink);letter-spacing:.01em}
.thm.conjecture{border-left-color:var(--rule-strong);background:var(--bg-raised);border:1px solid var(--rule);border-left:3px solid var(--rule-strong)}
.proof{margin:0 0 1.2rem;padding-left:1.1rem;border-left:1px solid var(--rule)}
.proof-head{font-style:italic;color:var(--ink-2)}
.qed{float:right;color:var(--ink-3)}
blockquote{margin:1rem 0 1rem 1.2rem;padding-left:1rem;border-left:2px solid var(--rule-strong);color:var(--ink-2);font-size:.97rem}
table{border-collapse:collapse;font-family:var(--sans);font-size:.85rem;line-height:1.4;margin:0 auto 1.2rem;font-variant-numeric:tabular-nums}
.tablewrap{overflow-x:auto;margin:1.4rem 0 1.6rem}
table caption{caption-side:bottom;text-align:left;font-size:.82rem;color:var(--ink-2);padding:.6rem 0 0;font-family:var(--serif)}
thead th{text-align:left;font-weight:600;border-top:1.5px solid var(--rule-strong);border-bottom:1px solid var(--rule-strong);padding:.45rem .7rem;vertical-align:bottom}
tbody td{padding:.35rem .7rem;vertical-align:top;border-bottom:1px solid var(--rule)}
tbody tr:last-child td{border-bottom:1.5px solid var(--rule-strong)}
td.right,th.right{text-align:right}
code{font-family:var(--mono);font-size:.85em;background:var(--code-bg);padding:.05em .3em;border-radius:3px}
pre{overflow-x:auto;background:var(--code-bg);padding:.8rem 1rem;border-radius:4px;font-size:.85rem}
.bibliography{list-style:none;padding:0;margin:0;font-size:.92rem}
.bibliography li{display:grid;grid-template-columns:3em 1fr;gap:.4rem;margin:0 0 .55rem}
.bibnum{font-family:var(--sans);color:var(--ink-3);font-variant-numeric:tabular-nums}
.bibliography li:target,.thm:target{outline:2px solid var(--accent);outline-offset:3px}
a.cite,a.ref{font-variant-numeric:tabular-nums}
mjx-container{overflow-x:auto;overflow-y:hidden;max-width:100%}
mjx-container[display="true"]{margin:.8rem 0 !important}
.colophon{margin-top:3rem;padding-top:1rem;border-top:1px solid var(--rule);font-family:var(--sans);font-size:.8rem;color:var(--ink-3)}
@media (max-width:640px){body{font-size:16px}.page{padding:2rem 1rem 4rem}nav#TOC > ul{columns:1}h1.doc-title{font-size:1.7rem}}
@media (prefers-reduced-motion: reduce){*{scroll-behavior:auto}}
</style>
<script>
window.MathJax = {tex:{inlineMath:[['\\(','\\)']],displayMath:[['\\[','\\]']],tags:'ams'},
  svg:{fontCache:'global'}, options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.js" async></script>
'''

masthead = '''<div class="page">
<header class="masthead">
<p class="eyebrow">Survey · Algebraic graph theory</p>
<h1 class="doc-title">The Alspach–Zhang conjecture: every cubic Cayley graph is 3-edge-colourable</h1>
<p class="subtitle">A survey of what is known, with self-contained proofs of the elementary reductions and a computer verification beyond the census range</p>
<p class="meta">September 2026 · LaTeX source, code and data: <code>survey.tex</code>, <code>code/</code> in the accompanying repository</p>
</header>
'''

# abstract: pandoc treats it as metadata, so convert it separately and insert it before the TOC
abs_tex = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, flags=re.S).group(1)
preamble = tex[:tex.index("\\begin{document}")]
abs_html = subprocess.run(["pandoc", "-f", "latex", "-t", "html5", "--mathjax", "--wrap=none"],
                          input=preamble + "\\begin{document}" + abs_tex + "\\end{document}",
                          capture_output=True, text=True, check=True).stdout
# table of contents (pandoc emits one only in standalone mode): build it from the numbered headings
toc_items = []
for level, hid, num, txt in re.findall(r'<h([23]) id="([^"]*)"><span class="secnum">([^<]*)</span>(.*?)</h\1>', body, flags=re.S):
    toc_items.append((int(level), hid, num, re.sub(r"<[^>]+>", "", txt)))
toc = ['<nav id="TOC" aria-label="Contents"><ul>']
open_sub = False
for level, hid, num, txt in toc_items:
    if level == 2:
        if open_sub:
            toc.append("</ul></li>"); open_sub = False
        else:
            if len(toc) > 1: toc.append("</li>")
        toc.append(f'<li><a href="#{hid}"><span class="tocnum">{num}</span>{txt}</a>')
    else:
        if not open_sub:
            toc.append("<ul>"); open_sub = True
        toc.append(f'<li><a href="#{hid}"><span class="tocnum">{num}</span>{txt}</a></li>')
if open_sub: toc.append("</ul>")
toc.append("</li></ul></nav>")
first_h2 = body.index('<h2 ')
body = body[:first_h2] + '<div class="abstract"><span class="label">Abstract</span>' + abs_html + '</div>\n' + "".join(toc) + "\n" + body[first_h2:]
# wrap tables for horizontal scrolling
body = re.sub(r'(<table>.*?</table>)', r'<div class="tablewrap">\1</div>', body, flags=re.S)
# right-align numeric columns in tables: mark cells that are pure numbers
body = re.sub(r'<td>(\d[\d,]*)</td>', r'<td class="right">\1</td>', body)

colophon = '''<p class="colophon">Prepared with the help of AI assistants (Claude and OpenAI Codex). Bibliographic details were checked against the cited sources where accessible; the computations are reproducible from the repository code.</p>
</div>
'''
open(OUT, "w", encoding="utf-8").write(head + masthead + body + colophon)
print("wrote", OUT, len(head + masthead + body + colophon), "bytes")

# standalone version for GitHub Pages (docs/survey.html)
DOCS = os.path.join(ROOT, "docs", "survey.html")
os.makedirs(os.path.dirname(DOCS), exist_ok=True)
standalone = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
              '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
              '<meta name="description" content="Survey of the Alspach-Zhang conjecture: every cubic Cayley graph is 3-edge-colourable.">\n'
              + head + '<style>img{max-width:100%}</style>\n</head>\n<body>\n'
              + '<p class="backlink"><a href="./">&larr; Project page</a></p>\n'
              + masthead + body + colophon + '</body>\n</html>\n')
standalone = standalone.replace('<p class="backlink">', '<p class="backlink" style="max-width:46rem;margin:1.2rem auto -2.5rem;padding:0 1.25rem;font-family:var(--sans);font-size:.85rem">')
open(DOCS, "w", encoding="utf-8").write(standalone)
print("wrote", DOCS)
