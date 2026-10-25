import glob, re
from pypdf import PdfReader
for f in sorted(glob.glob("*.pdf")):
    try:
        t = PdfReader(f).pages[0].extract_text() or ""
    except Exception as e:
        t = "ERR %s" % e
    t = re.sub(r"\s+", " ", t)[:120]
    print("%-38s | %s" % (f, t))
