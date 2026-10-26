"""Build the single-file app: dist/bearing-bounds.html.

Inlines the HiGHS solver and the Clarabel conic solver (each JavaScript + WebAssembly as base64) and
src/engine.js into src/app.html. Works on Windows, macOS and Linux (pathlib only).
Clarabel is vendored in vendor/clarabel (rebuild it with scripts/build_clarabel.py).

Usage:  python scripts/build.py
Needs:  `npm install` first (provides node_modules/highs).
"""
import base64
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
HIGHS = ROOT / "node_modules" / "highs" / "build"
CLARABEL = ROOT / "vendor" / "clarabel"
DIST = ROOT / "dist"


def main() -> int:
    if not (HIGHS / "highs.wasm").exists():
        print("node_modules/highs is missing. Run `npm install` first.", file=sys.stderr)
        return 1
    if not (CLARABEL / "clarabel_bb_bg.wasm").exists():
        print("vendor/clarabel is missing. Run `python scripts/build_clarabel.py`.", file=sys.stderr)
        return 1
    page = (SRC / "app.html").read_text(encoding="utf-8")
    for marker in ("/*__HIGHS__*/", "/*__ENGINE__*/", "__WASM_B64__", "/*__CLARABEL__*/", "__CLARABEL_B64__"):
        if marker not in page:
            print(f"Template marker {marker} not found in src/app.html", file=sys.stderr)
            return 1
    page = page.replace("/*__HIGHS__*/", (HIGHS / "highs.js").read_text(encoding="utf-8"))
    page = page.replace("/*__ENGINE__*/", (SRC / "engine.js").read_text(encoding="utf-8"))
    page = page.replace("__WASM_B64__", base64.b64encode((HIGHS / "highs.wasm").read_bytes()).decode("ascii"))
    page = page.replace("/*__CLARABEL__*/", (CLARABEL / "clarabel_bb.js").read_text(encoding="utf-8"))
    page = page.replace("__CLARABEL_B64__", base64.b64encode((CLARABEL / "clarabel_bb_bg.wasm").read_bytes()).decode("ascii"))
    DIST.mkdir(exist_ok=True)
    out = DIST / "bearing-bounds.html"
    out.write_text(page, encoding="utf-8")
    # A standalone copy with a document skeleton, for opening directly in a browser
    standalone = DIST / "bearing-bounds-standalone.html"
    standalone.write_text(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1,viewport-fit=cover\">"
        "<style>body{margin:0}[hidden]{display:none!important}</style></head><body>"
        + page + "</body></html>",
        encoding="utf-8",
    )
    print(f"Wrote {out.relative_to(ROOT)} ({out.stat().st_size / 1e6:.2f} MB) and {standalone.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
