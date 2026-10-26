"""Rebuild the vendored Clarabel WebAssembly module: vendor/clarabel/.

Compiles wasm/clarabel (a minimal binding of Clarabel.rs, pinned in Cargo.toml and Cargo.lock) in a
Docker container, so no local Rust toolchain is needed, and writes:
  vendor/clarabel/clarabel_bb.js        wasm-bindgen glue (target "no-modules": defines `wasm_bindgen`)
  vendor/clarabel/clarabel_bb_bg.wasm   the module
Only needed when wasm/clarabel changes; the build (scripts/build.py) uses the vendored files.

Usage:  python scripts/build_clarabel.py
Needs:  Docker (Docker Desktop on Windows) running.
"""
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
CRATE = ROOT / "wasm" / "clarabel"
OUT = ROOT / "vendor" / "clarabel"
IMAGE = "bb-clarabel-build"


def run(args):
    print(" ".join(str(a) for a in args))
    subprocess.run([str(a) for a in args], check=True)


def main() -> int:
    if shutil.which("docker") is None:
        print("Docker is not installed or not on PATH.", file=sys.stderr)
        return 1
    run(["docker", "build", "-t", IMAGE, str(CRATE)])
    script = (
        "cargo build --release --target wasm32-unknown-unknown && "
        "wasm-bindgen --target no-modules --no-typescript --out-dir /src/pkg "
        "target/wasm32-unknown-unknown/release/clarabel_bb.wasm"
    )
    run(["docker", "run", "--rm", "-v", f"{CRATE}:/src", "-v", "bb-clarabel-cargo:/usr/local/cargo/registry",
         IMAGE, "sh", "-c", script])
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ("clarabel_bb.js", "clarabel_bb_bg.wasm"):
        shutil.copyfile(CRATE / "pkg" / name, OUT / name)
        print(f"Wrote {(OUT / name).relative_to(ROOT)} ({(OUT / name).stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
