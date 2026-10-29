"""Build self-contained copies of the GEE scripts (no require()).

The scripts in gee/ load gee/config.js and gee/lib/common.js with require(),
which only works when the repository is cloned as an Earth Engine script
repository. This script writes copies with both modules pasted in, so each file
can be pasted straight into a new Code Editor script.

Usage
    python scripts/make_standalone.py --asset-root projects/<cloud-project>/assets/nigeria_fragmentation

Output
    gee_standalone/<script>.js   (git-ignored; regenerate after editing gee/)
"""
import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEE = ROOT / "gee"
OUT = ROOT / "gee_standalone"
REQUIRE_BLOCK = re.compile(
    r"var REPO = '[^']*';\s*\nvar cfg = require\(REPO \+ ':config'\);\s*\n"
    r"(?:var lib = require\(REPO \+ ':lib/common'\);\s*\n)?"
)


def module(path: Path, name: str) -> str:
    body = path.read_text()
    return (f"// ===== {name} (inlined) =====\n"
            f"var {name} = (function () {{\n  var exports = {{}};\n{body}\n  return exports;\n}})();\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--asset-root", default=None,
                    help="value for ASSET_ROOT in the inlined config")
    args = ap.parse_args()

    cfg = (GEE / "config.js").read_text()
    if args.asset_root:
        cfg, n = re.subn(r"var ASSET_ROOT = '[^']*';", f"var ASSET_ROOT = '{args.asset_root}';", cfg)
        assert n == 1, "ASSET_ROOT line not found in config.js"
    tmp_cfg = OUT / ".config.js"
    OUT.mkdir(exist_ok=True)
    tmp_cfg.write_text(cfg)
    cfg_block = module(tmp_cfg, "cfg")
    tmp_cfg.unlink()
    lib_block = module(GEE / "lib" / "common.js", "lib")

    for src in sorted(GEE.glob("*.js")):
        if src.name == "config.js":
            continue
        text = src.read_text()
        m = REQUIRE_BLOCK.search(text)
        if m:
            uses_lib = "lib/common" in m.group(0)
            inline = ("// Standalone version: config and lib/common are pasted below, so no require() is needed.\n"
                      "// Regenerate with scripts/make_standalone.py after editing gee/.\n"
                      + cfg_block + (lib_block if uses_lib else ""))
            text = text[:m.start()] + inline + text[m.end():]
        (OUT / src.name).write_text(text)
        print("wrote", (OUT / src.name).relative_to(ROOT))


if __name__ == "__main__":
    main()
