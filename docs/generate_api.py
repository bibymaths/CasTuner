import os
from pathlib import Path

# Root of the repository
ROOT = Path(__file__).resolve().parents[1]

SCRIPTS_DIR = ROOT / "scripts"
API_DIR = Path("docs/api")       # final Markdown will be placed here

API_DIR.mkdir(parents=True, exist_ok=True)


def python_modules(base=SCRIPTS_DIR):
    """
    Yield dotted Python module paths for all .py files in scripts/.
    """
    for root, dirs, files in os.walk(base):
        for f in files:
            if f.endswith(".py") and not f.startswith("__"):
                full = Path(root) / f
                rel = full.relative_to(ROOT)        # scripts/... path
                module = ".".join(rel.with_suffix("").parts)
                yield module


# -------------------------------
# Create API index
# -------------------------------
index = API_DIR / "index.md"
with index.open("w") as f:
    f.write("# API Reference\n\n")
    f.write("Automatically generated from all `.py` modules in `scripts/`.\n\n")
    f.write("## Modules\n\n")
    for module in sorted(python_modules()):
        name = module.split(".")[-1]
        f.write(f"- [{name}](./{module}.md)\n")


# -------------------------------
# Create one page per module
# -------------------------------
for module in python_modules():
    module_path = API_DIR / f"{module}.md"
    module_path.parent.mkdir(parents=True, exist_ok=True)

    with module_path.open("w") as f:
        f.write(f"# `{module}`\n\n")
        f.write(f"::: {module}\n")