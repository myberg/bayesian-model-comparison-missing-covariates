"""Static validation for the curated publication code bundle (stdlib only)."""
from pathlib import Path
import hashlib
import re

ROOT = Path(__file__).resolve().parents[1]
BANNED = re.compile("pa" + "tch" + "|" + "re" + "vision", re.IGNORECASE)
TEXT_SUFFIXES = {".R", ".jl", ".md", ".toml", ".cff", ".example", ".txt", ".py"}

files = [p for p in ROOT.rglob("*") if p.is_file() and "__pycache__" not in p.parts]
errors = []
for p in files:
    rel = p.relative_to(ROOT).as_posix()
    if BANNED.search(rel):
        errors.append(f"Old name in path: {rel}")
    if p.suffix.lower() in {s.lower() for s in TEXT_SUFFIXES} or p.name.startswith("."):
        text = p.read_text("utf-8", errors="replace")
        if BANNED.search(text):
            errors.append(f"Old name in contents: {rel}")
        for target in re.findall(r'(?m)^\s*source\(["\'](R/[^"\']+\.R)["\']\)', text):
            if not (ROOT / target).is_file():
                errors.append(f"Missing R source referenced by {rel}: {target}")

source_hashes = {}
for p in files:
    if p.suffix.lower() in {".r", ".jl"}:
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        if digest in source_hashes:
            errors.append(f"Duplicate source content: {source_hashes[digest]} and {p.relative_to(ROOT)}")
        source_hashes[digest] = p.relative_to(ROOT)

if not (ROOT / "julia/src/BayesREProbit.jl").is_file():
    errors.append("Shared Julia module missing")
if (ROOT / "neps_sc6_empirical/src/BayesREProbit.jl").exists():
    errors.append("Duplicate Julia module remains")

if errors:
    print("\n".join(errors))
    raise SystemExit(1)
print(f"Static checks passed: {len(files)} files examined")
