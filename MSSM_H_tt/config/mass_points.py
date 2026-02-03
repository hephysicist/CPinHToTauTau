# coding: utf-8
from __future__ import annotations

from pathlib import Path
from functools import lru_cache
from typing import Iterable, List
from columnflow.util import maybe_import

@lru_cache(maxsize=1)
def read_bdt_masses(path: str | Path | None = None) -> List[int]:
  """
  Load MSSM mass points from bdt_masses.yaml.
  Accepts either:
    - a dict with key 'masses'
    - a plain top-level YAML list
  Caches the result.
  """
  yaml = maybe_import("yaml")
  if yaml is None:
    raise RuntimeError("PyYAML is required to load mass points. Please `pip install pyyaml`.")

  p = Path(path) if path else Path(__file__).with_name("bdt_masses.yaml")
  if not p.exists():
    raise FileNotFoundError(f"Mass list file not found: {p}")

  data = yaml.safe_load(p.read_text())

  if isinstance(data, dict) and "masses" in data:
    masses = data["masses"]
  elif isinstance(data, list):
    masses = data
  else:
    raise ValueError("bdt_masses.yaml must be either a list or a dict with key 'masses'.")

  try:
    out = [int(m) for m in masses]
  except Exception as exc:
    raise ValueError("Failed to parse masses as integers from bdt_masses.yaml") from exc

  # optional: keep author order, but guard against accidental duplicates
  # (no sorting to preserve intended registration order)
  seen = set()
  uniq = []
  for m in out:
    if m not in seen:
      uniq.append(m)
      seen.add(m)
  return uniq
