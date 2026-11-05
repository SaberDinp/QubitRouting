"""Convenience entry point for IDEs and invocation from outside the repository."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from luteberget.__main__ import main


if __name__ == "__main__":
    raise SystemExit(main())
