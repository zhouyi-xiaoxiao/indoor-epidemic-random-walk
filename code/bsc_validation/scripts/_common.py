"""Shared paths and environment for the validation scripts."""
import os, sys, json, time
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
os.environ.setdefault("TMPDIR", os.path.join(ROOT, ".tmp"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(ROOT, ".mplcache"))
sys.path.insert(0, os.path.join(ROOT, "src"))
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")
LOGS = os.path.join(ROOT, "logs")
SEED = 20261001

def save_json(name, obj):
    p = os.path.join(DATA, name)
    with open(p + ".part", "w") as f:
        json.dump(obj, f, indent=1, default=float)
    os.replace(p + ".part", p)
    return p

def load_json(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)

def stamp():
    return time.strftime("%Y-%m-%d %H:%M:%S %Z")
