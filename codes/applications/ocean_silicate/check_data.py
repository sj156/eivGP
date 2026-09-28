"""Check the bundled cohort, folds, and calibration masks."""
from collections import Counter
from pathlib import Path
import csv
import math

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "ocean_silicate"

def read(name):
    with (DATA / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))

def band(row):
    u = float(row["U"])
    return 1 + sum(u >= edge for edge in (50, 100, 150))

cohort = read("cohort_manifest.csv")
assert len(cohort) == len({r["cohort_id"] for r in cohort}) == 200
assert len({r["profile"] for r in cohort}) == 200
assert Counter(map(band, cohort)) == {1: 50, 2: 50, 3: 50, 4: 50}
assert all(
    300 <= float(r["G2pressure"]) <= 1500
    and 0 < float(r["U"]) < 180
    and math.isclose(float(r["log_pressure"]), math.log(float(r["G2pressure"])), abs_tol=1e-10)
    and int(r["band"]) == band(r)
    for r in cohort
)
ids = {r["cohort_id"] for r in cohort}
for fold in range(1, 6):
    train = read(f"fold{fold:02d}_train.csv")
    test = read(f"fold{fold:02d}_test.csv")
    assert len(train) == 160 and len(test) == 40
    assert {r["cohort_id"] for r in train + test} == ids
    assert not ({r["profile"] for r in train} & {r["profile"] for r in test})
    assert Counter(map(band, test)) == {1: 10, 2: 10, 3: 10, 4: 10}
    for flag, count in (("calib_o20", 8), ("calib_o50", 20)):
        assert all(
            sum(int(r[flag]) for r in train if band(r) == klass) == count
            for klass in range(1, 5)
        )
    assert all(int(r["calib_o20"]) <= int(r["calib_o50"]) for r in train)
print("Data checks passed: 200 profiles, five profile-disjoint folds, nested masks.")
