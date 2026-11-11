"""Synthetic demo data for trying the app. Not a real biological sequence.

The template is random DNA with three deliberate edits, so every behaviour can be seen:
  * pair A  -> 500 bp product, plus an off-target site (a copy of the A forward
               primer site with 2 mismatches at its 5' end) that only binds at
               low annealing temperatures
  * pair B  -> 1000 bp product
  * pair D  -> 320 bp product, primers designed with complementary 3' ends
               (5 bp), so they form primer dimers
"""

import random

from .sequence_utils import reverse_complement

SWAP = {"A": "C", "C": "A", "G": "T", "T": "G"}


def make_demo():
    """Return (templates, primer_pairs) built from a fixed random sequence."""
    rng = random.Random(1)
    seq = [rng.choice("ACGT") for _ in range(1500)]

    # Off-target site for pair A: copy 100-120 to 300-320 with 2 mismatches near the 5' end
    off_target = seq[100:120]
    off_target[1] = SWAP[off_target[1]]
    off_target[4] = SWAP[off_target[4]]
    seq[300:320] = off_target

    # Dimer-prone pair D: the first 5 bases of the D reverse site copy the
    # last 5 bases of the D forward primer, so the two 3' ends pair
    seq[1000:1005] = seq[715:720]

    template = "".join(seq)
    templates = {"Demo_synthetic_1500bp": template}
    pairs = {
        "Pair": ["A", "B", "D"],
        "Forward": [template[100:120], template[200:220], template[700:720]],
        "Reverse": [reverse_complement(template[580:600]),
                    reverse_complement(template[1180:1200]),
                    reverse_complement(template[1000:1020])],
    }
    return templates, pairs