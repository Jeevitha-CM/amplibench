"""Primer-dimer check: can the 3' end of one primer pair with another primer?"""

from dataclasses import dataclass
from itertools import combinations_with_replacement

from .conditions import ReactionConditions
from .sequence_utils import reverse_complement


@dataclass
class Dimer:
    primer_a: str     # the primer whose 3' end pairs and gets extended
    primer_b: str     # the primer it pairs with (same name = self-dimer)
    overlap: int      # number of 3'-end bases that pair perfectly
    gc: int           # G/C bases in that overlap (stronger pairing)
    size: int         # approximate length of the dimer product (bp)


def three_prime_pairing(a: str, b: str, min_run: int) -> tuple[int, int] | None:
    """Longest 3' end of primer a that pairs perfectly with primer b.

    Returns (overlap length, position in b) or None.
    """
    for k in range(min(len(a), len(b)), min_run - 1, -1):
        target = reverse_complement(a[-k:])
        pos = b.find(target)
        if pos != -1:
            return k, pos
    return None


def find_dimers(primers: dict[str, str], conditions: ReactionConditions) -> list[Dimer]:
    """Check every primer against every primer (including itself)."""
    # Higher primer concentration lets weaker (shorter) 3' pairings form dimers
    min_run = 3 if conditions.primer_uM >= 0.8 else 4

    dimers = []
    for (name_1, seq_1), (name_2, seq_2) in combinations_with_replacement(primers.items(), 2):
        best = None
        for (na, a), (nb, b) in [((name_1, seq_1), (name_2, seq_2)),
                                 ((name_2, seq_2), (name_1, seq_1))]:
            hit = three_prime_pairing(a, b, min_run)
            if hit and (best is None or hit[0] > best.overlap):
                k, pos = hit
                run = a[-k:]
                best = Dimer(na, nb, k, run.count("G") + run.count("C"), len(a) + pos)
        if best:
            dimers.append(best)

    return sorted(dimers, key=lambda d: d.overlap, reverse=True)


def dimer_intensity(dimer: Dimer, conditions: ReactionConditions,
                    template_present: bool, has_products: bool) -> float:
    """Illustrative band brightness (0 to 0.75) for a primer dimer."""
    strength = 0.1 * (dimer.overlap - 2) + 0.04 * dimer.gc
    concentration = min(1.5, conditions.primer_uM / 0.4)
    value = strength * concentration
    if not template_present:
        value *= 1.6      # nothing else to bind: dimers dominate (typical NTC)
    elif has_products:
        value *= 0.6      # the real product competes for primers
    return round(min(0.75, value), 2)