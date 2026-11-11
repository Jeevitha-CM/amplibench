"""Find where a primer can bind on a template sequence."""

from dataclasses import dataclass

from .sequence_utils import reverse_complement


@dataclass
class BindingSite:
    primer_name: str
    strand: str        # "+" = binds as written, "-" = binds as reverse complement
    start: int         # 0-based start position on the top strand
    end: int           # end position (not included, like Python slicing)
    mismatches: int


def count_mismatches(a: str, b: str) -> int:
    """Count positions where two equal-length sequences differ."""
    return sum(1 for x, y in zip(a, b) if x != y)


def find_binding_sites(template: str, primer: str, primer_name: str,
                       max_mismatches: int = 2,
                       three_prime_exact: int = 3) -> list[BindingSite]:
    """Slide the primer along the template and return every valid binding site."""
    p_len = len(primer)
    if not 1 <= three_prime_exact <= p_len:
        raise ValueError("three_prime_exact must be between 1 and the primer length")

    primer_rc = reverse_complement(primer)
    sites = []

    for i in range(len(template) - p_len + 1):
        window = template[i:i + p_len]

        # "+" site: 3' end of the primer is at the RIGHT end of the window
        if window[-three_prime_exact:] == primer[-three_prime_exact:]:
            mm = count_mismatches(window, primer)
            if mm <= max_mismatches:
                sites.append(BindingSite(primer_name, "+", i, i + p_len, mm))

        # "-" site: 3' end of the primer is at the LEFT end of the window
        if window[:three_prime_exact] == primer_rc[:three_prime_exact]:
            mm = count_mismatches(window, primer_rc)
            if mm <= max_mismatches:
                sites.append(BindingSite(primer_name, "-", i, i + p_len, mm))

    return sites