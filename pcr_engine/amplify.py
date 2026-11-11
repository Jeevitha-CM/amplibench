"""Turn binding sites into PCR products and check the reaction setup."""

from dataclasses import dataclass, field

from Bio.SeqUtils import MeltingTemp as mt

from .binding import BindingSite, find_binding_sites
from .conditions import ReactionConditions
from .dimer import Dimer, find_dimers
from .program import ThermocyclerProgram
from .sequence_utils import reverse_complement


@dataclass
class Product:
    left: BindingSite      # "+" site, where extension goes to the right
    right: BindingSite     # "-" site, where extension goes to the left
    size: int
    sequence: str

    @property
    def mismatches(self) -> int:
        return self.left.mismatches + self.right.mismatches

    @property
    def specific(self) -> bool:
        return self.mismatches == 0


@dataclass
class PCRResult:
    amplified: bool
    template_present: bool = False
    start_copies: float = 0.0          # template copies in the tube (set by run_block)
    primer_lengths: list[int] = field(default_factory=list)   # primers present in the tube
    products: list[Product] = field(default_factory=list)
    dimers: list[Dimer] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)


def primer_tm(primer: str, conditions: ReactionConditions | None = None) -> float:
    """Nearest-neighbour Tm using the actual reaction conditions."""
    c = conditions or ReactionConditions()
    return round(mt.Tm_NN(primer, dnac1=c.primer_nM, dnac2=0, Na=50,
                          Mg=c.mg_mM, dNTPs=4 * c.dntp_mM_each), 1)


def allowed_mismatches(tm: float, anneal_temp: float) -> int:
    """Lower annealing temperature (relative to Tm) tolerates more mismatches.

    Within 2 C of Tm: 0, 2-6 C below: 1, 6-10 C: 2, 10-14 C: 3, more: 4.
    Illustrative rule, not a thermodynamic calculation.
    """
    gap = tm - anneal_temp
    return min(4, int((gap + 2) // 4))


def run_pcr(template: str, primers: dict[str, str], added: list[str],
            program: ThermocyclerProgram,
            conditions: ReactionConditions | None = None) -> PCRResult:
    """Simulate one PCR tube. `primers` holds only the primers actually added."""
    c = conditions or ReactionConditions()
    result = PCRResult(amplified=False, template_present=bool(template),
                       primer_lengths=[len(s) for s in primers.values()])

    # 1. Master mix carries buffer, Mg2+, dNTPs and polymerase
    if "master mix" not in added:
        result.messages.append("Reaction mix missing or incomplete (buffer, Mg2+, dNTPs, polymerase): no reaction.")
        return result

    # 2. Machine program and reaction mix
    problems = program.check() + c.check()
    if problems:
        result.messages.extend(problems)
        return result

    if not primers:
        result.messages.append("No primers added: nothing can be amplified.")
        return result

    # 3. Primer dimers can form with or without template
    result.dimers = find_dimers(primers, c)

    if not template:
        result.messages.append("No template: only primer dimers can form."
                               if result.dimers else "No template: no amplification.")
        return result

    # 4. Which primers can stay bound, and how many mismatches each tolerates
    usable = {}
    for name, seq in primers.items():
        tm = primer_tm(seq, c)
        if program.anneal_temp > tm:
            result.messages.append(f"{name} Tm is {tm} C, below the annealing temperature "
                                   f"({program.anneal_temp} C). It cannot stay bound.")
            continue
        usable[name] = allowed_mismatches(tm, program.anneal_temp)

    # 5. Binding sites
    sites = []
    for name, max_mm in usable.items():
        sites += find_binding_sites(template, primers[name], name, max_mismatches=max_mm)
    plus_sites = [s for s in sites if s.strand == "+"]
    minus_sites = [s for s in sites if s.strand == "-"]

    # 6. Pair sites that face each other and fit the extension time
    max_size = int(program.extension_sec / 60 * 1000)
    too_long = 0
    for left in plus_sites:
        for right in minus_sites:
            if right.start < left.end:
                continue
            size = right.end - left.start
            if size > max_size:
                too_long += 1
                continue
            product_seq = (primers[left.primer_name]
                           + template[left.end:right.start]
                           + reverse_complement(primers[right.primer_name]))
            result.products.append(Product(left, right, size, product_seq))

    if too_long:
        result.messages.append(f"{too_long} possible product(s) longer than {max_size} bp "
                               f"were skipped (extension time too short).")

    result.products.sort(key=lambda p: p.size)
    result.amplified = len(result.products) > 0
    off_target = sum(1 for p in result.products if not p.specific)
    result.messages.append(f"{len(result.products)} product(s) formed"
                           + (f", {off_target} from mismatched (off-target) binding." if off_target else "."))
    if result.dimers:
        result.messages.append(f"Primer dimer possible ({result.dimers[0].overlap} bp 3' overlap).")
    return result