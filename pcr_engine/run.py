"""Run several PCR tubes with one thermocycler program, like a real block."""

from dataclasses import dataclass, field, replace

from .amplify import PCRResult, run_pcr
from .conditions import ReactionConditions
from .dimer import dimer_intensity
from .gel import START_COPIES, band_intensity, overamplification, size_factor
from .program import ThermocyclerProgram

MAX_TUBES = 11   # our thermocycler has 11 tube holders


@dataclass
class Sample:
    name: str
    template: str
    forward: str
    reverse: str
    added: list[str] = field(default_factory=list)
    pair_name: str = "P"
    anneal_temp: float | None = None   # set in gradient mode: this tube's own Ta
    start_copies: float = START_COPIES  # template copies put in the tube


def run_block(samples: list[Sample], program: ThermocyclerProgram,
              conditions: ReactionConditions) -> dict[str, PCRResult]:
    """Run every tube with the same program and reaction mix."""
    if len(samples) > MAX_TUBES:
        raise ValueError(f"The block holds {MAX_TUBES} tubes; you loaded {len(samples)}.")

    names = [s.name for s in samples]
    if len(names) != len(set(names)):
        raise ValueError("Each tube needs a unique name.")

    results = {}
    for s in samples:
        primers = {}
        if "forward primer" in s.added:
            primers[f"{s.pair_name}_F"] = s.forward
        if "reverse primer" in s.added:
            primers[f"{s.pair_name}_R"] = s.reverse
        template = s.template if "template" in s.added else ""
        tube_program = program
        if s.anneal_temp is not None:
            tube_program = replace(program, anneal_temp=s.anneal_temp)
        result = run_pcr(template, primers, s.added, tube_program, conditions)
        result.start_copies = s.start_copies if template else 0
        results[s.name] = result
    return results


def results_to_lanes(results: dict[str, PCRResult], program: ThermocyclerProgram,
                     conditions: ReactionConditions) -> dict[str, list[tuple[int, float, str]]]:
    """Gel lanes: tube name -> [(size_bp, brightness, kind), ...]."""
    lanes = {}
    for name, result in results.items():
        bands = []
        for p in result.products:
            kind = "product" if p.specific else "off-target"
            copies = result.start_copies
            value = band_intensity(program.cycles, p.mismatches, copies) * size_factor(p.size)
            bands.append((p.size, round(min(1.0, value), 2), kind))

        # Smear from overamplification, under the main product
        if result.products:
            smear = overamplification(program.cycles, result.start_copies)
            if smear > 0:
                bands.append((max(p.size for p in result.products), smear, "smear"))

        # Unused primers: a faint cloud at the bottom, weaker when primers were used up
        if result.primer_lengths:
            cloud = min(0.35, 0.18 * conditions.primer_uM / 0.4)
            if result.amplified:
                cloud *= 0.6
            bands.append((max(result.primer_lengths), round(cloud, 2), "primer"))

        dimer_sizes = {}
        for d in result.dimers:
            value = dimer_intensity(d, conditions, result.template_present,
                                    bool(result.products))
            dimer_sizes[d.size] = max(value, dimer_sizes.get(d.size, 0))
        bands += [(size, value, "dimer") for size, value in dimer_sizes.items()]

        lanes[name] = bands
    return lanes