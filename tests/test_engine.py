"""Tests for the AmpliBench engine. Run with:  python -m pytest

Each test uses a template where the right answer is known in advance,
so it checks the science rules, not just that the code runs.
"""

import pytest

from pcr_engine.binding import find_binding_sites
from pcr_engine.conditions import ReactionConditions
from pcr_engine.dimer import find_dimers
from pcr_engine.gel import migration_cm, template_copies
from pcr_engine.program import ThermocyclerProgram
from pcr_engine.recipe import Recipe
from pcr_engine.run import Sample, run_block
from pcr_engine.sequence_utils import clean_sequence, parse_fasta, reverse_complement

FWD = "ATGCGTACGTTAGCATGCAA"
REV = "GGATCCTTAGGCATCCAGTT"
TEMPLATE = "GGGG" + FWD + "T" * 60 + reverse_complement(REV) + "CCCC"   # product = 100 bp
FULL = ["template", "forward primer", "reverse primer", "master mix"]
LOW_TA = ThermocyclerProgram(anneal_temp=40)                           # permissive annealing


# ---------- sequence handling ----------
def test_reverse_complement():
    assert reverse_complement("ATGCGTAAC") == "GTTACGCAT"


def test_clean_sequence_rejects_rna():
    with pytest.raises(ValueError):
        clean_sequence("AUGC")


def test_parse_fasta_multiple_records():
    records = parse_fasta(">a desc\nACGT\nAC\n>b\nGGG\n")
    assert records == {"a": "ACGTAC", "b": "GGG"}


# ---------- primer binding ----------
def test_perfect_binding_positions():
    fwd_sites = find_binding_sites(TEMPLATE, FWD, "F")
    rev_sites = find_binding_sites(TEMPLATE, REV, "R")
    assert [(s.strand, s.start, s.mismatches) for s in fwd_sites] == [("+", 4, 0)]
    assert rev_sites[0].strand == "-" and rev_sites[0].end == len(TEMPLATE) - 4


def test_three_prime_mismatch_blocks_binding():
    pos = 4 + len(FWD) - 1                           # last base of the forward site
    bad = TEMPLATE[:pos] + ("A" if TEMPLATE[pos] != "A" else "C") + TEMPLATE[pos + 1:]
    assert find_binding_sites(bad, FWD, "F") == []


def test_five_prime_mismatch_tolerated():
    bad = TEMPLATE[:4] + ("A" if TEMPLATE[4] != "A" else "C") + TEMPLATE[5:]
    sites = find_binding_sites(bad, FWD, "F", max_mismatches=1)
    assert sites and sites[0].mismatches == 1


# ---------- reaction ----------
def _run(added, program=LOW_TA, template=TEMPLATE):
    sample = Sample("T1", template, FWD, REV, added)
    return run_block([sample], program, ReactionConditions())["T1"]


def test_complete_reaction_gives_expected_size():
    result = _run(FULL)
    assert result.amplified
    assert [p.size for p in result.products] == [100]


def test_product_ends_carry_primer_sequences():
    product = _run(FULL).products[0]
    assert product.sequence.startswith(FWD)
    assert product.sequence.endswith(reverse_complement(REV))


def test_missing_master_mix_means_no_reaction():
    assert not _run(FULL[:3]).amplified


def test_no_template_control_has_no_product():
    assert _run(FULL[1:]).products == []


def test_annealing_above_tm_blocks_amplification():
    assert not _run(FULL, ThermocyclerProgram(anneal_temp=74)).amplified


def test_extension_time_limits_product_size():
    long_template = "GGGG" + FWD + "T" * 1500 + reverse_complement(REV) + "CCCC"
    short_ext = ThermocyclerProgram(anneal_temp=40, extension_sec=30)   # ~500 bp limit
    assert not _run(FULL, short_ext, long_template).amplified


# ---------- primer dimer ----------
def test_complementary_three_prime_ends_form_dimer():
    a = "ACGTACGTACGTAAGGCC"
    b = "TTTTTTTTTTTTTGGCCTT"[:-2]          # ends in GGCC, pairs with a's 3' GGCC
    dimers = find_dimers({"A": a, "B": b}, ReactionConditions())
    assert dimers and dimers[0].overlap >= 4


# ---------- recipe and gel ----------
def test_master_mix_recipe_concentrations():
    c = Recipe(mode="2X master mix", reaction_volume_ul=10, primer_pmol=5).conditions()
    assert c.primer_uM == 0.5 and c.buffer_x == 1.0 and c.mg_mM == 1.5


def test_recipe_overfilled_volume_is_rejected():
    assert Recipe(reaction_volume_ul=10, mm_ul=9, template_ul=2).check()


def test_smaller_fragments_run_further():
    assert migration_cm(100, 1.5, 40) > migration_cm(1000, 1.5, 40)


def test_longer_run_moves_bands_further():
    assert migration_cm(500, 1.5, 40) > migration_cm(500, 1.5, 20)


def test_template_copy_number():
    # 50 ng of a 1500 bp plasmid is about 3 x 10^10 copies
    assert 2.9e10 < template_copies(50, 1500) < 3.2e10


# ---------- gel realism ----------
def test_longer_products_glow_brighter():
    from pcr_engine.gel import size_factor
    assert size_factor(1000) > size_factor(200)


def test_overamplification_smear_only_with_excess_template():
    from pcr_engine.gel import overamplification
    plasmid = template_copies(50, 3000)          # ~1.5e10 copies: plateau early
    genomic = template_copies(50, 3.1e9)          # ~15,000 copies of a single-copy gene
    assert overamplification(35, plasmid) > 0
    assert overamplification(30, genomic) == 0


def test_primer_cloud_present_in_every_lane_with_primers():
    from pcr_engine.run import results_to_lanes
    conditions = ReactionConditions()
    results = run_block([Sample("T1", TEMPLATE, FWD, REV, FULL)], LOW_TA, conditions)
    kinds = [kind for _, _, kind in results_to_lanes(results, LOW_TA, conditions)["T1"]]
    assert "primer" in kinds and "product" in kinds


def test_gradient_holder_takes_one_tube():
    from pcr_engine.service import run_experiment
    full = ["template", "forward primer", "reverse primer", "master mix"]
    data = {"templates": {"t": TEMPLATE}, "pairs": {"P": {"forward": FWD, "reverse": REV}},
            "program": {"anneal_temp": 40}, "gradient": [40.0] * 11,
            "tubes": [{"name": "a", "template": "t", "pair": "P", "well": 3, "added": full},
                      {"name": "b", "template": "t", "pair": "P", "well": 3, "added": full}]}
    result = run_experiment(data)
    assert not result["ok"] and "holder 3" in result["errors"][0]


def test_block_holds_eleven_tubes():
    from pcr_engine.run import MAX_TUBES
    samples = [Sample(f"T{i}", TEMPLATE, FWD, REV, FULL) for i in range(MAX_TUBES + 1)]
    with pytest.raises(ValueError):
        run_block(samples, LOW_TA, ReactionConditions())