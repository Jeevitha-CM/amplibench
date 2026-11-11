import os
import time
import random

from pcr_engine.sequence_utils import reverse_complement
from pcr_engine.amplify import primer_tm
from pcr_engine.program import ThermocyclerProgram
from pcr_engine.run import Sample, run_block, results_to_lanes
from pcr_engine.gel import draw_gel

# Synthetic 1500 bp test template (fixed seed = same sequence every run)
random.seed(1)
template = "".join(random.choice("ACGT") for _ in range(1500))

# Pair A -> 500 bp product, Pair B -> 1000 bp product
pair_a = {"FwdA": template[100:120], "RevA": reverse_complement(template[580:600])}
pair_b = {"FwdB": template[200:220], "RevB": reverse_complement(template[1180:1200])}

everything = ["template", "forward primer", "reverse primer", "master mix"]

# Annealing temperature = lowest primer Tm minus 3 C (a common lab rule)
all_primers = {**pair_a, **pair_b}
for name, seq in all_primers.items():
    print(name, seq, primer_tm(seq))
anneal = round(min(primer_tm(s) for s in all_primers.values()) - 3, 1)
print("Annealing at", anneal, "C")

samples = [
    Sample("S1", template, pair_a, everything),
    Sample("S2", template, pair_b, everything),
    Sample("Mix", template, all_primers, everything),
    Sample("NTC", "", pair_a, ["forward primer", "reverse primer", "master mix"]),
    Sample("NoMM", template, pair_a, ["template", "forward primer", "reverse primer"]),
]

program = ThermocyclerProgram(anneal_temp=anneal, cycles=30)
results = run_block(samples, program)

for name, r in results.items():
    sizes = [p.size for p in r.products]
    print(name, "->", sizes, r.messages)

lanes = results_to_lanes(results, program)

os.makedirs("outputs", exist_ok=True)
stamp = time.strftime("%H%M%S")

fig = draw_gel(lanes, agarose=1.5, ladder="100 bp ladder", title="Block run 1")
fig.savefig(f"outputs/gel_100bp_{stamp}.png", dpi=150)

fig = draw_gel(lanes, agarose=0.7, ladder="100 bp ladder", show_sizes=False)
fig.savefig(f"outputs/gel_bad_choice_{stamp}.png", dpi=150)

print("Saved gels in the outputs folder with stamp", stamp)