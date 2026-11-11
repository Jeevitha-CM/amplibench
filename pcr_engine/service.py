"""One function that runs a whole experiment from plain data (dicts and lists).

The 3D bench sends JSON to the API, the API passes it here, and the answer goes
back as JSON. Keeping this separate from the web framework means it can be tested
without a server.
"""

import base64
import io

import matplotlib
matplotlib.use("Agg")              # draw images without a screen (server use)
import matplotlib.pyplot as plt

from .amplify import primer_tm
from .demo import make_demo
from .gel import (GEL_LENGTH_CM, LADDERS, STANDARD_LOAD_UL, _ladder_bands, draw_gel,
                  half_speed_size, migration_cm, template_copies)
from .program import GRADIENT_WELLS, ThermocyclerProgram
from .recipe import Recipe
from .run import Sample, results_to_lanes, run_block
from .sequence_utils import clean_sequence

# Reagents that make up the "reaction mix" in each setup mode
MIX_COMPONENTS = {
    "2X master mix": ["master mix"],
    "Manual components": ["buffer", "mgcl2", "dntps", "taq"],
}


def demo_payload() -> dict:
    """Demo templates and primer pairs, plus the choices the bench needs."""
    templates, pairs = make_demo()
    return {
        "templates": templates,
        "pairs": {p: {"forward": f, "reverse": r}
                  for p, f, r in zip(pairs["Pair"], pairs["Forward"], pairs["Reverse"])},
        "ladders": list(LADDERS),
        "gradient_wells": GRADIENT_WELLS,
    }


def primer_info(pairs: dict, recipe: Recipe) -> dict:
    """Tm and suggested annealing temperature of each pair under this recipe."""
    conditions = recipe.conditions()
    info = {}
    for name, p in pairs.items():
        tm_f = primer_tm(clean_sequence(p["forward"]), conditions)
        tm_r = primer_tm(clean_sequence(p["reverse"]), conditions)
        info[name] = {"tm_forward": tm_f, "tm_reverse": tm_r,
                      "suggested_ta": round(min(tm_f, tm_r) - 3, 1)}
    return info


def gel_layout(lanes: dict, agarose: float, run_minutes: float, ladder: str,
               ladder_position: str, load_ul: float) -> dict:
    """Where every band ends up (cm from the well), for drawing the 3D gel.

    Uses the same migration model as draw_gel, so the 3D gel and the image agree.
    """
    sizes, _, brightness = _ladder_bands(ladder, None)
    wells = [("ladder", "M")] + [("sample", name) for name in lanes]
    if ladder_position == "both":
        wells.append(("ladder", "M"))
    factor = (load_ul / STANDARD_LOAD_UL) ** 0.5

    out = []
    for kind, label in wells:
        if kind == "ladder":
            bands = [{"bp": bp, "alpha": round(brightness[bp], 2), "kind": "ladder"} for bp in sizes]
        else:
            bands = [{"bp": b[0], "alpha": round(min(1.0, b[1] * factor), 2), "kind": b[2]}
                     for b in lanes[label] if b[1] > 0]
        for band in bands:
            band["cm"] = round(migration_cm(band["bp"], agarose, run_minutes), 3)
        out.append({"kind": kind, "label": label, "bands": bands})

    dye_cm = migration_cm(0.7 * half_speed_size(agarose), agarose, run_minutes)
    return {"wells": out, "gel_length_cm": GEL_LENGTH_CM, "dye_cm": round(dye_cm, 3)}


def run_experiment(data: dict) -> dict:
    """Run tubes through PCR and a gel. `data` is the JSON sent by the bench."""
    recipe = Recipe(**data.get("recipe", {}))
    problems = recipe.check() + recipe.conditions().check()
    if problems:
        return {"ok": False, "errors": problems}

    program = ThermocyclerProgram(**data.get("program", {}))
    gradient = data.get("gradient")            # list of well temperatures, or None
    templates = {k: clean_sequence(v) for k, v in data["templates"].items()}
    pairs = data["pairs"]
    mix_parts = MIX_COMPONENTS[recipe.mode]

    samples, errors, tube_ta = [], [], {}
    used_wells = set()                         # one tube per gradient holder
    for i, tube in enumerate(data["tubes"]):
        name = tube["name"]
        # The forward and reverse primer can come from different pairs if that is what was pipetted
        fwd_name = tube.get("forward_pair", tube.get("pair"))
        rev_name = tube.get("reverse_pair", tube.get("pair"))
        if fwd_name not in pairs or rev_name not in pairs:
            errors.append(f"{name}: choose a primer pair.")
            continue
        label = fwd_name if fwd_name == rev_name else f"{fwd_name}/{rev_name}"

        # Translate what was pipetted into what the engine understands
        added_reagents = set(tube.get("added", []))
        added = []
        if "template" in added_reagents:
            added.append("template")
        if "forward primer" in added_reagents:
            added.append("forward primer")
        if "reverse primer" in added_reagents:
            added.append("reverse primer")
        if all(part in added_reagents for part in mix_parts):
            added.append("master mix")

        ta = program.anneal_temp
        if gradient:
            well = tube.get("well", i + 1)
            if not 1 <= well <= len(gradient):
                errors.append(f"{name}: gradient well must be 1-{len(gradient)}.")
                continue
            if well in used_wells:
                errors.append(f"{name}: gradient holder {well} already has a tube in it.")
                continue
            used_wells.add(well)
            ta = float(gradient[well - 1])
        tube_ta[name] = ta

        template = templates.get(tube.get("template"), "")
        copies = template_copies(recipe.template_ng, len(template)) if template else 0
        samples.append(Sample(name, template, clean_sequence(pairs[fwd_name]["forward"]),
                              clean_sequence(pairs[rev_name]["reverse"]), added, label,
                              anneal_temp=ta, start_copies=copies))

    if errors or not samples:
        return {"ok": False, "errors": errors or ["No tubes to run."]}

    conditions = recipe.conditions()
    results = run_block(samples, program, conditions)
    lanes = results_to_lanes(results, program, conditions)

    gel = data.get("gel", {})
    settings = {"agarose": gel.get("agarose", 1.5),
                "run_minutes": gel.get("run_minutes", 40),
                "ladder": gel.get("ladder", list(LADDERS)[0]),
                "ladder_position": gel.get("ladder_position", "first"),
                "load_ul": gel.get("load_ul", 5.0)}
    fig, notes = draw_gel(lanes, title=gel.get("title"), **settings)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    return {
        "ok": True,
        "results": [
            {"tube": name, "ta": tube_ta[name], "amplified": r.amplified,
             "products": [{"size": p.size, "specific": p.specific} for p in r.products],
             "dimers": sorted({d.size for d in r.dimers}),
             "messages": r.messages}
            for name, r in results.items()
        ],
        "gel_png": base64.b64encode(buffer.getvalue()).decode("ascii"),
        "gel_layout": gel_layout(lanes, **settings),
        "notes": notes,
        "run_minutes": program.total_minutes(),
    }