"""AmpliBench: virtual PCR bench and agarose gel simulator."""

import io

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from pcr_engine.amplify import primer_tm
from pcr_engine.demo import make_demo
from pcr_engine.dimer import find_dimers
from pcr_engine.gel import (AGAROSE_OPTIONS, LADDERS, MAX_WELLS, draw_gel,
                            template_copies, well_capacity_ul)
from pcr_engine.program import GRADIENT_WELLS, ThermocyclerProgram, gradient_temps
from pcr_engine.recipe import MODES, Recipe
from pcr_engine.run import MAX_TUBES, Sample, results_to_lanes, run_block
from pcr_engine.sequence_utils import clean_sequence, parse_fasta

st.set_page_config(page_title="AmpliBench", layout="wide")

TUBE_COLUMNS = ["Tube", "Template", "Primer pair", "Gradient well", "Template added",
                "Forward added", "Reverse added", "Reaction mix added"]
EMPTY_TUBES = pd.DataFrame({
    "Tube": ["T1"], "Template": [None], "Primer pair": [None], "Gradient well": [1],
    "Template added": [True], "Forward added": [True],
    "Reverse added": [True], "Reaction mix added": [True],
})
DEMO_TUBES = pd.DataFrame({
    "Tube": ["S1", "S2", "S3", "NTC-D", "NoMix"],
    "Template": ["Demo_synthetic_1500bp"] * 5,
    "Primer pair": ["A", "B", "D", "D", "A"],
    "Gradient well": [1, 2, 3, 4, 5],
    "Template added": [True, True, True, False, True],
    "Forward added": [True] * 5,
    "Reverse added": [True] * 5,
    "Reaction mix added": [True, True, True, True, False],
})
EMPTY_PAIRS = pd.DataFrame({"Pair": [""], "Forward": [""], "Reverse": [""]})

st.session_state.setdefault("templates", {})
st.session_state.setdefault("pairs_df", EMPTY_PAIRS)
st.session_state.setdefault("tubes_df", EMPTY_TUBES)


def is_checked(value) -> bool:
    """Checkbox cells can be True, False or empty (NaN) for new rows."""
    return bool(value) if pd.notna(value) else False


def step_row(label, temp_default, time_default, key, temp_range=(30.0, 100.0)):
    """One thermocycler step laid out as: step name | temperature | time."""
    c1, c2, c3 = st.columns([1.3, 1, 1])
    c1.markdown(f"**{label}**")
    temp = c2.number_input("C", temp_range[0], temp_range[1], temp_default, step=0.5,
                           key=f"{key}_t", label_visibility="collapsed")
    secs = c3.number_input("s", 0, 1200, time_default, key=f"{key}_s",
                           label_visibility="collapsed")
    return temp, secs


# =================== LEFT PANEL ===================
with st.sidebar:
    # ---------- reaction setup ----------
    st.header("Reaction setup (per tube)")
    setup_mode = st.radio("Setup", MODES)
    v = float(st.selectbox("Reaction volume (uL)", [10, 20, 25, 50], index=0))

    c1, c2 = st.columns(2)
    primer_pmol = c1.number_input("Each primer (pmol)", 0.5, 100.0, round(0.5 * v, 1),
                                  step=0.5, key=f"pp{v}")
    primer_stock = c2.number_input("Primer stock (pmol/uL)", 1.0, 100.0, 10.0, step=1.0)
    c1, c2 = st.columns(2)
    template_ul = c1.number_input("Template (uL)", 0.0, v, round(0.1 * v, 1), step=0.5,
                                  key=f"tu{v}")
    template_conc = c2.number_input("Template (ng/uL)", 0.0, 5000.0, 50.0, step=10.0)

    recipe = Recipe(mode=setup_mode, reaction_volume_ul=v, primer_pmol=primer_pmol,
                    primer_stock_pmol_per_ul=primer_stock, template_ul=template_ul,
                    template_ng_per_ul=template_conc)
    if setup_mode == "2X master mix":
        recipe.mm_ul = st.number_input("2X master mix (uL)", 0.0, v, v / 2, step=0.5,
                                       key=f"mm{v}")
    else:
        c1, c2 = st.columns(2)
        recipe.buffer_ul = c1.number_input("10X buffer (uL)", 0.0, v, v / 10, step=0.1,
                                           key=f"bu{v}")
        recipe.mg_ul = c2.number_input("25 mM MgCl2 (uL)", 0.0, v, round(1.5 * v / 25, 2),
                                       step=0.1, key=f"mg{v}")
        c1, c2 = st.columns(2)
        recipe.dntp_ul = c1.number_input("10 mM dNTPs (uL)", 0.0, v, round(0.02 * v, 2),
                                         step=0.05, key=f"dn{v}")
        recipe.taq_ul = c2.number_input("Taq 5 U/uL (uL)", 0.0, v, round(0.005 * v, 3),
                                        step=0.01, format="%.3f", key=f"tq{v}")

    conditions = recipe.conditions()
    st.dataframe(pd.DataFrame(recipe.components(), columns=["Component", "Stock", "uL"]),
                 hide_index=True, use_container_width=True)
    st.caption(f"In the tube: {conditions.primer_pmol:g} pmol of each primer "
               f"({conditions.primer_uM} uM), {recipe.template_ng:g} ng template, "
               f"Mg2+ {round(conditions.mg_mM, 2)} mM, dNTPs "
               f"{round(conditions.dntp_mM_each, 2)} mM each.")
    setup_ok = True
    for problem in recipe.check() + conditions.check():
        st.error(problem)
        setup_ok = False
    for note in recipe.notes() + conditions.notes():
        st.warning(note)

    # ---------- thermocycler ----------
    st.header("Thermocycler")
    pcr_mode = st.radio("PCR type", ["Conventional", "Gradient"], horizontal=True)
    st.caption("Step | Temperature (C) | Time (s)")
    init_t, init_s = step_row("Initial denaturation", 95.0, 300, "init", (80.0, 100.0))
    den_t, den_s = step_row("Denaturation", 95.0, 30, "den", (70.0, 100.0))
    if pcr_mode == "Conventional":
        ann_t, ann_s = step_row("Annealing", 58.0, 30, "ann", (30.0, 75.0))
        well_temps = None
    else:
        c1, c2, c3 = st.columns([1.3, 1, 1])
        c1.markdown("**Annealing**")
        c2.markdown("gradient")
        ann_s = c3.number_input("s", 0, 1200, 30, key="ann_s", label_visibility="collapsed")
        c1, c2 = st.columns(2)
        low = c1.number_input("Gradient low (C)", 30.0, 75.0, 50.0, step=0.5)
        high = c2.number_input("Gradient high (C)", 30.0, 75.0, 62.0, step=0.5)
        start = pd.DataFrame({"Well": range(1, GRADIENT_WELLS + 1),
                              "Ta (C)": gradient_temps(low, high)})
        grad = st.data_editor(start, key=f"grad_{low}_{high}", hide_index=True,
                              disabled=["Well"], use_container_width=True)
        well_temps = [float(t) for t in grad["Ta (C)"]]
        ann_t = well_temps[0]
    ext_t, ext_s = step_row("Extension", 72.0, 60, "ext", (40.0, 80.0))
    cycles = st.number_input("Cycles (denaturation, annealing, extension)", 1, 45, 30)
    _, fin_s = step_row("Final extension", 72.0, 300, "fin", (40.0, 80.0))
    c1, c2, c3 = st.columns([1.3, 1, 1])
    c1.markdown("**Hold**")
    c2.markdown("4")
    c3.markdown("until removed")

    program = ThermocyclerProgram(init_t, init_s, den_t, den_s, ann_t, ann_s,
                                  ext_t, ext_s, int(cycles), fin_s)
    st.caption(f"Run time: about {program.total_minutes()} min (excluding ramp time)")

# =================== MAIN PAGE ===================
st.title("AmpliBench")
st.caption("Virtual PCR bench and agarose gel simulator. Outputs are rule-based "
           "simulations of primer binding, not guarantees of wet-lab results.")

# ---------- 1. templates ----------
st.subheader("1. Templates")
left, right = st.columns([3, 1])
with left:
    uploaded = st.file_uploader("Upload FASTA", type=["fa", "fasta", "fna", "txt"])
    pasted = st.text_area("Or paste FASTA", height=110)
with right:
    st.write("")
    if st.button("Load demo data (synthetic)"):
        templates, pairs = make_demo()
        st.session_state.templates = templates
        st.session_state.pairs_df = pd.DataFrame(pairs)
        st.session_state.tubes_df = DEMO_TUBES
        for key in ["pairs_editor", "tubes_editor", "results"]:
            st.session_state.pop(key, None)
        st.rerun()

fasta_text = uploaded.getvalue().decode("utf-8", errors="replace") if uploaded else pasted
if fasta_text.strip():
    try:
        st.session_state.templates = parse_fasta(fasta_text)
    except ValueError as err:
        st.error(f"FASTA problem: {err}")

templates = st.session_state.templates
if templates:
    st.dataframe(pd.DataFrame([
        {"Template": name, "Length (bp)": len(seq),
         "GC (%)": round(100 * (seq.count("G") + seq.count("C")) / len(seq), 1),
         f"Copies in {recipe.template_ng:g} ng": f"{template_copies(recipe.template_ng, len(seq)):.2e}"}
        for name, seq in templates.items()
    ]), hide_index=True, use_container_width=True)
else:
    st.info("Add a template FASTA or load the demo data.")

# ---------- 2. primer pairs ----------
st.subheader("2. Primer pairs")
pairs_table = st.data_editor(st.session_state.pairs_df, num_rows="dynamic",
                             key="pairs_editor", use_container_width=True)

pairs, suggested_ta, tm_rows = {}, {}, []
for _, row in pairs_table.iterrows():
    cells = [row["Pair"], row["Forward"], row["Reverse"]]
    if not all(isinstance(c, str) and c.strip() for c in cells):
        continue
    try:
        fwd, rev = clean_sequence(row["Forward"]), clean_sequence(row["Reverse"])
    except ValueError as err:
        st.error(f"Pair {row['Pair']}: {err}")
        continue
    name = row["Pair"].strip()
    pairs[name] = {"forward": fwd, "reverse": rev}
    tm_f, tm_r = primer_tm(fwd, conditions), primer_tm(rev, conditions)
    suggested_ta[name] = round(min(tm_f, tm_r) - 3, 1)
    dimers = find_dimers({"F": fwd, "R": rev}, conditions)
    tm_rows.append({
        "Pair": name, "Fwd length": len(fwd), "Fwd Tm (C)": tm_f,
        "Rev length": len(rev), "Rev Tm (C)": tm_r,
        "Suggested Ta (C)": suggested_ta[name],
        "3' dimer risk": f"{dimers[0].overlap} bp overlap" if dimers else "none found",
    })
if tm_rows:
    st.dataframe(pd.DataFrame(tm_rows), hide_index=True, use_container_width=True)
    st.caption("Tm uses the reaction setup on the left (nearest-neighbour method).")

# ---------- 3. tubes ----------
st.subheader(f"3. Load the tubes (up to {MAX_TUBES})")
tubes_table = st.data_editor(
    st.session_state.tubes_df, num_rows="dynamic", key="tubes_editor",
    use_container_width=True, column_order=TUBE_COLUMNS,
    column_config={
        "Template": st.column_config.SelectboxColumn(options=list(templates)),
        "Primer pair": st.column_config.SelectboxColumn(options=list(pairs)),
        "Gradient well": st.column_config.NumberColumn(
            min_value=1, max_value=GRADIENT_WELLS, step=1, help="Used only in gradient PCR"),
        "Reaction mix added": st.column_config.CheckboxColumn(
            help="Master mix, or buffer + Mg2+ + dNTPs + polymerase in manual setup"),
    },
)


def tube_ta(row) -> float | None:
    """Annealing temperature this tube actually sees."""
    if well_temps is None:
        return ann_t
    well = row["Gradient well"]
    if pd.isna(well) or not 1 <= int(well) <= GRADIENT_WELLS:
        return None
    return well_temps[int(well) - 1]


for pair_name in sorted({p for p in tubes_table["Primer pair"] if p in suggested_ta}):
    best = suggested_ta[pair_name]
    temps = {tube_ta(r) for _, r in tubes_table.iterrows() if r["Primer pair"] == pair_name}
    temps.discard(None)
    off = [t for t in temps if abs(t - best) > 1.0]
    if off:
        hint = ("use gradient PCR or run this pair separately" if well_temps is None
                else "place it in a well closer to that temperature")
        st.warning(f"Pair {pair_name} works best near {best} C but sees "
                   f"{', '.join(map(str, sorted(off)))} C. Within +/-1 C is ideal; {hint}.")

n_tubes = sum(1 for t in tubes_table["Tube"] if isinstance(t, str) and t.strip())
if n_tubes:
    with st.expander(f"Mix to prepare for {n_tubes} tubes"):
        st.dataframe(pd.DataFrame(recipe.mix_table(n_tubes)), hide_index=True,
                     use_container_width=True)
        st.caption("Add the mix to each tube, then the template (water for the NTC).")

if st.button("Run PCR", type="primary", disabled=not setup_ok):
    samples, problems, tube_temps = [], [], {}
    for _, row in tubes_table.iterrows():
        tube = row["Tube"].strip() if isinstance(row["Tube"], str) else ""
        if not tube:
            continue
        added = [label for label, col in [
            ("template", "Template added"), ("forward primer", "Forward added"),
            ("reverse primer", "Reverse added"), ("master mix", "Reaction mix added"),
        ] if is_checked(row[col])]
        template_seq = templates.get(row["Template"], "")
        if "template" in added and not template_seq:
            problems.append(f"{tube}: choose a template.")
            continue
        pair_name = row["Primer pair"]
        if pair_name not in pairs:
            problems.append(f"{tube}: choose a primer pair.")
            continue
        ta = tube_ta(row)
        if ta is None:
            problems.append(f"{tube}: give a gradient well between 1 and {GRADIENT_WELLS}.")
            continue
        tube_temps[tube] = ta
        copies = template_copies(recipe.template_ng, len(template_seq)) if template_seq else 0
        samples.append(Sample(tube, template_seq, pairs[pair_name]["forward"],
                              pairs[pair_name]["reverse"], added, pair_name,
                              anneal_temp=ta, start_copies=copies))
    for p in problems:
        st.warning(p)
    if samples:
        try:
            st.session_state.results = run_block(samples, program, conditions)
            st.session_state.run_program = program
            st.session_state.run_conditions = conditions
            st.session_state.tube_temps = tube_temps
            st.session_state.run_volume = v
        except ValueError as err:
            st.error(str(err))

# ---------- 4. results ----------
if "results" in st.session_state:
    results = st.session_state.results
    st.subheader("4. PCR results")
    st.dataframe(pd.DataFrame([
        {"Tube": tube, "Ta (C)": st.session_state.tube_temps.get(tube),
         "Amplified": "Yes" if r.amplified else "No",
         "Products (bp)": ", ".join(f"{p.size}" + ("" if p.specific else " (off-target)")
                                    for p in r.products) or "-",
         "Primer dimer": ", ".join(sorted({f"~{d.size} bp" for d in r.dimers})) or "-",
         "Notes": " ".join(r.messages)}
        for tube, r in results.items()
    ]), hide_index=True, use_container_width=True)

    # ---------- 5. gel ----------
    st.subheader("5. Agarose gel")
    lanes = results_to_lanes(results, st.session_state.run_program,
                             st.session_state.run_conditions)

    g1, g2, g3 = st.columns(3)
    agarose = g1.selectbox("Agarose (%)", AGAROSE_OPTIONS, index=AGAROSE_OPTIONS.index(1.5))
    run_choice = g2.radio("Run time", ["Short (20 min)", "Long (40 min)", "Custom"], index=1)
    run_minutes = (g2.slider("Minutes", 5, 90, 30) if run_choice == "Custom"
                   else 20 if run_choice.startswith("Short") else 40)
    ladder_side = g3.radio("Ladder position", ["First well", "Both ends"])
    ladder_position = "both" if ladder_side == "Both ends" else "first"
    needed = len(lanes) + (2 if ladder_position == "both" else 1)

    g1, g2, g3 = st.columns(3)
    wells = g1.number_input("Wells on the gel", needed, MAX_WELLS, needed,
                            help=f"{len(lanes)} samples + ladder = {needed} lanes")
    load_ul = g2.number_input("Load per well (uL of PCR product)", 1.0,
                              float(st.session_state.run_volume),
                              min(5.0, float(st.session_state.run_volume)), step=1.0,
                              help=f"Wells hold about {well_capacity_ul(wells):g} uL each")
    ladder = g3.selectbox("Ladder", list(LADDERS) + ["Custom ladder"])
    custom_ladder = None
    if ladder == "Custom ladder":
        sizes_text = g3.text_input("Band sizes (bp), comma separated", "100, 200, 500, 1000")
        try:
            custom_ladder = [int(s) for s in sizes_text.split(",") if s.strip()]
        except ValueError:
            st.error("Band sizes must be whole numbers separated by commas.")
    else:
        g3.caption(f"Source: {LADDERS[ladder]['source']}. Check your lot's datasheet.")

    g1, g2, g3 = st.columns(3)
    show_sizes = g1.checkbox("Label product sizes", value=True)
    lane_numbers = g2.checkbox("Number lanes (key below gel)", value=True)
    gel_title = g3.text_input("Gel title (optional)")

    try:
        fig, gel_notes = draw_gel(lanes, agarose=agarose, run_minutes=run_minutes,
                                  ladder=ladder, custom_ladder=custom_ladder,
                                  wells=int(wells), ladder_position=ladder_position,
                                  load_ul=load_ul, show_sizes=show_sizes,
                                  lane_numbers=lane_numbers, title=gel_title or None)
    except ValueError as err:
        st.error(str(err))
        st.stop()

    st.pyplot(fig)
    for note in gel_notes:
        st.info(note)

    png, svg = io.BytesIO(), io.BytesIO()
    fig.savefig(png, format="png", dpi=300, bbox_inches="tight")
    fig.savefig(svg, format="svg", bbox_inches="tight")
    plt.close(fig)

    fasta_out = "".join(f">{tube}_{p.size}bp{'' if p.specific else '_offtarget'}\n{p.sequence}\n"
                        for tube, r in results.items() for p in r.products)
    d1, d2, d3 = st.columns(3)
    d1.download_button("Download gel (PNG)", png.getvalue(), "amplibench_gel.png", "image/png")
    d2.download_button("Download gel (SVG)", svg.getvalue(), "amplibench_gel.svg", "image/svg+xml")
    if fasta_out:
        d3.download_button("Download products (FASTA)", fasta_out,
                           "amplibench_products.fasta", "text/plain")