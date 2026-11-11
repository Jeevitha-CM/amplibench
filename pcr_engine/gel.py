"""Simulated agarose gel run: wells, run time, migration, ladders and band brightness.

The migration model is illustrative. It reproduces the qualitative behaviour of real
gels (small DNA runs fast, separation grows with run time, higher agarose resolves
smaller fragments, small fragments can run off the gel) but is not calibrated to any
specific gel box, buffer or stain.
"""

import math

import matplotlib.pyplot as plt

# ---------- gel setup (a typical mini gel) ----------
GEL_LENGTH_CM = 7.5            # distance from the wells to the bottom edge
FIELD_V_PER_CM = 5.0           # about 100 V on a mini gel tank (stated, not varied)
MAX_SPEED_CM_PER_MIN = 0.26    # speed of very small DNA at this field (illustrative)
SHAPE = 0.75                   # how steeply speed falls with size (illustrative)

AGAROSE_OPTIONS = [0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
MAX_WELLS = 20


def well_capacity_ul(wells: int) -> float:
    """Approximate volume one well holds on a mini gel: more wells = narrower wells."""
    if wells <= 8:
        return 25
    if wells <= 10:
        return 20
    if wells <= 12:
        return 15
    return 10
STANDARD_LOAD_UL = 5.0   # brightness model assumes 5 uL of product loaded

# ---------- ladders ----------
# Commercial values are taken from manufacturer documentation (see README).
# "mass_ng" = DNA per band from the manufacturer's figure, used for brightness.
LADDERS = {
    "100 bp DNA Ladder (NEB)": {
        "bands": [1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100],
        "reference": [1000, 517, 500],
        "mass_ng": {1517: 45, 1200: 35, 1000: 95, 900: 27, 800: 24, 700: 21, 600: 18,
                    517: 48.5, 500: 48.5, 400: 38, 300: 29, 200: 25, 100: 48},
        "recommended": (1.2, 3.0),
        "source": "New England Biolabs, N3231 product figure (0.5 ug/lane; 500/517 doublet = 97 ng, split equally)",
    },
    "1 kb Plus DNA Ladder (NEB)": {
        "bands": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200,
                  1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        "reference": [3000, 1000, 500],
        "mass_ng": {10000: 40, 8000: 40, 6000: 48, 5000: 40, 4000: 32, 3000: 120,
                    2000: 40, 1500: 57, 1200: 45, 1000: 122, 900: 34, 800: 31,
                    700: 27, 600: 23, 500: 124, 400: 49, 300: 37, 200: 32, 100: 61},
        "recommended": (1.0, 1.4),
        "source": "New England Biolabs, N3200 product sheet (1 ug/lane; sizes rounded)",
    },
    "GeneRuler 50 bp (Thermo Scientific)": {
        "bands": [1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50],
        "reference": [500, 250],
        "mass_ng": None,
        "recommended": None,
        "source": "Thermo Scientific product information, SM0371",
    },
    "GeneRuler 100 bp (Thermo Scientific)": {
        "bands": [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        "reference": [500],
        "mass_ng": None,
        "recommended": (1.0, 2.0),
        "source": "Thermo Scientific product information, SM0241/SM0243",
    },
    "GeneRuler 1 kb (Thermo Scientific)": {
        "bands": [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500,
                  2000, 1500, 1000, 750, 500, 250],
        "reference": [6000, 3000, 1000],
        "mass_ng": None,
        "recommended": (1.0, 1.0),
        "source": "Thermo Scientific product information, SM0311",
    },
}

# ---------- band brightness from PCR cycles (illustrative) ----------
START_COPIES = 1e4
PLATEAU_COPIES = 1e12
VISIBLE_COPIES = 1e10


def template_copies(template_ng: float, template_bp: int) -> float:
    """Copies of template: mass / (length x 650 g/mol per bp) x Avogadro's number."""
    if template_ng <= 0 or template_bp <= 0:
        return 0.0
    return template_ng * 1e-9 / (template_bp * 650) * 6.022e23


def band_intensity(cycles: int, mismatches: int = 0, start_copies: float = START_COPIES) -> float:
    """Return band brightness from 0 (invisible) to 1 (saturated). Illustrative model."""
    efficiency = max(0.3, 0.9 - 0.15 * mismatches)
    growth = (1 + efficiency) ** cycles
    if start_copies <= 0:
        return 0.0
    start = min(start_copies, PLATEAU_COPIES)
    copies = PLATEAU_COPIES / (1 + (PLATEAU_COPIES / start - 1) / growth)

    if copies < VISIBLE_COPIES:
        return 0.0
    fraction = ((math.log10(copies) - math.log10(VISIBLE_COPIES))
                / (math.log10(PLATEAU_COPIES) - math.log10(VISIBLE_COPIES)))
    return round(0.15 + 0.85 * fraction, 2)


def size_factor(bp: int) -> float:
    """DNA stains bind along the length: at equal copies, longer products glow brighter.

    Brightness is scaled to a 500 bp product, with a gentle (square-root) curve because
    staining and camera exposure saturate. Illustrative.
    """
    return min(1.4, max(0.35, (bp / 500) ** 0.5))


def overamplification(cycles: int, start_copies: float) -> float:
    """Smear strength (0 to 0.45) when the reaction sits at its plateau for many cycles.

    Too much template or too many cycles: products re-prime on each other and form a smear.
    """
    if start_copies <= 0:
        return 0.0
    cycles_to_plateau = math.log(PLATEAU_COPIES / min(start_copies, PLATEAU_COPIES / 2)) / math.log(1.9)
    extra = cycles - cycles_to_plateau
    return round(max(0.0, min(1.0, (extra - 8) / 10)) * 0.45, 2)


# ---------- migration ----------
def half_speed_size(agarose: float) -> float:
    """Fragment size (bp) that moves at half the maximum speed. Denser gel = smaller."""
    return 500 * agarose ** -1.75


def migration_cm(bp: float, agarose: float, minutes: float) -> float:
    """Distance travelled from the well: speed x time, speed falls with size."""
    relative_speed = 1 / (1 + (bp / half_speed_size(agarose)) ** SHAPE)
    return MAX_SPEED_CM_PER_MIN * minutes * relative_speed


def _ladder_bands(ladder: str, custom_ladder: list[int] | None):
    """Return (sizes, reference sizes, {size: brightness}) for the chosen ladder."""
    if ladder == "Custom ladder":
        if not custom_ladder or any(bp <= 0 for bp in custom_ladder):
            raise ValueError("Enter at least one positive band size for the custom ladder.")
        sizes = sorted(set(custom_ladder), reverse=True)
        return sizes, [], {bp: 0.7 for bp in sizes}
    if ladder not in LADDERS:
        raise ValueError(f"Unknown ladder '{ladder}'. Choose from: {list(LADDERS)}")

    info = LADDERS[ladder]
    refs = info["reference"]
    if info["mass_ng"]:
        top = max(info["mass_ng"].values())
        brightness = {bp: 0.35 + 0.65 * info["mass_ng"][bp] / top for bp in info["bands"]}
    else:
        brightness = {bp: 1.0 if bp in refs else 0.65 for bp in info["bands"]}
    return info["bands"], refs, brightness


def _draw_band(ax, x, y, height, alpha, fuzzy=False):
    """A band = a bright core plus a soft halo. Primer dimers are drawn as a diffuse smear."""
    if fuzzy:
        layers = [(3.0, 0.10), (2.2, 0.15), (1.5, 0.22), (1.0, 0.30)]
    else:
        layers = [(1.6, 0.20), (1.0, 1.0)]
    for scale, strength in layers:
        h = height * scale
        ax.add_patch(plt.Rectangle((x - 0.3, y - h / 2), 0.6, h, color="white",
                                   alpha=alpha * strength, linewidth=0))


def draw_gel(lanes: dict[str, list[tuple]],
             agarose: float = 1.0,
             run_minutes: float = 40,
             ladder: str = "GeneRuler 1 kb (Thermo Scientific)",
             custom_ladder: list[int] | None = None,
             wells: int | None = None,
             ladder_position: str = "first",
             load_ul: float = STANDARD_LOAD_UL,
             show_sizes: bool = True,
             lane_numbers: bool = True,
             title: str | None = None):
    """Draw a simulated gel. Returns (figure, list of notes about the run).

    lanes = {"Sample": [(size_bp, brightness) or (size_bp, brightness, kind), ...]}
    kind is "product", "off-target" or "dimer".
    """
    notes = []
    sizes, refs, ladder_brightness = _ladder_bands(ladder, custom_ladder)
    names = list(lanes.keys())

    # ---- assign wells ----
    ladder_lanes = 2 if ladder_position == "both" else 1
    needed = len(names) + ladder_lanes
    comb = needed if wells is None else wells           # by default: exactly the lanes needed
    if needed > comb:
        raise ValueError(f"{needed} lanes needed (samples + ladder), but the gel has {comb} wells.")
    if comb > MAX_WELLS:
        raise ValueError(f"At most {MAX_WELLS} wells on one gel.")

    wells = [None] * comb                     # None = empty well
    wells[0] = "LADDER"
    for i, name in enumerate(names):
        wells[i + 1] = name
    if ladder_position == "both":
        wells[len(names) + 1] = "LADDER"

    # ---- recommended agarose for this ladder ----
    if ladder in LADDERS and LADDERS[ladder]["recommended"]:
        low, high = LADDERS[ladder]["recommended"]
        if not low <= agarose <= high:
            notes.append(f"{ladder} is recommended on {low}-{high}% agarose; "
                         f"you are using {agarose}%.")

    if load_ul > well_capacity_ul(comb):
        notes.append(f"{load_ul:g} uL is more than a {comb}-well comb holds "
                     f"(about {well_capacity_ul(comb):g} uL): sample would overflow into "
                     f"neighbouring wells.")
    # more DNA loaded = brighter band, but staining saturates, so the gain is gentle
    load_factor = (load_ul / STANDARD_LOAD_UL) ** 0.5

    band_height = 0.05 + 0.015 * math.sqrt(run_minutes / 20)   # bands broaden with time

    fig, ax = plt.subplots(figsize=(0.75 * comb + 3.2, 7.5))
    fig.patch.set_facecolor("white")
    ax.add_patch(plt.Rectangle((-0.5, -0.15), comb, GEL_LENGTH_CM + 0.15, color="black"))

    for i in range(1, comb):
        ax.plot([i - 0.5, i - 0.5], [0, GEL_LENGTH_CM], color="#444444",
                linestyle="--", linewidth=0.5)

    lost = {}
    visible_product_sizes = set()
    dimer_on_gel = False

    for x, content in enumerate(wells):
        ax.add_patch(plt.Rectangle((x - 0.3, -0.12), 0.6, 0.12, color="#2b2b2b"))

        if content == "LADDER":
            bands = [(bp, ladder_brightness[bp], "ladder") for bp in sizes]
            header = "M"
        elif content is None:
            bands = []
            header = str(x + 1) if lane_numbers else ""
        else:
            bands = [(b[0], min(1.0, b[1] * load_factor), b[2] if len(b) == 3 else "product")
                     for b in lanes[content]]
            header = str(x + 1) if lane_numbers else content

        for bp, alpha, kind in bands:
            if alpha <= 0:
                continue
            y = migration_cm(bp, agarose, run_minutes)
            if kind == "smear":                      # haze from the well down past the product
                bottom = min(GEL_LENGTH_CM, y * 1.15)
                for k in range(24):
                    top = bottom * k / 24
                    ax.add_patch(plt.Rectangle((x - 0.28, top), 0.56, bottom / 24, color="white",
                                               alpha=alpha * 0.35 * (0.3 + 0.7 * k / 24), linewidth=0))
                continue
            if y > GEL_LENGTH_CM:
                if kind not in ("primer", "dimer"):
                    lost.setdefault("Ladder" if kind == "ladder" else content, []).append(bp)
                continue
            fuzzy = kind in ("dimer", "primer")
            height = band_height * (2.4 if kind == "primer" else 2 if kind == "dimer" else 1)
            _draw_band(ax, x, y, height, alpha, fuzzy=fuzzy)
            if fuzzy:
                dimer_on_gel = min(y, dimer_on_gel) if dimer_on_gel else y
            elif kind != "ladder":
                visible_product_sizes.add(bp)

        ax.text(x, -0.35, header, ha="center", va="bottom", fontsize=11, fontweight="bold")

    # ---- ladder sizes on the left (reference bands in bold, crowded labels skipped) ----
    last_y = None
    for bp in sorted(sizes, reverse=True):
        y = migration_cm(bp, agarose, run_minutes)
        if y > GEL_LENGTH_CM or (last_y is not None and y - last_y < 0.2):
            continue
        ax.text(-0.6, y, str(bp), ha="right", va="center", fontsize=8.5,
                fontweight="bold" if bp in refs else "normal")
        last_y = y

    # ---- product sizes on the right ----
    if show_sizes:
        last_y = None
        for bp in sorted(visible_product_sizes, reverse=True):
            y = migration_cm(bp, agarose, run_minutes)
            if last_y is not None and y - last_y < 0.22:
                continue
            ax.text(comb - 0.35, y, f"{bp} bp", ha="left", va="center",
                    fontsize=9.5, fontweight="bold")
            last_y = y
        if dimer_on_gel:
            ax.text(comb - 0.35, dimer_on_gel, "primers /\ndimers",
                    ha="left", va="bottom", fontsize=8, style="italic", color="#555555")

    # ---- lane key under the gel ----
    if lane_numbers:
        key = [f"{x + 1}: {'M' if c == 'LADDER' else c}"
               for x, c in enumerate(wells) if c is not None]
        ax.text((comb - 1) / 2, GEL_LENGTH_CM + 0.35,
                "   ".join(key) + f"\nM = {ladder}", ha="center", va="top", fontsize=8.5)

    for where, bps in lost.items():
        notes.append(f"{where}: {', '.join(map(str, sorted(bps)))} bp ran off the gel.")

    ax.set_title(title or f"{agarose}% agarose, {run_minutes:g} min at "
                          f"{FIELD_V_PER_CM:g} V/cm (simulated)", fontsize=11, pad=12)
    ax.set_xlim(-1.5, comb + 2.6)
    ax.set_ylim(GEL_LENGTH_CM + 1.3, -0.6)
    ax.axis("off")
    return fig, notes