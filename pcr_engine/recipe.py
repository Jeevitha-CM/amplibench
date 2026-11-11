"""Reaction setup: turn what you pipette into final amounts in the tube.

Two ways to set up a PCR:
  * "2X master mix": ready mix (buffer, Mg2+, dNTPs, polymerase) + primers + template + water
  * "Manual components": each component added separately from its stock
"""

from dataclasses import dataclass

from .conditions import ReactionConditions

MODES = ["2X master mix", "Manual components"]

# Typical 1X composition of a 2X Taq master mix. Kits differ; the mix is used as supplied.
MM_MG_MM_AT_1X = 1.5
MM_DNTP_MM_EACH_AT_1X = 0.2
MM_TAQ_UNITS_PER_UL_AT_1X = 0.025


@dataclass
class Recipe:
    mode: str = "2X master mix"
    reaction_volume_ul: float = 10.0

    # both modes
    primer_pmol: float = 5.0                  # pmol of EACH primer in the tube
    primer_stock_pmol_per_ul: float = 10.0    # 10 pmol/uL = 10 uM working stock
    template_ul: float = 1.0
    template_ng_per_ul: float = 50.0

    # 2X master mix mode
    mm_ul: float = 5.0

    # manual mode
    buffer_ul: float = 1.0          # 10X buffer without Mg
    mg_ul: float = 0.6              # 25 mM MgCl2
    dntp_ul: float = 0.2            # 10 mM each dNTP
    taq_ul: float = 0.05            # 5 U/uL Taq

    @property
    def primer_ul(self) -> float:
        """Volume of each primer to add = pmol wanted / stock pmol per uL."""
        return self.primer_pmol / self.primer_stock_pmol_per_ul

    @property
    def template_ng(self) -> float:
        return self.template_ul * self.template_ng_per_ul

    def components(self) -> list[tuple[str, str, float]]:
        """(component, stock, volume per reaction in uL), water last."""
        if self.mode == "2X master mix":
            rows = [("2X master mix", "2X", self.mm_ul)]
        else:
            rows = [("PCR buffer (no Mg)", "10X", self.buffer_ul),
                    ("MgCl2", "25 mM", self.mg_ul),
                    ("dNTP mix", "10 mM each", self.dntp_ul),
                    ("Taq polymerase", "5 U/uL", self.taq_ul)]
        stock = f"{self.primer_stock_pmol_per_ul:g} pmol/uL"
        rows += [("Forward primer", stock, round(self.primer_ul, 2)),
                 ("Reverse primer", stock, round(self.primer_ul, 2)),
                 ("Template DNA", f"{self.template_ng_per_ul:g} ng/uL", self.template_ul),
                 ("Nuclease-free water", "-", round(self.water_ul(), 2))]
        return rows

    def water_ul(self) -> float:
        """Water makes up the rest of the reaction volume."""
        used = 2 * self.primer_ul + self.template_ul
        if self.mode == "2X master mix":
            used += self.mm_ul
        else:
            used += self.buffer_ul + self.mg_ul + self.dntp_ul + self.taq_ul
        return self.reaction_volume_ul - used

    def conditions(self) -> ReactionConditions:
        """Final concentrations = stock x volume added / total volume."""
        v = self.reaction_volume_ul
        if self.mode == "2X master mix":
            strength = 2 * self.mm_ul / v          # 1.0 when the mix is half the volume
            return ReactionConditions(
                v, self.primer_pmol,
                mg_mM=MM_MG_MM_AT_1X * strength,
                dntp_mM_each=MM_DNTP_MM_EACH_AT_1X * strength,
                buffer_x=strength,
                taq_units=MM_TAQ_UNITS_PER_UL_AT_1X * strength * v)
        return ReactionConditions(
            v, self.primer_pmol,
            mg_mM=25 * self.mg_ul / v,
            dntp_mM_each=10 * self.dntp_ul / v,
            buffer_x=10 * self.buffer_ul / v,
            taq_units=5 * self.taq_ul)

    def check(self) -> list[str]:
        """Setup problems that make the recipe impossible."""
        if self.water_ul() < -1e-9:
            used = self.reaction_volume_ul - self.water_ul()
            return [f"Components add up to {round(used, 2)} uL, more than the "
                    f"{self.reaction_volume_ul:g} uL reaction."]
        return []

    def notes(self) -> list[str]:
        notes = []
        if self.mode == "2X master mix" and abs(2 * self.mm_ul - self.reaction_volume_ul) > 0.01:
            notes.append(f"A 2X mix should be half the reaction: {self.reaction_volume_ul / 2:g} uL.")
        volumes = [vol for _, _, vol in self.components()[:-1] if vol > 0]
        if any(vol < 0.5 for vol in volumes):
            notes.append("Volumes below 0.5 uL cannot be pipetted accurately: prepare a mix "
                         "for all tubes (see the mix table) and add it to each tube.")
        return notes

    def mix_table(self, n_tubes: int, extra: float = 0.1) -> list[dict]:
        """Volumes for a shared mix of everything except template, plus 10% extra."""
        factor = n_tubes * (1 + extra)
        return [{"Component": name, "Stock": stock, "Per reaction (uL)": round(vol, 2),
                 f"Mix for {n_tubes} tubes + 10% (uL)": round(vol * factor, 2)}
                for name, stock, vol in self.components() if name != "Template DNA"]