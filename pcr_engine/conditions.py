"""Final reaction conditions in the tube (what the chemistry actually sees)."""

from dataclasses import dataclass


@dataclass
class ReactionConditions:
    reaction_volume_ul: float = 10.0
    primer_pmol: float = 5.0           # pmol of EACH primer in the tube
    mg_mM: float = 1.5                 # final MgCl2
    dntp_mM_each: float = 0.2          # final concentration of EACH dNTP
    buffer_x: float = 1.0              # final buffer strength (1.0 = 1X)
    taq_units: float = 0.25            # polymerase units in the tube

    @property
    def primer_uM(self) -> float:
        """Final primer concentration. pmol divided by microlitres gives micromolar."""
        return round(self.primer_pmol / self.reaction_volume_ul, 3)

    @property
    def primer_nM(self) -> float:
        return self.primer_uM * 1000

    @property
    def free_mg_mM(self) -> float:
        """dNTPs bind Mg2+ roughly one to one, so only the rest is free for the enzyme."""
        return round(self.mg_mM - 4 * self.dntp_mM_each, 2)

    def check(self) -> list[str]:
        """Problems that stop the reaction."""
        problems = []
        if self.reaction_volume_ul <= 0:
            problems.append("Reaction volume must be greater than zero.")
        if self.buffer_x <= 0:
            problems.append("No reaction buffer: the polymerase cannot work.")
        if self.taq_units <= 0:
            problems.append("No polymerase in the reaction.")
        if self.dntp_mM_each <= 0:
            problems.append("No dNTPs: nothing to build new DNA from.")
        if self.primer_pmol <= 0:
            problems.append("Primer amount must be greater than zero.")
        if self.free_mg_mM <= 0:
            problems.append("No free Mg2+ left after dNTPs bind it: polymerase is inactive.")
        return problems

    def notes(self) -> list[str]:
        """Warnings that do not stop the reaction but change the result."""
        notes = []
        if self.primer_uM > 1.0:
            notes.append(f"Primer at {self.primer_uM} uM is high: primer dimers become more likely.")
        elif self.primer_uM < 0.1:
            notes.append(f"Primer at {self.primer_uM} uM is low: yield may be poor.")
        if self.mg_mM > 3.0:
            notes.append("Mg2+ above 3 mM tends to reduce specificity.")
        if self.buffer_x > 0 and not 0.8 <= self.buffer_x <= 1.2:
            notes.append(f"Buffer is at {round(self.buffer_x, 2)}X instead of 1X.")
        units_per_ul = self.taq_units / self.reaction_volume_ul if self.reaction_volume_ul else 0
        if 0 < units_per_ul < 0.01:
            notes.append("Very little polymerase: yield may be weak.")
        elif units_per_ul > 0.1:
            notes.append("Excess polymerase can increase non-specific products.")
        return notes