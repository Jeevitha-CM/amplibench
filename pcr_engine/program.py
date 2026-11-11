"""Thermocycler program: the settings you enter on the PCR machine."""

from dataclasses import dataclass

GRADIENT_WELLS = 11   # gradient positions on our machine


def gradient_temps(low: float, high: float, columns: int = GRADIENT_WELLS) -> list[float]:
    """Starting gradient: temperatures evenly spaced from low to high (editable per well)."""
    if high < low:
        low, high = high, low
    step = (high - low) / (columns - 1)
    return [round(low + i * step, 1) for i in range(columns)]


@dataclass
class ThermocyclerProgram:
    initial_denat_temp: float = 95.0
    initial_denat_sec: float = 180
    denat_temp: float = 95.0
    denat_sec: float = 30
    anneal_temp: float = 58.0
    anneal_sec: float = 30
    extension_temp: float = 72.0
    extension_sec: float = 60
    cycles: int = 30
    final_extension_sec: float = 300
    hold_temp: float = 4.0            # after the run, held indefinitely

    def check(self) -> list[str]:
        """Return problems that would stop amplification."""
        problems = []
        if self.denat_temp < 90:
            problems.append("Denaturation below 90 C: DNA strands do not separate.")
        if not 65 <= self.extension_temp <= 75:
            problems.append("Extension temperature outside Taq's working range (65-75 C).")
        if self.anneal_temp >= self.denat_temp:
            problems.append("Annealing temperature must be lower than denaturation temperature.")
        if self.cycles < 1:
            problems.append("At least one cycle is needed.")
        if self.anneal_sec <= 0 or self.extension_sec <= 0 or self.denat_sec <= 0:
            problems.append("Step times must be greater than zero.")
        return problems

    def total_minutes(self) -> float:
        """Total run time, ignoring ramp time between temperatures."""
        per_cycle = self.denat_sec + self.anneal_sec + self.extension_sec
        total = self.initial_denat_sec + self.cycles * per_cycle + self.final_extension_sec
        return round(total / 60, 1)