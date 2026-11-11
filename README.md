# AmpliBench

### Interactive In Silico PCR & Gel Electrophoresis

**A virtual molecular biology laboratory for PCR experimentation and DNA visualization.**

AmpliBench is an interactive scientific application designed to bring PCR experimentation and agarose gel electrophoresis into a computational environment.

Combining an interactive 3D laboratory, configurable PCR workflows, amplification simulation, and gel visualization, AmpliBench provides a virtual environment for exploring molecular biology experiments without consuming laboratory reagents.

---

## Live Application

**[Launch AmpliBench →](https://amplibench.onrender.com)**

Access the interactive virtual laboratory directly through your browser.

---

## Overview

AmpliBench integrates experimental configuration, computational simulation, and scientific visualization into a unified laboratory experience.

### Key Capabilities

| Module | Description |
|---|---|
| PCR Configuration | Prepare reaction setups, configure templates and primers, and evaluate reaction conditions. |
| Primer Analysis | Explore primer characteristics and reaction-dependent parameters. |
| Virtual Thermocycler | Configure and visualize conventional and gradient PCR workflows. |
| Amplification Simulation | Explore computational amplification outcomes and potential products. |
| Gel Electrophoresis | Simulate DNA migration and visualize electrophoretic band patterns. |
| 3D Workbench | Interact with a virtual laboratory through a guided experimental workflow. |
| Scientific Export | Export simulated products and gel visualizations for further analysis. |

---

## Interface Preview

### Virtual Laboratory

![AmpliBench 3D Workbench](docs/screenshots/bench_pipetting.png)

### Gel Visualization

![AmpliBench Simulated Gel](docs/screenshots/gel_result.png)

---

## Features

### Reaction Preparation
- FASTA template support, including multi-FASTA.
- Primer pair configuration.
- Flexible reaction preparation workflows.
- Master mix and manual reagent configuration.
- Automated reaction volume and concentration calculations.
- Reaction composition validation.

### PCR Simulation
- Computational exploration of primer-dependent amplification.
- Visualization of potential amplification products.
- Reaction controls and simulated experimental outcomes.
- Product sequence generation and export.

### Thermocycling
- Programmable PCR cycling conditions.
- Conventional and gradient PCR workflows.
- Interactive thermocycler visualization.
- Configurable experimental parameters.

### Gel Electrophoresis
- Customizable agarose gel conditions.
- DNA ladder selection and configuration.
- Sample loading and electrophoresis simulation.
- Visualization of DNA fragment migration.
- Adjustable run conditions and gel rendering.

### Interactive 3D Workbench
- Virtual pipetting and reaction preparation.
- PCR block loading and thermocycler visualization.
- Well-by-well gel loading.
- Electrophoresis and gel imaging.
- Experimental modification and result re-visualization.

### Export & Visualization
- Gel images in PNG and SVG formats.
- Simulated PCR products in FASTA format.

---

## Getting Started

### Online Access

**[Open AmpliBench →](https://amplibench.onrender.com)**

### Local Development

**Requirements:** Python 3.10 or newer.

```bash
git clone https://github.com/Jeevitha-CM/amplibench.git
cd amplibench

python -m pip install -r requirements.txt
```

Launch the application:

```bash
python -m uvicorn api:app --reload
```

The development server enables local access for testing and development.

### Testing

Install development dependencies:

```bash
python -m pip install -r requirements-dev.txt
```

Execute the test suite:

```bash
python -m pytest
```

---

## Technology Stack

| Technology | Application |
|---|---|
| Python | Scientific computing and application logic |
| FastAPI | Backend framework |
| Uvicorn | ASGI server |
| Gunicorn | Production process management |
| Three.js | Interactive 3D visualization |
| Biopython | Biological sequence processing |
| Matplotlib | Scientific visualization |
| pandas | Data handling |
| pytest | Automated testing |

---

## Scientific Scope

AmpliBench is designed for molecular biology education, experimental planning, and computational exploration.

The application provides simulated representations of PCR amplification and agarose gel electrophoresis. These simulations are intended to support conceptual understanding and experimental preparation rather than replace laboratory validation.

### Limitations

Simulated outcomes are not guaranteed to reproduce experimental results. Actual amplification and electrophoretic behaviour depend on biological, biochemical, and experimental conditions that may not be fully represented computationally.

**Experimental validation remains essential.**

---

## Developer

**Jeevitha C M**  
PhD Research Scholar  
SRM Institute of Science and Technology

AmpliBench is developed as part of the **21-Day Scientific Coding Challenge**, exploring the intersection of scientific computing, molecular biology, and interactive software development.

---

## License

Distributed under the MIT License.

See [`LICENSE`](LICENSE) for details.

---

<p align="center">
  <strong>AmpliBench</strong>
  <br>
  An interactive approach to computational molecular biology.
</p>
