\# AmpliBench



\*\*An interactive \*in silico\* PCR and gel electrophoresis bench\*\*



\*\*Try it:\*\* \_live demo link coming soon\_



AmpliBench lets you plan a PCR and see the likely result on an agarose gel before you spend reagents in the lab. Add your template and primers, set up the reactions, run them on a virtual thermocycler and load the products on a simulated gel. It shows where primers bind, which products form, when off-target bands or primer dimers appear, and what the gel would look like.



It has two views of the same experiment:



\- \*\*3D Lab\*\*: a guided virtual bench. Pipette each tube, load the thermocycler, run the PCR, then load, run and image the gel. Made for learning and teaching.

\- \*\*Classic Workbench\*\*: the same experiment as clear, numbered forms and tables, with no animation. Made for quick planning and checking primers.



> AmpliBench is a teaching and planning tool. Its results are rule-based simulations of primer binding and gel migration, not guarantees of wet-lab outcomes. Always confirm with a real experiment.



!\[3D Lab](docs/screenshots/bench\_pipetting.png)

!\[Simulated gel](docs/screenshots/gel\_result.png)



\## How to use



1\. \*\*Open AmpliBench\*\* and choose \*\*3D Lab\*\* or \*\*Classic Workbench\*\*. You can switch at any time; your experiment stays the same.

2\. \*\*Add templates\*\*: upload or paste FASTA, or press \*\*Load demo data\*\* to explore with a ready-made example.

3\. \*\*Enter primer pairs\*\* (5' to 3'). AmpliBench shows each primer's Tm and a suggested annealing temperature.

4\. \*\*Set up the reactions\*\*: 2X master mix or individual components, reaction volume, primer amount (pmol) and template amount (ng). Water is calculated for you.

5\. \*\*Fill the tubes\*\* (up to 11). Leave out the template in one tube for a no-template control (NTC).

6\. \*\*Program the thermocycler\*\*: conventional, or gradient to test several annealing temperatures at once.

7\. \*\*Run the PCR and the gel\*\*: choose agarose %, run time and ladder. Change gel settings afterwards and the gel redraws instantly.

8\. \*\*Keep your work\*\*: download the gel image or product sequences, print the \*\*Protocol\*\* sheet for the bench, or \*\*Save\*\* the experiment and \*\*Open\*\* it later.



The 3D Lab opens with a short tour; press \*\*Tour\*\* in the header to see it again. \*\*Undo\*\*, \*\*Redo\*\* and \*\*Clear all\*\* are always in the header.



\## Features



\- \*\*Templates and primers\*\*: upload or paste FASTA (multi-FASTA supported), enter primer pairs.

\- \*\*Reaction setup\*\*: 2X master mix or manual components (10X buffer, MgCl2, dNTPs, Taq). Primers in pmol, template in ng. Water and final concentrations are calculated, with checks for overfilled tubes, missing buffer, polymerase or free Mg2+.

\- \*\*Primer analysis\*\*: nearest-neighbour Tm under the actual reaction conditions, suggested annealing temperature, 3' dimer risk.

\- \*\*Thermocycler\*\*: 11 tube holders; staged program (initial denaturation, cycling, final extension, hold), conventional or \*\*gradient PCR\*\* (each holder gets its own editable annealing temperature; moving a tube into an occupied holder swaps the two).

\- \*\*PCR simulation\*\* per tube:

&#x20; - binding with a strict 3'-end rule and mismatch tolerance that depends on how far the annealing temperature is below Tm (off-target products),

&#x20; - product size limited by extension time,

&#x20; - product sequence built from the primers (primer-introduced changes are carried into the product),

&#x20; - \*\*primer dimers\*\* from 3'-end complementarity, stronger at higher primer concentration and in no-template controls,

&#x20; - band brightness from template copy number and cycle number (logistic amplification with a plateau).

\- \*\*Pipetting the way labs do it\*\*: each component into each tube, or a \*\*cocktail per primer pair\*\* (master mix + primers + water for all tubes of that gene plus one extra), then template or NTC water per tube. Separate tubes for every primer and template; adding the wrong one gives the wrong result, as at the bench.

\- \*\*Gel electrophoresis\*\*: agarose percentage, run time (short or long run), comb with wells matching your samples, ladder on one or both sides, load volume. Small fragments can run off a long run.

\- \*\*Realistic gel appearance\*\*: longer products glow brighter (stain binds along the length), a faint cloud of unused primers at the bottom of each lane, and a smear when the reaction is overamplified (too much template or too many cycles).

\- \*\*Ladders\*\*: NEB 100 bp, NEB 1 kb Plus (with per-band DNA mass), GeneRuler 50 bp, 100 bp and 1 kb, generic ladders, or a custom ladder.

\- \*\*3D bench\*\*: a short guided tour on first use, click the machines themselves to work with them, pipette each tube with fresh tips, load the block, watch the run, load the gel well by well, run it and image it. Go back to any step and change it; re-view the gel at other agarose percentages or run times instantly.

\- \*\*Undo / Redo / Clear all\*\* for every change in either view (Ctrl+Z / Ctrl+Y).

\- \*\*Exports\*\*: gel as PNG or SVG, products as FASTA, a printable \*\*protocol sheet\*\* (primers with Tm, per-tube and cocktail volumes, tube plan, program, gel loading order, simulated result), and \*\*save / open\*\* of the whole experiment as a JSON file.



\## Scientific model and limitations



| Part | How it is modelled | Limitation |

|---|---|---|

| Tm | Biopython nearest-neighbour Tm with primer, Na+, Mg2+ and dNTP concentrations | Mismatched sites use the perfect-match Tm |

| Binding | Exact match of the last 3 bases at the 3' end; mismatches allowed elsewhere, more at lower annealing temperature | Rule-based thresholds, not thermodynamics; no secondary structure or GC-clamp scoring |

| Amplification | Logistic growth from template copies to a plateau | Illustrative brightness, not a yield prediction |

| Primer dimer | Perfect 3'-end complementarity of at least 4 bases (3 at high primer concentration) | Does not model hairpins or partial complementarity |

| Template | Linear | Products spanning the origin of circular plasmids are not detected |

| Gel brightness | Product copies x size factor (square root of length); primer cloud; smear after many cycles at plateau | Illustrative, not quantitative |

| Gel | Speed falls with fragment size (sigmoid in log size), distance = speed x time | Not calibrated to a specific gel box, buffer or stain |

| Extension | About 1 kb per minute (Taq) | High-fidelity or fast enzymes differ |



Templates can be up to 200 kb each (genes, plasmids, cosmids, BACs). Whole genomes would need an indexed search, as used by BLAST or UCSC In-Silico PCR.



\## Roadmap



\- Primer quality check before the run (length, GC%, 3' clamp, Tm difference, hairpins, self-dimers)

\- Mg2+ concentration affecting specificity

\- Touchdown PCR programs

\- Circular (plasmid) templates

\- Multiplex PCR (more than one primer pair per tube)

\- Indexed primer search for genome-scale templates



\## Data



The demo template is \*\*synthetic\*\* (random sequence with deliberate edits). It is not a biological sequence. Ladder band sizes and masses come from manufacturer documentation (below); check the datasheet of your lot.



\## References



\- Cock PJA et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. \*Bioinformatics\* 25(11):1422-1423.

\- Allawi HT, SantaLucia J Jr (1997). Thermodynamics and NMR of internal G-T mismatches in DNA. \*Biochemistry\* 36(34):10581-10594 (Biopython's default nearest-neighbour table).

\- von Ahsen N, Wittwer CT, Schutz E (2001). Oligonucleotide melting temperatures under PCR conditions. \*Clinical Chemistry\* 47(11):1956-1961 (Mg2+ and dNTP correction used by Biopython).

\- New England Biolabs: 100 bp DNA Ladder (N3231) and 1 kb Plus DNA Ladder (N3200) product information.

\- Thermo Scientific: GeneRuler 50 bp (SM0371), 100 bp (SM0241/SM0243) and 1 kb (SM0311) DNA Ladder product information.

\- Three.js (https://threejs.org), FastAPI (https://fastapi.tiangolo.com), Matplotlib.



\## Author



Jeevitha C M, PhD Research Scholar, SRM Institute of Science and Technology. 



\## License



MIT License. See \[LICENSE](LICENSE).



\---



Running AmpliBench yourself or contributing? See the \[developer guide](docs/DEVELOPMENT.md).

