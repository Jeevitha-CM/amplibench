"""AmpliBench web API: serves the 3D bench and runs the PCR engine.

Start it with:   python -m uvicorn api:app --reload
Then open:       http://127.0.0.1:8000
"""

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse

from pcr_engine.recipe import Recipe
from pcr_engine.sequence_utils import parse_fasta
from pcr_engine.service import demo_payload, primer_info, run_experiment

app = FastAPI(title="AmpliBench API")
BENCH_PAGE = Path(__file__).parent / "bench3d" / "index.html"


@app.get("/")
def bench():
    """The 3D bench page."""
    return FileResponse(BENCH_PAGE)


@app.get("/api/demo")
def demo():
    """Demo templates, primer pairs and ladder names."""
    return demo_payload()


@app.post("/api/primers")
def primers(data: dict):
    """Tm and suggested Ta for each primer pair under the chosen recipe."""
    try:
        return primer_info(data["pairs"], Recipe(**data.get("recipe", {})))
    except (KeyError, TypeError, ValueError) as err:
        raise HTTPException(status_code=400, detail=str(err))


@app.post("/api/run")
def run(data: dict):
    """Run the experiment the bench describes and return results plus the gel image."""
    try:
        return run_experiment(data)
    except (KeyError, TypeError, ValueError) as err:
        raise HTTPException(status_code=400, detail=str(err))


@app.post("/api/fasta")
def fasta(data: dict):
    """Read FASTA text with the same parser and checks as the Streamlit app."""
    try:
        return {"records": parse_fasta(data.get("text", ""))}
    except ValueError as err:
        raise HTTPException(status_code=400, detail=str(err))