import json
from pathlib import Path

from qecsim import app
from qecsim.models.color import Color666Code, Color666MPSDecoder
from qecsim.models.generic import BitFlipErrorModel


def run(distance: int, error_rate: float, nr_sims: int = 10000, results_dir: str = "./results_tn") -> None:
    code = Color666Code(distance)
    error_model = BitFlipErrorModel()
    decoder = Color666MPSDecoder(chi=8)
    data = app.run(code, error_model, decoder, error_rate, max_runs=nr_sims)
    filename = f"distance={distance},p={round(error_rate, 4)}.json"
    path = Path(results_dir)
    path.mkdir(parents=True, exist_ok=True)
    with (path / filename).open("w") as out:
        out.write(json.dumps(data))
