# event_histogram

Generate one monthly histogram per annotation CSV in `event_data/` for manual DAS event identifications along a 50 km interrogated fiber.

## Setup

```bash
python -m pip install -r requirements.txt
```

## Usage

```bash
python event_histograms.py --input-dir event_data --output-dir figures --log-file figures/event_metrics.log --bin-km 5
```

Options:

- `--bin-km`: programmable X-axis spacing/chunk size in kilometers. The default is `5` km.
- `--total-length-km`: total fiber length. The default is `50` km.
- `--output-dir`: folder for the monthly histogram figures.
- `--log-file`: metrics report path.
- `--format`: figure format (`png`, `pdf`, or `svg`).

Each figure is an overlapped histogram with transparent bars for the legend classes `other`, `noise`, `vessel`, and `seismic`, plus a smooth distribution curve for each class. The X axis is fiber length in kilometers, and the Y axis is the number of events recorded.
