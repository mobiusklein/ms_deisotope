import io
import sys
import logging
import base64
try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
except ImportError:
    plt = None

from flask import Flask, request, jsonify, abort, Response, send_file

from ms_deisotope.data_source.proxi import S3MassSpectraDataArchive

app = Flask(__name__)

BLOCK_SIZE = 2 ** 20
logger = logging.getLogger("s3_proxi_server")


backend: S3MassSpectraDataArchive = None


def figax(*args, **kwargs):
    fig = Figure(*args, **kwargs)
    _canvas = FigureCanvasAgg(fig)
    return fig.add_subplot(1, 1, 1)


# Begin PROXI Services

@app.route(r"/api/proxi/v<int:version>/spectra")
def handle_usi(version):
    if version > 1:
        return abort(404)
    args = request.args
    usi = args['usi']
    app.logger.info("Loading USI %s", usi)
    response = backend.handle_usi(usi)
    return jsonify(**response)


@app.route(r"/api/proxi/v<int:version>/datasets/")
def handle_datasets(version):
    if version > 1:
        return abort(404)
    app.logger.info("Listing all %d datasets", len(backend.index))
    result = backend.enumerate_datasets()
    return jsonify(result)


@app.route(r"/api/proxi/v<int:version>/datasets/<identifier>")
def handle_describe_dataset(version, identifier):
    if version > 1:
        return abort(404)
    app.logger.info("Describing Dataset %r", identifier)
    result = backend.describe_dataset(identifier)
    return jsonify(result)


# Begin Non-PROXI services

@app.route(r"/api/datasets/")
def list_datasets():
    logger.info("Listing all datasets (%d)", len(backend.index))
    return jsonify(datasets=backend.list_datasets())


@app.route(r"/api/draw/<usi>")
def draw_usi(usi):
    scan = backend.handle_usi(usi, convert_json=False)
    ax = figax()
    ax.figure.set_dpi(160)
    ax.figure.set_figwidth(10)
    ax.figure.set_figheight(4)
    mz_start = request.args.get('start_mz')
    mz_end = request.args.get('end_mz')
    intensity_start = request.args.get("start_intensity")
    intensity_end = request.args.get("end_intensity")
    if scan.peak_set is None:
        scan.pick_peaks()
    scan.plot.centroided(ax=ax, lw=0.5, color='black')
    if scan.deconvoluted_peak_set is not None:
        scan.plot.deconvoluted(ax=ax, lw=0.5, color='teal')
    ax.set_title(usi)
    lo, hi = ax.get_xlim()
    if mz_start:
        lo = float(mz_start)
    if mz_end:
        hi = float(mz_end)
    ax.set_xlim(lo, hi)
    lo, hi = ax.get_ylim()
    if intensity_start:
        lo = float(intensity_start)
    if intensity_end:
        hi = float(intensity_end)
    ax.set_ylim(lo, hi)
    fig = ax.figure

    message = ''
    if scan.ms_level > 1:
        message = "Prec m/z: %0.3f\nPrec z: %r" % (
            scan.precursor_information.mz, scan.precursor_information.charge)

    fig.text(0.85, 0.85, "MS Level: %d\n%s" % (scan.ms_level, message), ha='left')
    buffer = io.BytesIO()
    ax.figure.savefig(buffer, format='png', bbox_inches='tight')
    return Response(b"data:image/png;base64," + base64.encodebytes(buffer.getvalue()), mimetype='image/png')


@app.route("/")
def home():
    return Response("""<!doctype html5>
<html>
    <head>
        <title>PROXI Server</title>
        <style>
        * {
            font-family: sans-serif
        }
        </style>
    </head>
    <body>
        <h1>PROXI Server</h1>
        <div>
            <select id="select-datasets" name="datasets">
            </select>
            <select id="select-datafiles" name="datafiles">
            </select>
            <input type="number" id="scan-number" name="scan-number" placeholder="Scan Number"/>
        </div>
        <div>
            <span>m/z Range</span>
            <input type="number" id="start-mz" name="start-mz" placeholder="Minimum m/z"/>
            <input type="number" id="end-mz" name="end-mz" placeholder="Maximum m/z"/>
        </div>
        <div>
            <span>Intensity Range</span>
            <input type="number" id="start-intensity" name="start-intensity" placeholder="Minimum Intensity"/>
            <input type="number" id="end-intensity" name="end-intensity" placeholder="Maximum Intensity"/>
        </div>
        <div>
            <img id="spectrum-image" />
        </div>
        <script>
            function init() {
                const datasetSelect = document.getElementById("select-datasets")
                const datafileSelect = document.getElementById("select-datafiles")
                const scanNumberInput = document.getElementById("scan-number")
                const spectrumImage = document.getElementById("spectrum-image")
                const startMzInput = document.getElementById("start-mz")
                const endMzInput = document.getElementById("end-mz")
                const startIntensityInput = document.getElementById("start-intensity")
                const endIntensityInput = document.getElementById("end-intensity")

                let datasets = null
                fetch("/api/datasets/").then((data) => data.json().then((j) => {
                    datasets = j.datasets
                    Object.keys(datasets).map((key) => {
                        datasetSelect.append(new Option(key, key, false, false))
                    })
                }))

                datasetSelect.onchange = function(event) {
                    datafileSelect.textContent = ''
                    if(datasetSelect.selectedOptions[0] === undefined) {
                        return
                    }
                    let nextDataset = datasetSelect.selectedOptions[0].value
                    datasets[nextDataset].map((item) => {
                        datafileSelect.append(new Option(item, item, false, false))
                    })
                }
                datasetSelect.onchange()

                const updatePlot = function(event) {
                    const scanNumber = scanNumberInput.value
                    const usi = `mzspec:${datasetSelect.selectedOptions[0].value}:${datafileSelect.selectedOptions[0].value}:scan:${scanNumber}`
                    console.log(usi)
                    fetch(`/api/draw/${usi}?start_mz=${startMzInput.value}&end_mz=${endMzInput.value}&start_intensity=${startIntensityInput.value}&end_intensity=${endIntensityInput.value}`).then((resp) => resp.text().then((data) => {
                        spectrumImage.src = data
                    }))
                }
                console.log("Registering Hooks")
                scanNumberInput.onblur = updatePlot

                startMzInput.onchange = updatePlot
                endMzInput.onchange = updatePlot
                startIntensityInput.onchange = updatePlot
                endIntensityInput.onchange = updatePlot
            }
            init()
        </script>
    </body>
</html>""")


def main(s3_bucket, host=None, port=None):
    if host is None:
        host = "0.0.0.0"
    if port is None:
        port = 8000
    logging.basicConfig(level="INFO", format=(
        '[%(asctime)s] %(levelname)s: %(message)s'))
    global backend
    logger.info("Building Backend For %r", s3_bucket)
    backend = S3MassSpectraDataArchive(s3_bucket)
    logger.info("Starting Server")
    app.run(host=host, port=port, debug=True, threaded=True)


if __name__ == "__main__":
    bucket = sys.argv[1]
    main(bucket)
