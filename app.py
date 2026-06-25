from flask import Flask, request, render_template, jsonify, Response
import os
import subprocess
import threading
import queue
import csv
import io
import json
import glob

app = Flask(__name__)

BASE_DIR = os.getcwd()
PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(PIPELINE_DIR, "CONFIG")

os.makedirs(CONFIG_DIR, exist_ok=True)

# Store for job status and results
job_store = {}


@app.route("/get_analysis/<name>")
def get_analysis(name):
    path = os.path.join(CONFIG_DIR, name)
    if not os.path.exists(path):
        return jsonify({"error": "Analysis not found"}), 404

    fasta_path = os.path.join(path, "seq.fasta")
    params_path = os.path.join(path, "params.make")

    fasta_content = ""
    if os.path.exists(fasta_path):
        with open(fasta_path) as f:
            fasta_content = f.read()[:800]

    params = {}
    if os.path.exists(params_path):
        with open(params_path) as f:
            for line in f:
                line = line.strip()
                if "=" in line:
                    key, val = line.split("=", 1)
                    params[key.strip()] = val.strip()

    return jsonify({"fasta": fasta_content, "params": params})


def find_mut_file(analysis_name, fastq_basename):
    search_dirs = [
        os.path.join(PIPELINE_DIR, "30-pipeline", "blast", "030-mut"),
        os.path.join(PIPELINE_DIR, "30-pipeline"),
        PIPELINE_DIR,
    ]
    for d in search_dirs:
        if not os.path.exists(d):
            continue
        for pattern in [f"*{fastq_basename}*.mut", f"*{analysis_name}*.mut", "*.mut"]:
            matches = glob.glob(os.path.join(d, pattern))
            if matches:
                return max(matches, key=os.path.getmtime)
    return None


def parse_mut_file(filepath):
    """
    Parse .mut file which has two sections separated by blank lines:
    1. Single mutations  (Gene, Prot Mut, Mut Reads, WT Reads, VAF)
    2. Compound mutations (Gene, Prot Mut 1, Prot Mut 2, Double Mutated Transcripts, Total Transcripts)
    Returns {"single": [...], "compound": [...]}
    """
    result = {"single": [], "compound": []}
    if not filepath or not os.path.exists(filepath):
        return result

    with open(filepath) as f:
        content = f.read()

    # Split into blocks separated by one or more blank lines
    blocks = []
    current = []
    for line in content.splitlines():
        if line.strip() == "":
            if current:
                blocks.append(current)
                current = []
        else:
            current.append(line)
    if current:
        blocks.append(current)

    for idx, block in enumerate(blocks):
        if not block:
            continue
        header = [h.strip() for h in block[0].split("\t")]
        if all(h == "" for h in header):
            continue

        # Detect compound section by header content
        is_compound = any("prot mut 1" in h.lower() or "prot mut 2" in h.lower()
                          or "double mutated" in h.lower() or "total transcripts" in h.lower()
                          for h in header)

        rows = []
        for line in block[1:]:
            parts = [p.strip() for p in line.split("\t")]
            if all(p == "" for p in parts):
                continue
            while len(parts) < len(header):
                parts.append("")
            row = dict(zip(header, parts[:len(header)]))
            if any(v for v in row.values()):
                rows.append(row)

        if is_compound:
            result["compound"].extend(rows)
        else:
            result["single"].extend(rows)

    return result


def stream_subprocess(cmd, cwd, q, label):
    """Run subprocess and stream output line-by-line to queue."""
    def emit(msg, level="info"):
        q.put(json.dumps({"type": "log", "level": level, "msg": msg}))
    try:
        proc = subprocess.Popen(
            cmd, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1
        )
        for line in proc.stdout:
            line = line.rstrip()
            if line:
                emit(f"  {line}")
        proc.wait()
        return proc.returncode
    except Exception as e:
        emit(f"❌ Erro ao executar {label}: {str(e)}", "error")
        return -1


def run_pipeline(job_id, is_new, analysis_name, fastq_r1_filename,
                 config_analysis_path, params=None):
    q = job_store[job_id]["queue"]

    def emit(msg, level="info"):
        q.put(json.dumps({"type": "log", "level": level, "msg": msg}))

    try:
        if is_new:
            emit("⚙️ Criando params.make...")
            params_path = os.path.join(config_analysis_path, "params.make")
            with open(params_path, "w") as f:
                f.write(f"MAX_E_VALUE={params['evalue']}\n")
                f.write(f"MIN_ALIGN_LEN={params['min_len']}\n")
                f.write(f"MIN_PERC_IDENT={params['identity']}\n")
                f.write(f"MIN_DEPTH={params['depth']}\n")
                f.write(f"MIN_VAF={params['vaf']}\n")

            emit("🗄️ Construindo banco de dados BLAST...")
            rc = stream_subprocess(
                ["bash", "database.bash", analysis_name],
                PIPELINE_DIR, q, "database.bash"
            )
            if rc != 0:
                emit(f"❌ Erro no database.bash (código {rc})", "error")
                emit("💡 Verifique se o arquivo FASTA está correto e não está vazio.", "warn")
                q.put(json.dumps({"type": "done", "success": False,
                                   "error": "Erro no database.bash"}))
                return

        emit("🔄 Convertendo FASTQ → FASTA...")
        emit("🔗 Realizando joining das reads pareadas...")
        emit("💥 Executando BLAST (pode demorar alguns minutos)...")
        emit("🔍 Identificando mutações e calculando VAF...")

        rc2 = stream_subprocess(
            ["bash", "mutation.bash", f"00-fastq/{fastq_r1_filename}", analysis_name],
            PIPELINE_DIR, q, "mutation.bash"
        )

        if rc2 != 0:
            emit(f"❌ Erro no mutation.bash (código {rc2})", "error")
            q.put(json.dumps({"type": "done", "success": False,
                               "error": "Erro no mutation.bash"}))
            return

        emit("✅ Pipeline concluído! Buscando resultados...")

        fastq_basename = (fastq_r1_filename
                          .replace(".fastq.gz", "")
                          .replace("_R1_001", "")
                          .replace("_R1_", ""))
        mut_path = find_mut_file(analysis_name, fastq_basename)

        parsed = {"single": [], "compound": []}
        if mut_path:
            emit(f"📄 Arquivo de resultado: {os.path.basename(mut_path)}")
            parsed = parse_mut_file(mut_path)
            total = len(parsed["single"]) + len(parsed["compound"])
            n_single   = len(parsed["single"])
            n_compound = len(parsed["compound"])
            emit(f"🎯 {n_single} mutação(ões) simples, {n_compound} composta(s) encontrada(s)", "success")
        else:
            emit("⚠️ Arquivo .mut não encontrado automaticamente.", "warn")

        job_store[job_id]["result"] = parsed
        job_store[job_id]["mut_path"] = mut_path

        # Send all rows flattened for the frontend renderTable (it splits by header)
        all_rows = parsed["single"] + parsed["compound"]
        q.put(json.dumps({"type": "done", "success": True,
                           "rows": all_rows,
                           "single": parsed["single"],
                           "compound": parsed["compound"],
                           "count": len(parsed["single"])}))

    except Exception as e:
        emit(f"❌ Erro inesperado: {str(e)}", "error")
        q.put(json.dumps({"type": "done", "success": False, "error": str(e)}))


@app.route("/run", methods=["POST"])
def run():
    import uuid, re
    fastq_r1 = request.files.get("fastq_r1")
    fastq_r2 = request.files.get("fastq_r2")  # optional

    if not fastq_r1 or not fastq_r1.filename:
        return jsonify({"error": "Arquivo FASTQ R1 não enviado"}), 400

    create_new = request.form.get("create_new") == "true"

    if create_new:
        analysis_name = request.form.get("analysis_new", "").strip()
        fasta = request.files.get("fasta")
        if not analysis_name:
            return jsonify({"error": "Nome da análise vazio"}), 400
        if not re.match(r'^[A-Za-z0-9_\-]+$', analysis_name):
            return jsonify({"error": "Nome da análise deve conter apenas letras, números, hífens (-) e underscores (_). Evite espaços e caracteres especiais."}), 400
        if not fasta:
            return jsonify({"error": "FASTA obrigatório para nova análise"}), 400
        params = {
            "evalue": request.form.get("evalue", "1e-5"),
            "min_len": request.form.get("min_len", "49"),
            "identity": request.form.get("identity", "90"),
            "depth": request.form.get("depth", "500"),
            "vaf": request.form.get("vaf", "3"),
        }
    else:
        analysis_name = request.form.get("analysis_existing", "").strip()
        if not analysis_name:
            return jsonify({"error": "Nenhuma análise selecionada"}), 400
        fasta = None
        params = None

    # Validate R1 has the expected pattern
    if "_R1_" not in fastq_r1.filename:
        return jsonify({"error": "Arquivo R1 não segue o padrão esperado (_R1_)"}), 400

    # If R2 provided, validate that it matches R1
    has_r2 = fastq_r2 and fastq_r2.filename and fastq_r2.filename != ""
    if has_r2:
        if "_R2_" not in fastq_r2.filename:
            return jsonify({"error": "Arquivo R2 não segue o padrão esperado (_R2_)"}), 400
        prefix_r1 = fastq_r1.filename.replace(".fastq.gz", "")
        prefix_r2 = fastq_r2.filename.replace(".fastq.gz", "").replace("_R2_", "_R1_")
        if prefix_r1 != prefix_r2:
            return jsonify({"error": "R1 e R2 não correspondem à mesma amostra"}), 400

    fastq_dir = os.path.join(PIPELINE_DIR, "00-fastq")
    os.makedirs(fastq_dir, exist_ok=True)
    fastq_r1.save(os.path.join(fastq_dir, fastq_r1.filename))
    if has_r2:
        fastq_r2.save(os.path.join(fastq_dir, fastq_r2.filename))

    config_analysis_path = os.path.join(CONFIG_DIR, analysis_name)

    if create_new:
        os.makedirs(config_analysis_path, exist_ok=True)
        fasta_file_path = os.path.join(config_analysis_path, "seq.fasta")
        fasta.save(fasta_file_path)
    else:
        if not os.path.exists(config_analysis_path):
            return jsonify({"error": "Análise não encontrada"}), 400

    job_id = str(uuid.uuid4())
    job_store[job_id] = {"queue": queue.Queue(), "result": None, "mut_path": None}

    thread = threading.Thread(
        target=run_pipeline,
        args=(job_id, create_new, analysis_name, fastq_r1.filename,
              config_analysis_path, params),
        daemon=True
    )
    thread.start()

    return jsonify({"job_id": job_id})


@app.route("/status/<job_id>")
def status(job_id):
    if job_id not in job_store:
        return Response(
            'data: {"type":"error","msg":"Job não encontrado"}\n\n',
            mimetype="text/event-stream"
        )

    def generate():
        q = job_store[job_id]["queue"]
        while True:
            try:
                msg = q.get(timeout=60)
                yield f"data: {msg}\n\n"
                data = json.loads(msg)
                if data.get("type") == "done":
                    break
            except queue.Empty:
                yield 'data: {"type":"ping"}\n\n'

    return Response(
        generate(),
        mimetype="text/event-stream",
        headers={"X-Accel-Buffering": "no", "Cache-Control": "no-cache"}
    )


@app.route("/export/<job_id>")
def export_csv(job_id):
    if job_id not in job_store:
        return "Job não encontrado", 404
    parsed = job_store[job_id].get("result")
    if not parsed:
        return "Nenhum dado para exportar", 404

    # Support both old flat-list format and new dict format
    if isinstance(parsed, list):
        single, compound = parsed, []
    else:
        single   = parsed.get("single", [])
        compound = parsed.get("compound", [])

    if not single and not compound:
        return "Nenhum dado para exportar", 404

    output = io.StringIO()
    # "sep=;" tells Excel (especially in Brazilian locale) to use semicolon as delimiter
    output.write("sep=;\n")

    if single:
        writer = csv.DictWriter(output, fieldnames=list(single[0].keys()),
                                extrasaction='ignore', delimiter=';')
        writer.writeheader()
        writer.writerows(single)

    if compound:
        if single:
            output.write("\n")   # blank separator line between sections
        writer2 = csv.DictWriter(output, fieldnames=list(compound[0].keys()),
                                 extrasaction='ignore', delimiter=';')
        writer2.writeheader()
        writer2.writerows(compound)

    output.seek(0)
    return Response(
        output.getvalue(),
        mimetype="text/csv",
        headers={"Content-Disposition": "attachment; filename=mutations.csv"}
    )


@app.route("/", methods=["GET"])
def index():
    analyses = [
        d for d in os.listdir(CONFIG_DIR)
        if os.path.isdir(os.path.join(CONFIG_DIR, d))
    ]
    return render_template("index.html", analyses=analyses)


if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', threaded=True)
