from flask import Flask, render_template, url_for, request, send_file
import os
from gsea_pipeline import initialize, read_annotation_data, meta_filter_and_adata_append, cell_filter, control_filter, run_diffexp, process_diffexp, run_gseapy, ai_analysis
from dotenv import load_dotenv
import markdown
import pandas as pd

load_dotenv()

groups = ""
cell_subtypes = ""

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.abspath("data")

@app.route("/")
def home():
    return render_template('input.html')


@app.route("/graphs")
def graphs():
    initialize()
    adata, meta = read_annotation_data()
    meta_filter_and_adata_append(adata, meta)
    cell_subtype = cell_filter(adata, cell_subtypes)
    enhanced_cell_subtype = control_filter(cell_subtype, groups)
    results = run_diffexp(enhanced_cell_subtype, groups)
    results = process_diffexp(results)
    out = run_gseapy(results)
    ai_analysis(out, cell_subtypes)

    return render_template('graphs.html')


@app.route('/show-text')
def show_text():
    with open(os.path.abspath('outputs/llmoutput.txt'), 'r') as file:
        text_content = file.read()
    return render_template('ai-results.html', text_content=text_content)


@app.route("/submit", methods=["POST", "GET"])
def submit():
    if request.method == "POST":
        meta, counts = request.files['metadata-file'], request.files['counts-file']
        meta.save(os.path.join(app.config['UPLOAD_FOLDER'], meta.filename))
        counts.save(os.path.join(app.config['UPLOAD_FOLDER'], counts.filename))
        subtypes = pd.read_csv(os.path.join(app.config['UPLOAD_FOLDER'], meta.filename))['Cell_subtype'].unique().tolist()
        groups = pd.read_csv(os.path.join(app.config['UPLOAD_FOLDER'], meta.filename))['group'].unique().tolist()

    return render_template('secondinput.html', subtypes=subtypes, groups=groups)

@app.route("/second_submit", methods=["POST", "GET"])
def second_submit():
    if request.method == "POST":
        global groups, cell_subtypes
        groups, cell_subtypes = request.form['groups'].split(".")[1], request.form['cell-subtypes']
        print(groups, cell_subtypes)

    return graphs()

