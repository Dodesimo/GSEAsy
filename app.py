import glob

import matplotlib
from flask import Flask, render_template, url_for, request, send_file
import os
from gsea_pipeline import initialize, read_annotation_data, meta_filter_and_adata_append, cell_filter, control_filter, \
    run_diffexp, process_diffexp, run_gseapy, ai_analysis, filter_control
from dotenv import load_dotenv
import markdown
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from flask import abort, render_template, request
from werkzeug.exceptions import HTTPException

matplotlib.use('Agg')

load_dotenv()

cell_subtypes = ""
control_group = ""
knockout_group = ""
experimental_description = ""
control_genes = ""

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.abspath("data")


@app.route("/")
def home():
    return render_template('input.html')

@app.route("/graphs")
def graphs():
    if os.path.exists("static/heatmap.png"):
        return render_template('graphs.html')
    else:
        initialize()
        adata, meta = read_annotation_data()

        if adata is None or meta is None:
            abort(400)

        meta_filter_and_adata_append(adata, meta)
        cell_subtype = cell_filter(adata, cell_subtypes)

        if cell_subtype is None:
            abort(400)

        enhanced_cell_subtype = control_filter(cell_subtype, control_group, knockout_group)

        if enhanced_cell_subtype is None:
            abort(400)

        results = run_diffexp(enhanced_cell_subtype, control_group, knockout_group)
        results = process_diffexp(results)

        if results is None:
            abort(400)

        out, table = run_gseapy(results)
        ai_analysis(out, cell_subtypes, experimental_description)

        # Render heat map and save it to disk
        gene_expression_data = adata.to_df()  # Assuming adata is an AnnData object
        experimental_or_control = adata.obs['group']  # Experimental/control group info
        gene_condition_data = gene_expression_data.groupby(experimental_or_control).mean()

        # Plot heat map
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(gene_condition_data.T[:10], cmap='viridis', vmin=0, vmax=1)
        plt.title('Heat Map of Gene Expression Across Experimental and Control Groups')

        # Save the heat map to disk
        heatmap_path = os.path.join("static", 'heatmap.png')
        plt.savefig(heatmap_path)
        plt.close(fig)

        #Filter, generate table for log fold.
        filtered = filter_control(control_genes)


        return render_template('graphs.html', table = table, filter_table = filtered)


@app.route('/show-text')
def show_text():
    with open(os.path.abspath('outputs/llmoutput.txt'), 'r') as file:
        text_content = file.read()
    return render_template('ai-results.html', text_content=markdown.markdown(text_content))


@app.route("/submit", methods=["POST", "GET"])
def submit():
    if request.method == "POST":
        files = glob.glob('static/*')
        if len(files) != 0:
            for f in files:
                os.remove(f)

        meta, counts = request.files['metadata-file'], request.files['counts-file']

        if meta is None or counts is None:
            abort(400)
        else:
            meta.save(os.path.join(app.config['UPLOAD_FOLDER'], meta.filename))
            counts.save(os.path.join(app.config['UPLOAD_FOLDER'], counts.filename))
            subtypes = pd.read_csv(os.path.join(app.config['UPLOAD_FOLDER'], meta.filename))[
                'Cell_subtype'].unique().tolist()
            groups = pd.read_csv(os.path.join(app.config['UPLOAD_FOLDER'], meta.filename))['group'].unique().tolist()

    return render_template('secondinput.html', subtypes=subtypes, groups=groups)


@app.route("/second_submit", methods=["POST", "GET"])
def second_submit():
    if request.method == "POST":
        global experimental_description, control_group, knockout_group, cell_subtypes, control_genes
        control_genes, experimental_description, control_group, knockout_group, cell_subtypes = request.form.get("cg"), request.form.get("ed"), request.form['control_group'], request.form['knockout_group'], request.form['cell-subtypes']

    return graphs()



#Error handlers
@app.errorhandler(400)
def page_not_found(error):
    return render_template('400.html'), 400
