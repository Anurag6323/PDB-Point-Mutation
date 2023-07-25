from flask import Flask, render_template, request, send_file
from point_mutator import point_mutator

# To aid in naming the download file
global pdb_id

# Initialising the App
app = Flask(__name__)


# Tool Page
@app.route('/', methods=["GET", "POST"])
def home():
    global pdb_id
    if request.method == "POST":

        # Getting input from the webpage
        pdb_id = request.form['pdb_name']
        residue = request.form['residue_num']
        req_mut = request.form['mutation']
        chain = request.form['chain_name']

        # Data Processing
        point_mutator(pdb_id, chain, residue, req_mut)

        # Rendering the Result
        return render_template("result.html")

    else:
        return render_template("tool.html")


# Download the mutated pdb file
@app.route('/download1')
def download_file1():
    path = "./downloads/mutated.pdb"
    return send_file(path, as_attachment=True, download_name=pdb_id+"_mutated.pdb")


# Download the pymol session file
@app.route('/download2')
def download_file2():
    path = "./downloads/mutated.pse"
    return send_file(path, as_attachment=True, download_name=pdb_id+"_mutated.pse")


# Download the interaction Summary file
@app.route('/download3')
def download_file3():
    path = "./downloads/interaction_summary/mutated_ligplot.sum"
    return send_file(path, as_attachment=True, download_name=pdb_id+"_mutated.sum")


# Running the App
if __name__ == '__main__':
    app.run(debug=True)
