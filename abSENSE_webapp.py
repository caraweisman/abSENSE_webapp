from flask import Flask, render_template, request, flash, redirect, url_for, session
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from Custom_Plot_Homolog_Detectability import *
from Precomp_Plot_Homolog_Detectability import *
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_svg import FigureCanvasSVG
import matplotlib.pyplot as plt
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from flask_wtf import Form
from wtforms import TextField, IntegerField, StringField, BooleanField, validators, SubmitField, FieldList, DecimalField, FloatField, RadioField, SelectField

app = Flask(__name__)
app.secret_key = '\x0f\xba\x04\xfa\xd0\x95\x9co\xf0\t\xbd-\x90\xf7t\xa9\x19\r\xb6\x0f(s\xf0,'

hcnumspec = 5

class numspecs_ask(Form):
    numspecs = IntegerField("How many species do you have (including your focal species)?", [validators.NumberRange(min=3, message = 'An integer number of at least 3 species is required!'), validators.DataRequired(message='Required!')])
    submit = SubmitField("Continue")
    
class fungi_gene_ask(Form):
    gene_desc = TextField('S. cerevisiae gene name, ID, description, NCBI protein accession (eg NP_XXXX)?')
    submit = SubmitField("Search")

class insects_gene_ask(Form):
    gene_desc = TextField('D. melanogaster gene name, ID, description, NCBI protein accession (eg NP_XXXX)?')
    submit = SubmitField("Search")

#### home
    
@app.route('/')
def home():
    
    return render_template('home.html')

@app.route('/about')
def about():
    return render_template('about.html')

#### precomputed fungi
    
@app.route('/browse_fungi', methods = ['GET', 'POST'])
def browse_fungi_search():
        form = fungi_gene_ask()
        if request.method == 'GET':
            return render_template('browse_fungi.html', form = form)
        elif request.method == 'POST':
            session['fungiquery'] = form.gene_desc.data
            #return form.gene_desc.data
            return redirect('/browse_fungi_searchresults')
        
@app.route('/browse_fungi_searchresults', methods = ['GET', 'POST'])
def browse_fungi_select():
    acc_descs = np.genfromtxt('Fungi_Data/Scer_accessionids_descriptions', dtype=str, delimiter='\n')
    choices = []
    for i in range(0, len(acc_descs)):
        if session['fungiquery'].lower() in acc_descs[i].lower():
            choices.append((acc_descs[i],acc_descs[i]))
    if len(choices) == 0:
        choices.append(('No results found!', 'No results found!'))
    # add back button
    class fungi_gene_confirm(Form):
        fungi_gene_choices = RadioField('Did you mean:', choices = choices)
        submit = SubmitField("Select gene and see results")
    form = fungi_gene_confirm()
    if request.method == 'GET':
        return render_template('browse_fungi_searchresults.html', form = form)
    elif request.method == 'POST':
        session['geneid'] = re.split('\s', form.fungi_gene_choices.data)[0][1:]
        session['fullgenedesc'] = form.fungi_gene_choices.data[1:]
        return redirect('/browse_fungi_generesults')  
    
@app.route('/browse_fungi_generesults', methods = ['GET', 'POST'])
def browse_fungi_generesults():
    bitscores = np.genfromtxt('Fungi_Data/Fungi_Bitscores', dtype=str, delimiter='\t')
    pred_specs = []
    species = bitscores[0][1:]
    for i in range(0, len(bitscores)):
        if session['geneid'] in bitscores[i][0]:
            for j in range(0, len(bitscores[i][1:])):
                if bitscores[i][j+1] != 'N/A' and bitscores[i][j+1] != '0':
                    pred_specs.append(species[j]) # is this indexing right?
    if len(pred_specs) < 3:
        return render_template('not_enough_points.html')
    else:
        fig, aparam, bparam, rsq, notinfitspecs, mlpreds, highnns, lownns, undets, ambigspecs, absentspecs, species, scores = fungi_make_pred_plot(session['geneid'], pred_specs, 'fungi')
        show_bitscores = []
        for i in range(0, len(scores)):
            if scores[i] == 'N/A':
                show_bitscores.append('Orthology ambiguous')
            elif scores[i] == '0':
                show_bitscores.append("No homolog detected")
            else:
                show_bitscores.append(scores[i])
        figfile = BytesIO()
        plt.tight_layout()
        plt.savefig(figfile, format='png')
        figfile.seek(0)
        figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
        plt.close()
        return render_template('fungi_output.html', fullname = session['fullgenedesc'], plot=figdata_png, species = species, preds = len(absentspecs), notinfitspecs = notinfitspecs, absents = absentspecs, mlpreds = mlpreds, highnns = highnns, lownns = lownns, undets = undets, numscores = len(show_bitscores), bitscores=show_bitscores)

@app.route('/browse_fungi_random', methods = ['GET', 'POST'])
def fungi_show_random():
    bitscores = np.genfromtxt('Fungi_Data/Fungi_Bitscores', dtype=str, delimiter='\t')
    species = bitscores[0][1:]
    pred_specs = []
    showconds = False
    while showconds == False: 
        pred_specs = []
        randgenenum = np.random.randint(1, len(bitscores))
        gene = bitscores[randgenenum][0]
        #print gene
        numabs = 0
        for j in range(0, len(bitscores[randgenenum][1:])):
            if bitscores[randgenenum][j+1] != 'N/A' and bitscores[randgenenum][j+1] != '0':
                pred_specs.append(species[j]) # is this indexing right?
            elif bitscores[randgenenum][j] == '0':
                numabs = numabs + 1
        if len(pred_specs) > 3: #and numabs > 0:
            showconds = True
    acc_descs = np.genfromtxt('Fungi_Data/Scer_accessionids_descriptions', dtype=str, delimiter='\n')
    for i in range(0, len(acc_descs)):
        if gene.lower() in acc_descs[i].lower():
           session['fullgenedesc'] = acc_descs[i][1:]
    fig, aparam, bparam, rsq, notinfitspecs, mlpreds, highnns, lownns, undets, ambigspecs, absentspecs, species, scores = fungi_make_pred_plot(gene, pred_specs, 'fungi')
    show_bitscores = []
    for i in range(0, len(scores)):
        if scores[i] == 'N/A':
            show_bitscores.append('Orthology ambiguous')
        elif scores[i] == '0':
            show_bitscores.append("No homolog detected")
        else:
            show_bitscores.append(scores[i])
    figfile = BytesIO()
    plt.tight_layout()
    plt.savefig(figfile, format='png')
    figfile.seek(0)
    figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
    plt.close()
    return render_template('fungi_output_random.html', fullname = session['fullgenedesc'], plot=figdata_png, species = species, preds = len(absentspecs), notinfitspecs = notinfitspecs, absents = absentspecs, mlpreds = mlpreds, highnns = highnns, lownns = lownns, undets = undets, numscores = len(show_bitscores), bitscores=show_bitscores)

#### precomputed insects
    
    
@app.route('/browse_insects', methods = ['GET', 'POST'])
def browse_insects_search():
        form = insects_gene_ask()
        if request.method == 'GET':
            return render_template('browse_insects.html', form = form)
        elif request.method == 'POST':
            session['insectsquery'] = form.gene_desc.data
            #return form.gene_desc.data
            return redirect('/browse_insects_searchresults')
        
@app.route('/browse_insects_searchresults', methods = ['GET', 'POST'])
def browse_insects_select():
    acc_descs = np.genfromtxt('Insect_Data/Dmel_accessionids_descriptions', dtype=str, delimiter='\n')
    choices = []
    for i in range(0, len(acc_descs)):
        if session['insectsquery'].lower() in acc_descs[i].lower():
            choices.append((acc_descs[i],acc_descs[i]))
    if len(choices) == 0:
        choices.append(('No results found!', 'No results found!'))
    #print choices
    # add back button
    class insects_gene_confirm(Form):
        insects_gene_choices = RadioField('Did you mean:', choices = choices)
        submit = SubmitField("Select gene and see results")
    form = insects_gene_confirm()
    if request.method == 'GET':
        return render_template('browse_insects_searchresults.html', form = form)
    elif request.method == 'POST':
        session['geneid'] = re.split('\s', form.insects_gene_choices.data)[0][1:]
        session['fullgenedesc'] = form.insects_gene_choices.data[1:]
        return redirect('/browse_insects_generesults')  
    
    
@app.route('/browse_insects_generesults', methods = ['GET', 'POST'])
def browse_insects_generesults():
    bitscores = np.genfromtxt('Insect_Data/Insect_Bitscores', dtype=str, delimiter='\t')
    pred_specs = []
    species = bitscores[0][1:]
    for i in range(0, len(bitscores)):
        if session['geneid'] in bitscores[i][0]:
            for j in range(0, len(bitscores[i][1:])):
                if bitscores[i][j+1] != 'N/A' and bitscores[i][j+1] != '0':
                    pred_specs.append(species[j]) # is this indexing right?
                    #print bitscores[i][j], species[j]
    #print len(pred_specs), pred_specs
    if len(pred_specs) < 3:
        return render_template('not_enough_points.html')
    else:
        fig, aparam, bparam, rsq, notinfitspecs, mlpreds, highnns, lownns, undets, ambigspecs, absentspecs, species, scores = fungi_make_pred_plot(session['geneid'], pred_specs, 'insects')
        show_bitscores = []
        for i in range(0, len(scores)):
            if scores[i] == 'N/A':
                show_bitscores.append('Orthology ambiguous')
            elif scores[i] == '0':
                show_bitscores.append("No homolog detected")
            else:
                show_bitscores.append(scores[i])
        figfile = BytesIO()
        plt.tight_layout()
        plt.savefig(figfile, format='png')
        figfile.seek(0)
        figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
        plt.close()
        return render_template('insects_output.html',  fullname = session['fullgenedesc'], plot=figdata_png, species = species, preds = len(absentspecs), notinfitspecs = notinfitspecs, absents = absentspecs, mlpreds = mlpreds, highnns = highnns, lownns = lownns, undets = undets, numscores = len(show_bitscores), bitscores=show_bitscores)


@app.route('/browse_insects_random', methods = ['GET', 'POST'])
def insects_show_random():
    bitscores = np.genfromtxt('Insect_Data/Insect_Bitscores', dtype=str, delimiter='\t')
    species = bitscores[0][1:]
    pred_specs = []
    showconds = False
    while showconds == False: 
        pred_specs = []
        randgenenum = np.random.randint(1, len(bitscores))
        gene = bitscores[randgenenum][0]
        #print gene
        numabs = 0
        for j in range(0, len(bitscores[randgenenum][1:])):
            if bitscores[randgenenum][j+1] != 'N/A' and bitscores[randgenenum][j+1] != '0':
                pred_specs.append(species[j]) # is this indexing right?
            elif bitscores[randgenenum][j] == '0':
                numabs = numabs + 1
        if len(pred_specs) > 3: #and numabs > 0:
            showconds = True
    acc_descs = np.genfromtxt('Insect_Data/Dmel_accessionids_descriptions', dtype=str, delimiter='\n')
    for i in range(0, len(acc_descs)):
        if gene.lower() in acc_descs[i].lower():
           session['fullgenedesc'] = acc_descs[i][1:]
    fig, aparam, bparam, rsq, notinfitspecs, mlpreds, highnns, lownns, undets, ambigspecs, absentspecs, species, scores = fungi_make_pred_plot(gene, pred_specs, 'insects')
    show_bitscores = []
    for i in range(0, len(scores)):
        if scores[i] == 'N/A':
            show_bitscores.append('Orthology ambiguous')
        elif scores[i] == '0':
            show_bitscores.append("No homolog detected")
        else:
            show_bitscores.append(scores[i])
    figfile = BytesIO()
    plt.tight_layout()
    plt.savefig(figfile, format='png')
    figfile.seek(0)
    figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
    plt.close()
    return render_template('insects_output_random.html', fullname = session['fullgenedesc'], plot=figdata_png, species = species, preds = len(absentspecs), notinfitspecs = notinfitspecs, absents = absentspecs, mlpreds = mlpreds, highnns = highnns, lownns = lownns, undets = undets, numscores = len(show_bitscores), bitscores=show_bitscores)


#### custom

@app.route('/custom_numspec', methods = ['GET', 'POST'])
def form_test():
    form = numspecs_ask()
    if request.method == 'POST':
        if form.validate() == False:
            flash('All fields are required.')
            return render_template('custom_numspec.html', form=form)
        else:
            session['numspecs'] = form.numspecs.data
            numspecs = form.numspecs.data
            return redirect('/custom_enterspec') 
    elif request.method == 'GET':
         return render_template('custom_numspec.html', form = form)

@app.route('/custom_enterspec',  methods = ['GET', 'POST'])
def custom_enterspec2():
    class specs_dists_ask(Form):
        genename = TextField("Gene name:")
        genesize = IntegerField('Gene length (aa):', default=400)
        ecutoff = FloatField('E-value cutoff:', default=0.001)
        focalspec = TextField("Name of focal species:")
        #predspecfocal = SelectField('Predict bitscore in this species <br> (using species for which "no" is selected)?', choices=[('no','no'),('yes','yes')])
        bsfocal = FloatField('Bitscore in focal species <br> (score of gene vs itself; must be larger than 0):',  [validators.DataRequired(),validators.NumberRange(min=1)])
        dbsize = FloatField('Database size (per species) (aa):', default=4000000)
        otherspecs = FieldList(TextField('Name of additional species: <br> ', [validators.DataRequired()]), min_entries=(session['numspecs']-1) , max_entries=(session['numspecs']-1))
        otherspecdists = FieldList(FloatField('Distance from focal species: <br> ', [validators.DataRequired()]), min_entries=(session['numspecs']-1) , max_entries=(session['numspecs']-1))
        #predspecsothers = FieldList(SelectField('Predict bitscore in this species <br> (using species for which "no" is selected)?', choices=[('no','no'),('yes','yes')]), min_entries=(session['numspecs']-1) , max_entries=(session['numspecs']-1))
        bsothers = FieldList(FloatField('Bitscore in this species <font color="fc8123"> (Enter 0 if none detected): </font> <br>'),  min_entries=(session['numspecs']-1) , max_entries=(session['numspecs']-1))
        submit = SubmitField("Analyze")
    form = specs_dists_ask()
    if request.method == 'POST':
        genename = form.genename.data
        session['genename'] = form.genename.data
        genesize = float(form.genesize.data)
        ethresh = float(form.ecutoff.data)
        dbsize = float(form.dbsize.data)
        specnames = []
        predspeclocs = []

        predspeclocs.append(0) # focal spec is in species list first (done below)
        for i in range(0, len(form.bsothers.data)):
            if form.bsothers.data[i] != 0 :
                predspeclocs.append(i+1)
            
        specnames.append(str(form.focalspec.data))
        for i in range(0, len(form.otherspecs.data)):
            specnames.append(str(form.otherspecs.data[i]))

                
        distances = []
        distances.append(float(0))
        for i in range(0, len(form.otherspecdists.data)):
            distances.append(float(form.otherspecdists.data[i]))
            
        bitscores = []
        bitscores.append(float(form.bsfocal.data))
        for i in range(0, len(form.bsothers.data)):
            bitscores.append(float(form.bsothers.data[i]))
            
        if form.validate() == False:
            return render_template('custom_enterspec.html', form=form)
        else:
            
            fig, aparam, bparam, rsq, notinfitspecs, mlpreds, highnns, lownns, undets, ambigspecs, absentspecs, species, scores = custom_make_pred_plot(genename, genesize, ethresh, specnames, distances, dbsize, predspeclocs, bitscores)
            show_bitscores = []

            for i in range(0, len(scores)):
                if scores[i] == 'N/A':
                    show_bitscores.append('Orthology ambiguous')
                elif scores[i] == 0 :
                    show_bitscores.append("No homolog detected")
                else:
                    show_bitscores.append(str(scores[i]))
            figfile = BytesIO()
            plt.tight_layout()
            plt.savefig(figfile, format='png')
            figfile.seek(0)
            figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
            plt.close()
            return render_template('custom_output.html', fullname = genename, plot=figdata_png, species = species, preds = len(absentspecs), notinfitspecs = notinfitspecs, absents = absentspecs, mlpreds = mlpreds, highnns = highnns, lownns = lownns, undets = undets, numscores = len(show_bitscores), bitscores=show_bitscores)      
    elif request.method == 'GET':
        a = session['numspecs']
        return render_template('custom_enterspec.html', form=form)

