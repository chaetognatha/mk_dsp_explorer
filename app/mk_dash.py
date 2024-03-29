# Run this app with `python app/mk_dash.py` and
# visit http://127.0.0.1:8050/ in your web browser.
from spatialomics_toolkit import statistics
from statsmodels.stats.multitest import fdrcorrection as fdr
import numpy as np
import dash
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import dash_bio
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
from umap import UMAP
from lifelines import CoxPHFitter
import scipy.stats as stats
from pathlib import Path
from hashlib import md5
import pickle
import multiprocessing as mp
import base64
import datetime
import io
import re
from configparser import ConfigParser
import os
import sys
import json
import time

# GLOBALS

# debug
dbg = True
# plot div style
pltstyle = {'margin': 'auto', 'width': '80%'}
# buttonstyle
btnstyle = {'width': '40%'}
# sliderstyle
sldstyle = {'margin': 'auto', 'width': '80%'}

# relative path do cache directory
cache = Path('.cache/')
# boolean to clear cache
clear_cache = False

# make sure it exists
cache.mkdir(exist_ok=True)


# GLOBAL FUNCTIONS
def debug_log(string):
    """
    This function writes a string to the debug log
    :param string: string to write to log
    :return: nothing
    -------
    """
    dbgfile = Path('debug.log')
    with dbgfile.open('a') as f: 
        f.write(f'{time.asctime(time.localtime())}:\t{string}\n')
   # if dbg:
    #    debug_log(string)

def cache_hash(string):
    """
    This little function creates a unique identifier string
    :param string: string that consists of a combination of all ingoing
    parameters
    :return: unique reproducible identifier string
    """
    bstr = bytes(string, 'UTF-8')
    return md5(bstr).hexdigest()

class PrettyDict:
    """
    This class is used to make dictionaries more readable
    """
    def __init__(self, d):
        self.d = d

    def __repr__(self):
        for k, v in self.d.items():
            # single column
            if type(v) is list:
                my_str = f'{k}:\n' + '\n'.join(v)
            else:
                my_str = f'{k}:\n{v}\n'
        return my_str



def restart():
    """
    Not in use atm, but may be useful
    Returns - nothing, this function just restarts the app
    -------

    """
    sys.stdout.flush()
    os.execv(sys.argv[0], sys.argv)

def get_config():
    config_path = Path('conf.json')
    with open(config_path, 'r') as f:
        conf = json.load(f)
    # illegal characters
    illegal_chars = re.compile(r'[#@&\s,.\-+]')
    # DEFINE A SET OF FACTORS
    # Should default to a list of all columns in global df
    # has to start as empty list
    # can be replaced with factors from settings popup!
    factors = conf['factors']
    # also empty list, update after uploading data or in settings
    biomarkers = conf['biomarkers']
    # another really important global
    pat_id = conf['patient_id']
    # controls need to be scrubbed from dataframe
    controls = conf['controls']
    factors = [re.sub(illegal_chars, '', f) for f in factors]
    biomarkers = [re.sub(illegal_chars, '', b) for b in biomarkers]
    pat_id = re.sub(illegal_chars, '', pat_id)
    controls = [re.sub(illegal_chars, '', c) for c in controls]
    biomarkers = [b for b in biomarkers if b not in controls]
    auto_encode = conf['encode_factors']
    auto_transform = conf['transform_biomarkers']
    encoding = conf['encoding']
    # init state make dict from config file
    d = {
        'factors': factors,
        'biomarkers': biomarkers,
        'controls': controls,
        'patient_id': pat_id,
        'encode_factors': auto_encode,
        'transform_biomarkers': auto_transform,
        'encoding': encoding,
    }
    if dbg:
        debug_log(f'get_config out: {d}')

    return d

def set_config(db):
    config_file = Path('conf.json')
    with config_file.open('w') as f:
        json.dump(db, f, indent=4)
    if dbg:
        debug_log(f'set_config written to file: {db}')

def load_data():
    """
    This function reads in the data from the dataframe
    :return: dataframe
    -------
    """
    df = pd.DataFrame()
    # read in dataframe
    try:
        df = pd.read_csv('.cache/data.csv')
    except:
        debug_log('No dataframe found')
    if dbg:
        debug_log(f'get_data out: {type(df)}')
    return df

def save_data(df):
    """
    This function writes the dataframe to a csv file
    :param df: dataframe
    :return: nothing
    -------
    """
    df.to_csv('.cache/data.csv', index=False)
    if dbg:
        debug_log(f'set_data out: {type(df)}')


def parse_contents(contents, filename, date, enc='utf-8'):
    illegal_chars = re.compile(r'[#@&\s,.\-+]')
    """

    Parameters
    ----------
    contents
    filename
    date
    enc: encoding of the file, default is utf-8 probably best untouched
    it seems dash cannot use anything but utf-8 internally

    Returns: pandas dataframe
    -------

    """
    # we want to change the globals based on uploaded data
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)

    try:
        if filename.endswith('.csv'):
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode(enc, 'ignore')),
                encoding=enc)
            # CLEAN
            df.columns = df.columns.str.replace(illegal_chars, '',
                                                regex=True)

        elif filename.endswith('.tsv'):
            # Assume that the user uploaded a TSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode(enc, 'ignore')),
                delimiter='\t', encoding=enc)
            # CLEAN
            df.columns = df.columns.str.replace(illegal_chars, '',
                                                regex=True)

        elif filename.endswith('.xlsx'):
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
            df.columns = df.columns.str.replace(illegal_chars, '',
                                                regex=True)

        else:
            # try to identify the filetype
            fh = io.StringIO(decoded.decode(enc, 'ignore'))
            my_str = fh.read()
            result = {
                '\t' : 0,
                ',' : 0,
                ';' : 0
            }
            for separation_type in result.keys():
                result[separation_type] = len(my_str.split(separation_type))
            # find the most common separation character
            sep = max(result, key=lambda x: result[x])
            # read in the data
            df = pd.read_csv(io.StringIO(decoded.decode(enc, 'ignore')),
                             sep=sep, encoding=enc)
            # CLEAN
            df.columns = df.columns.str.replace(illegal_chars, '',
                                                regex=True)


    except Exception as e:
        debug_log('parse_contents: There was an error reading the uploaded file:\n',
              e)

    if dbg:
        debug_log('parse_contents out:')
        debug_log(f'df is a {type(df)}')
        debug_log(df.shape)
    return df


# DASH

app = Dash(__name__,
           external_stylesheets=[dbc.themes.MATERIA],
           suppress_callback_exceptions=True)


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& APP LAYOUT &&&&&&&&&&&&&&&&&&&&&&&&&
app.layout = html.Div(children=[
    # HEADER 1
    html.H1(children=' ', style={'textAlign': 'center'}),

    ################################### SETTINGS POPUP ########################
    html.Div([
        dcc.Loading([
            html.Div([
                dbc.Button("Settings", id="open", n_clicks=0),
                dbc.Modal(
                    [
                        dbc.ModalHeader(dbc.ModalTitle("Settings")),
                        dbc.ModalBody([

                            dcc.Upload(id='data_upload',
                                       children=html.Div([
                                           'Drag and Drop or ',
                                           html.A('Select a File')
                                       ],
                                           style={'width': '100%',
                                                  'height': '60px',
                                                  'lineHeight': '60px',
                                                  'borderWidth': '1px',
                                                  'borderStyle': 'dashed',
                                                  'borderRadius': '5px',
                                                  'textAlign': 'center'}),
                                       ),
                            html.Div([
                            ]),
                            
                            dcc.Loading(html.Div(id='settings-out')),
                            # auto encode factors output
                            html.Div(id='auto-encode-out'),
                        ]),
                        dbc.Row([
                            dbc.Col(
                            [html.Label('Transform biomarkers: '),
                            dcc.RadioItems(['Yes', 'No'], 'Yes', id='transform-biomarkers'),
                            html.Label('Encode factors: '),
                            dcc.RadioItems(['Yes', 'No'], 'Yes', id='encode-factors'),]
                            )
                        ]),

                        dbc.ModalFooter(
                            dbc.Row([
                                dbc.Col(
                                    html.Div([

                            dbc.Button("OK", id="close", className="ms-auto", n_clicks=0),
                            dbc.Button("Reload", id="reload", n_clicks=0),
                                                                            # download data
                            dbc.Button("Download CSV", id="btn_csv"),
                            dcc.Download(id="download-dataframe-csv"),
                                    ])
                        )]),
                        ),
                    ],
                    id="modal",
                    is_open=False,
                    size='lg',
                    fullscreen="lg-down",
                    centered=True,

                ),
                # here's how to format dbc elements
            ], className="d-grid gap-2 d-md-flex justify-content-md-end",
            ),
        ])
    ]),

    html.Hr(),
    html.Div(id='main'),
    # In this section we define variables that are stored and then shared
    # across callback functions
    # the advantage of this is that the app has a state that can be unique for
    # a user and is ready to be deployed in a situation where
    # several people use it at once
    # dataframe is the main dataframe, it is kept locally
    dcc.Store(id='dataframe', storage_type='session'),
    # settings are initialized with conf.ini but are then
    # updated via the settings in the app
    # it is a dictionary
    # settings moved to conf.json
    # keeping this to keep callback alive
    dcc.Store(id='settings', storage_type='session'),
], style={'textAlign': 'center'})


# CALLBACKS
# MODAL CALLBACK ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
@app.callback(
    Output("modal", "is_open"),
    Input("open", "n_clicks"), 
    Input("close", "n_clicks"),
    State("modal", "is_open"),
)
def toggle_modal(n1, n2, is_open):
    return not is_open if n1 or n2 else is_open

@app.callback(Output('settings', 'data'),
              Input("reload", "n_clicks"),
              Input('factors', 'value'),
              Input('biomarkers', 'value'),
              Input('pat-id', 'value'),
              Input('controls', 'value'),
              Input('transform-biomarkers', 'value'),
              Input('encode-factors', 'value'),
              State('dataframe', 'data'),
              State('settings', 'data'),
              #prevent_initial_call=True
              )
def update_config(n, factor, biomarker, patid, 
                    controls,transbio, encfac, df, db):
    """

    Parameters
    ----------
    f : factors
    b : biomarkers
    p : patient id
    c : controls
    d : settings data

    Modifies json file with settings
    also returns a json object with the settings but this is broken
    -------

    """
    df = load_data()
    if not db:
        # store preset values in a config file
        # be conservative of what to keep here!
        db = get_config()
        if dbg:
            debug_log('update_config: no db in memory calling get_config:')
            debug_log(db)
    else:
        try :

            db = json.loads(db)
            if dbg:
                debug_log('update_config: db from json:')
                debug_log(db)
        except Exception as e:
            debug_log('update_config: error loading json:', e)
            db = get_config()
            if dbg:
                debug_log(e)
                debug_log('update_config: no db in memory calling get_config:')
                debug_log(db)
    # check which changed
    if factor:
        for f in factor:
            if f not in df.columns:
                factor.remove(f)
                debug_log('update_config: factor not in dataframe:', f)
        # update locally stored app settings
        db['factors'] = factor
        if dbg:
            debug_log(f'update_config: factors: {factor}')
    if biomarker:
        for b in biomarker:
            if b not in df.columns:
                biomarker.remove(b)
                debug_log('update_config: biomarker not in dataframe:', b)
        db['biomarkers'] = biomarker
        if dbg:
            debug_log(f'update_config: biomarkers: {biomarker}')
    if patid:
        db['patient_id'] = patid
        if dbg:
            debug_log(f'update_config: patient_id: {patid}')
    if controls:
        if db['controls'] == 'None':
            db['controls'] = []
        else:
            db['controls'] = controls
        if dbg:
            debug_log(f'update_config: controls: {controls}')
    if transbio == 'Yes':
        db['transform_biomarkers'] = True
        if dbg:
            debug_log(f'update_config: transform_biomarkers: {transbio}')
    else:
        db['transform_biomarkers'] = False
        if dbg:
            debug_log(f'update_config: transform_biomarkers: {transbio}')
    if encfac == 'Yes':
        db['encode_factors'] = True
        if dbg:
            debug_log(f'update_config: encode_factors: {encfac}')
    else:
        db['encode_factors'] = False
        if dbg:
            debug_log(f'update_config: encode_factors: {encfac}')
    if dbg:
        debug_log(f'update_config (end): settings: {db}')

    # write to file
    set_config(db)

    return json.dumps(db)


@app.callback(Output('settings-out', 'children'),             
              Input("reload", "n_clicks"),
              Input('dataframe', 'data'),
              State('settings', 'data'),
)
def init_settings(df, n, db):
    db = get_config()
    df = load_data()
    if isinstance(df, pd.DataFrame) and db:
        # cleaning up the settings
        for factor in db['factors']:
            if factor not in df.columns:
                db['factors'].remove(factor)
                debug_log(f'init_settings: factor {factor} not in df')
        for biomarker in db['biomarkers']:
            if biomarker not in df.columns:
                db['biomarkers'].remove(biomarker)
                debug_log(f'init_settings: biomarker {biomarker} not in df')
        set_config(db)
        if dbg:
            debug_log(f'init_settings: df: {df.shape}')
            debug_log(f'init_settings: db: {db}')
        # dataframe is stored as json
        
        return html.Div([
            html.Div(['Select biomarkers']),
            dcc.Dropdown(list(df), db['biomarkers'],
                         id='biomarkers',
                         multi=True, persistence=True, persistence_type='memory'),
            html.Div(['Select factor columns']),
            dcc.Dropdown(list(df), db['factors'], id='factors',
                         multi=True, persistence=True, persistence_type='memory'),
            html.Div(['Select column with patient ID']),
            dcc.Dropdown(list(df), db['patient_id'], id='pat-id',
                         persistence=True, persistence_type='memory'),
            html.Div(['Select controls']),
            dcc.Dropdown(db['biomarkers'] + db['controls'], db['controls'],
                         id='controls',
                         persistence=True, persistence_type='memory', multi=True),


        ])
    else:
        if dbg:
            debug_log(f'init_settings: no df or db: {type(df)} {db}')
        return html.Div(['Upload some data!'])

# xxxxxxxxxxxxxxxxxxxxxxxxxxx DATA UPLOAD xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
@app.callback(
              Output('dataframe', 'data'),
              Input('data_upload', 'contents'),
              State('data_upload', 'filename'),
              State('data_upload', 'last_modified'),
              State('settings', 'data'),
              
              )
def update_df(list_of_contents, list_of_names, list_of_dates, db):
    if list_of_contents:
        # this is the first step in the event chain
        if dbg:
            debug_log('SENDING UPLOAD TO PARSE CONTENT')
        db = get_config()
        # empty cache
        for p in cache.glob('*'):
            p.unlink()
        # parse contents
        df = parse_contents(list_of_contents, list_of_names, list_of_dates,
                            db['encoding'])
        if dbg:
            debug_log('UPLOAD PARSED')
        # encode factors
        if db['encode_factors']:
            if dbg:
                debug_log('update_df: encoding factors')
            # just do all columns
            for factor in df.columns:
                if df[factor].dtype == 'object' and factor != db['patient_id']:
                    df[factor] = df[factor].astype('category')
        if db['transform_biomarkers']:
            if dbg:
                debug_log('update_df: transforming biomarkers')
            try:
                df.loc[:, db['biomarkers']] = \
                    df.loc[:, db['biomarkers']].apply(lambda x: np.log2(x))
                if dbg:
                    debug_log('TRANSFORM DONE')
            except:
                # this just means config was not properly edited
                debug_log('update_df: '
                      'You may have incorrectly '
                      'specified biomarkers in conf.ini')
        # store in cache
        save_data(df)

        return df.to_json()




@app.callback(Output('auto-encode-out', 'children'),
              Input('settings', 'data'),
              Input('dataframe', 'data'),
              )
def auto_encode(db, df):
    """
    this callback generates a key for factor encoding
    Parameters
    ----------
    db: database of settings
    df: dataframe from dcc.Store (json)

    Returns
    -------

    """
    db = get_config()
    df = load_data()
    if isinstance(df, pd.DataFrame) and db['encode_factors']:
                # cleaning up the settings
        for factor in db['factors']:
            if factor not in df.columns:
                db['factors'].remove(factor)
                debug_log(f'init_settings: factor {factor} not in df')
        if dbg:
            debug_log('auto_encode: encoding factors')
        offenders = []
        #df = pd.read_json(df)
        for factor in db['factors']:
            if df[factor].dtype == 'object' and factor != db['patient_id']:
                offenders.append(factor)
        if offenders:
            my_str = f"""
### Results of automatic factor encoding:
             
            """
            # auto-encode
            for factor in offenders:
                cats = pd.Categorical(df[factor])
                # this actually happens when the data is loaded
                #df[factor] = cats.codes
                expand = '\n\n'.join([f'{i}:\t{cat}' for i,
                                                    cat in enumerate(
                        cats.categories
                    )
                                      ])
                # using markdown double returns to make it look nice
                my_str += f"""
#####                          {factor}\n\n{expand}\n\n
"""

        else: my_str = "No categorical factors found."

        return dcc.Markdown(my_str)

#Download callback
@app.callback(
    Output("download-dataframe-csv", "data"),
    Input("btn_csv", "n_clicks"),
    State("dataframe", "data"),
    prevent_initial_call=True,
    
)
def func(n_clicks, df):
    # replace the stored dataframe with the cached one
    df = load_data()
    return dcc.send_data_frame(df.to_csv, "data.csv")

#lmm results download
@app.callback(
    Output('download-dataframe-lmm', 'data'),
    Input('btn_lmm', 'n_clicks')
)
def get_lmm_res(n):
    if n:
        lmm_res = pd.read_csv(cache.joinpath('lmm.csv'))
        return dcc.send_data_frame(lmm_res.to_csv, "lmm.csv")


# DYNAMIC LAYOUT CALLBACK
@app.callback(Output('main', 'children'),
              Input('dataframe', 'data'),
              Input('reload', 'n_clicks'),
              Input('settings', 'data'),
              )
def dynamic_layout(df, n, db):
    db = get_config()
    df = load_data()
    if dbg and isinstance(df, pd.DataFrame):
        debug_log(f'dynamic_layout: df: {df.shape}')
        debug_log(f'dynamic_layout: db: {db}')
    if isinstance(df, pd.DataFrame):
        # cleaning settings again
        for factor in db['factors']:
            if factor not in df.columns:
                db['factors'].remove(factor)
                debug_log(f'init_settings: factor {factor} not in df')
        for biomarker in db['biomarkers']:
            if biomarker not in df.columns:
                db['biomarkers'].remove(biomarker)
                debug_log(f'init_settings: biomarker {biomarker} not in df')
        return html.Div([
            html.Div(children=[
                # Parallel plot
                html.H2(children='Parallel plot'),
                # dbc layout grid
                dbc.Row([
                    dbc.Col(html.Div([
                                html.Label('Select all factors to show:'),
                                dcc.Dropdown(db['factors'],
                                db['factors'][0], id='dd-pp',
                                multi=True, persistence=True,
                                    persistence_type='memory'),
                    ]), width="auto"),

                    dbc.Col(html.Div([
                            html.Label('Select which factor to use for color:'),
                            dcc.Dropdown(db['biomarkers'] + db['factors'],
                            db['factors'][0], id='color-pp', persistence=True,
                                         persistence_type='memory'),
                    ]), width="auto"),
                    ]),


                dcc.Loading(dcc.Graph(id='pp')),
                html.Hr(),
            ]),
            html.Div(children=[
        # Scatterplot #########################################################
                html.Div([
                    # TITLE / TEXT
                    html.H2(children='Interactive scatterplot'),
                    dbc.Row([
                        dbc.Col(
                            html.Div(children=[

                                html.Label('Select factor'),
                                # DROPDOWN 1
                                dcc.Dropdown(db['biomarkers'],
                                             db['biomarkers'][0],
                                             id='1', persistence=True,
                                             persistence_type='memory'),


                            ]), width=3
                        ),
                        dbc.Col(
                          html.Div([
                              html.Label('Select factor'),
                              # DROPDOWN 2
                              dcc.Dropdown(db['biomarkers'],
                                           db['biomarkers'][1],
                                           id='2', persistence=True,
                                           persistence_type='memory'),
                          ]), width=3
                        ),
                    ]),


                    # PLOT
                    dcc.Loading(
                        html.Div([dcc.Graph(id='scatter-plot'), ],
                                 style=pltstyle)),

                ]),
                html.Hr(),
          # UMAP #############################################################
                html.Div([
                    html.H2(children='UMAP'),
                    dbc.Row([
                        dbc.Col(
                            html.Div([
                            html.Label('Select grouping factor:'),
                            dcc.Dropdown(db['factors'], id='umap_factor',
                                         persistence=True,
                                         persistence_type='memory')
                    ]), width=3)
                    ]),
                    ]),


                        html.Label('Select a minimum distance:'),
                        html.Div([
                            dcc.Slider(id='min_dist',
                                       min=0.1,
                                       max=1,
                                       step=0.1,
                                       value=0.1),
                        ], style=sldstyle),

                        html.Label('Select how many nearest neighbors:'),
                        html.Div([
                            dcc.Slider(id='n_neighbor',
                                        min=2,
            # max should be 1/4 * number of rows
                                        max=df.shape[0]//4,
                                        step=1,
                                        value=2), ], style=sldstyle) ]),

                    html.Div([dcc.Graph(id='umap'), ], style=pltstyle),


                html.Hr(),

            # VOLCANO PLOT #####################################################
            html.Div([
                html.H2(children='Linear Mixed Model Regression'),
                dbc.Row([
                    dbc.Col(
                        html.Div([
                            html.Label('Select Fixed Effect'),
                            dcc.Dropdown(db['factors'], id='fixed_effect',
                                         persistence=True,
                                         persistence_type='memory'),
                        ]), width=3
                    ),
                    dbc.Col(
                        html.Div(id='dyn-lmm-level')
                    ),
                    dbc.Col(
                        html.Div([
                            html.Label('Adjust p-value?'),
                            dcc.RadioItems(
                                ['Unadj. p', 'FDR adj. p'],
                                'FDR adj. p', id='use-fdr',
                                inline=False),
                        ]), width=3
                    ),
                    dbc.Col(
                        html.Div([
                            dbc.Button("Download Results", id="btn_lmm"),
                            dcc.Download(id="download-dataframe-lmm"),
                        ])
                    ),
                ]),
                html.Div([
                    html.Label('Select p-value limit'),
                    html.Div(
                        [dcc.Slider(0.001, 1, 0.05, value=0.05, id='lmm-plim',
                                    marks={0.05: '0.05*',
                                           0.5: '0.5',
                                           1: '1'
                                           }), ], style=sldstyle),
                ], style={'textAlign': 'center'}),

                'Effect sizes',
                html.Div([dcc.RangeSlider(
                    id='default-volcanoplot-input',
                    min=-3,
                    max=3,
                    step=0.05,
                    marks={i: {'label': str(i)} for i in range(-3, 3)},
                    value=[-1, 1]
                ), ], style=sldstyle),
                dcc.Loading(
                    html.Div([dcc.Graph(id='volcano_plot')], style=pltstyle)),
            ]),
            html.Hr(),
     # COX REGRESSION ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            html.Div([
                html.H2(children='Cox regression results'),
                html.Br(),

                dbc.Row([
                    dbc.Col(
                        html.Div([
                            # TIME
                            html.Div(['Time']),
                            dcc.Dropdown(db['factors'], id='time',
                                         persistence=True,
                                         persistence_type='memory'),
                        ]), width=3
                    ),
                    dbc.Col(
                        html.Div([
                            # EVENT
                            html.Div(['Event']),
                            dcc.Dropdown(db['factors'], id='event',
                                         persistence=True, persistence_type='memory'),
                        ]), width=3
                    ),
                    dbc.Col(
                        html.Div([
                            html.Label('Select covariants '
                                       'for univariate cox regression:'),
                            dcc.Dropdown(db['factors'] + db['biomarkers'],
                                         db['biomarkers'],
                                         id='covariates', multi=True,
                                         style=sldstyle,
                                         persistence=True, persistence_type='memory'),
                        ]), width=6
                    ),
                ], justify='evenly'),


                html.Br(),
                html.Label('Select p-value limit'),
                html.Div([dcc.Slider(0.001, 1, 0.05, value=0.05, id='cox-plim',
                                     marks={0.05: '0.05*',
                                            0.5: '0.5',
                                            1: '1'
                                            }), ], style=sldstyle),

                'Effect sizes',
                html.Div([dcc.RangeSlider(
                    id='cox-effects',
                    min=-3,
                    max=3,
                    step=0.05,
                    marks={i: {'label': str(i)} for i in range(-3, 3)},
                    value=[-0.1, 0.1],
                    allowCross=False
                ), ], style=sldstyle),
               dcc.Loading(
                   html.Div([dcc.Graph(id='cox_plot')], style=pltstyle)),

            ])
        ])
    else:
        return html.Div(['No DATA (check settings)'])

#DYNAMIC LMM LEVEL BINARIZER
@app.callback(Output('dyn-lmm-level', 'children'),
              Input('fixed_effect', 'value'))
def dyn_level_lmm(fe):
    df = load_data()
    db = get_config()
    if fe:
        this = [
                    html.Label('Test this/these level/s (1)'
                         ' against all others (0)'),
                    dcc.Dropdown(df[fe].unique(), 'None', id='bin-lvl',
                    persistence=False, multi=True),
                    ]
    else: this = 'Select a fixed effect'

    return this



######################## Parallel plot callback
@app.callback(Output('pp', 'figure'),
              Input('dd-pp', 'value'),
              Input('dataframe', 'data'),
              State('settings', 'data'),
              Input('color-pp', 'value'))
def parallel_plot(dd, df, db, c):
    fig = {}
    df = load_data()
    if isinstance(df, pd.DataFrame)and len(dd) > 1 and c:
        db = get_config()
        #df = pd.read_json(df)
        all_together = f"{dd}{df}{db}{c}"
        cached_path = cache.joinpath(cache_hash(all_together))
        if cached_path.exists():
            fig = pickle.load(open(cached_path, 'rb'))
        else:
            # the map here is not a lovely solution (it works), but I blame
            # plotly for not making this function behave like
            # literally all other plotting functions and just
            # map colors without complaining :(
            if df[c].dtype == object:
                df[c] = df[c].map({name: i for i,
                                               name in enumerate(
                    df[c].unique()
                )})
            for d in dd:
                if df[d].dtype == object:
                    df[d] = df[d].map({name: i for i,
                                                   name in enumerate(
                        df[d].unique()
                    )})
        # works well for more than 2 dimensions: px.colors.qualitative.Light24
            fig = px.parallel_categories(df, dimensions=dd,
                                          color=c,
                                          color_continuous_scale = \
                                              px.colors.diverging.Spectral,
                                          )

    return fig


# ~~~~~~~~~~~~~~~~~~~~~~ SCATTERPLOT CALLBACK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
@app.callback(Output(component_id='scatter-plot', component_property='figure'),
              [Input(component_id='1', component_property='value'),
               Input(component_id='2', component_property='value'),
               Input('dataframe', 'data'),
               State('settings', 'data')])
def scatterplot(value_1, value_2, df, db):
    fig = {}
    df = load_data()
    df = df.dropna(subset=[value_1, value_2])
    if isinstance(df, pd.DataFrame) and value_1 and value_2:
        db = get_config()
        #df = pd.read_json(df)
        all_together = f"{value_1}{value_2}{df}{db}"
        cached_path = cache.joinpath(cache_hash(all_together))
        if cached_path.exists():
            fig = pickle.load(open(cached_path, 'rb'))
        else:
            # the figure/plot created using the data filtered above
            r = stats.spearmanr(df[value_1], df[value_2], nan_policy='omit')
            txt = f"Spearman \u03C1:\t{round(r[0], 4)}\tp-value: {r[1]}"
            fig = px.scatter(df, x=value_1, y=value_2, title=txt,
                             # trendline='ols',
                             color=db['patient_id'],
                             color_continuous_scale = \
                                 px.colors.qualitative.Light24,
                             color_discrete_sequence = \
                                 px.colors.qualitative.Light24,
                             hover_name=db['patient_id'],
                             )

            pickle.dump(fig, open(cached_path, 'wb'))

    return fig


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   UMAP CALLBACK   ~~~~~~~~~~~~~~~~~~
@app.callback(Output('umap', 'figure'),
              Input('umap_factor', 'value'),
              Input('min_dist', 'value'),
              Input('n_neighbor', 'value'),
              Input('dataframe', 'data'),
              State('settings', 'data'))
def get_umap(factor, min_dist, n_neighbor, df, db):
    fig = {}
    db = get_config()
    df = load_data()
    df = df.dropna(subset=[factor])
    # in case biomarkers happen to be infinite
    df.loc[:,db['biomarkers']].replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=db['biomarkers'])
    if isinstance(df, pd.DataFrame) and factor:
        #df = pd.read_json(df)
        all_together = f"{factor}{min_dist}{n_neighbor}{df}{db}"
        cached_path = cache.joinpath(cache_hash(all_together))
        if cached_path.exists():
            proj_3d = pickle.load(open(cached_path, 'rb'))
        else:
            umap_3d = UMAP(n_components=3, init='random',
                           n_neighbors=n_neighbor, min_dist=min_dist)
            proj_3d = umap_3d.fit_transform(df[db['biomarkers']])
            pickle.dump(proj_3d, open(cached_path, 'wb'))

        fig = px.scatter_3d(
            proj_3d, x=0, y=1, z=2,
            color=df[factor], labels={'color': factor},
            color_continuous_scale=px.colors.diverging.Spectral,
            color_discrete_sequence=px.colors.qualitative.Alphabet,
            opacity=0.7
        )

        fig.update_traces(marker_size=6)

    return fig


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LMM REGRESSION PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~
@app.callback(Output('volcano_plot', 'figure'),
              [Input('fixed_effect', 'value'),
               Input('lmm-plim', 'value'),
               Input('default-volcanoplot-input', 'value'),
               Input('dataframe', 'data'),
               State('settings', 'data'),
               Input('use-fdr', 'value'),
               Input('bin-lvl', 'value'),
               ])
def volcano(fixed_effect, plim, effects, df, db, fdr, levels):
    fig = {}
    df = load_data()
    df = df.dropna(subset=[fixed_effect])
    my_title = f'{fixed_effect} as fixed effect'
    if isinstance(df, pd.DataFrame) and fixed_effect:
        db = get_config()
        if isinstance(levels,list):
            my_title += ' with '
            debug_log(type(levels))
            unique_levels = df[fixed_effect].unique()
            m = {}
            exclude = ''
            #expect levels to be a list
            for lvl in unique_levels:
                if lvl in levels:
                    m[lvl] = 1
                    my_title += f'{lvl} '
                else:
                    m[lvl] = 0
                    exclude += f'{lvl} '


            df[fixed_effect] = df[fixed_effect].map(m)
            my_title += f'vs {exclude}'
            debug_log('Recoded column temporarily:')
            debug_log(df[fixed_effect])
            debug_log('LMM mapping:')
            debug_log(m)
        my_p = 'FDR' if fdr == 'FDR adj. p' else 'p-values'
        #df = pd.read_json(df)
        # to avoid errors caused by getting inf in neglog10
        if plim == 1:
            plim = 0.9999
        all_together = f"{fixed_effect}{plim}{effects}{df}{db}{fdr}{levels}"
        cached_path = cache.joinpath(cache_hash(all_together))
        if cached_path.exists():
            fig = pickle.load(open(cached_path, 'rb'))
        else:
            # in this context we want a list of biomarker names
            biomarkers = db['biomarkers']
            # DEALING WITH NANS
            df.dropna(subset=biomarkers, inplace=True)
            df.dropna(subset=fixed_effect, inplace=True)
            results = statistics.lmm(df, biomarkers, FE=fixed_effect,
                                     RE=db['patient_id'])
            fig = dash_bio.VolcanoPlot(
                dataframe=results,
                point_size=8,
                effect_size_line=effects,
                effect_size_line_width=1,
                effect_size_line_color='black',
                genomewideline_width=1,
                genomewideline_value=-1 * np.log10(plim),
                genomewideline_color='red',
                effect_size='coefficients',
                p=my_p,
                snp=None,
                gene=None,
                annotation='biomarker',
                highlight_color='red',
                col='gray'
            )
            fig.update_layout(showlegend=False, title=my_title)
            fig.update_xaxes(title='LMM coefficient')
            # add annotation
            for i, p in enumerate(results[my_p]):
                if p <= plim and \
                        (abs(effects[0]) <= abs(results.coefficients[i]) or
                         results.coefficients[i] >= effects[1]):
                    fig.add_annotation(dict(font=dict(color='black', size=10),
                                            x=results.coefficients[i],
                                            y=-1 * np.log10(p),
                                            showarrow=False,
                                            text=results.biomarker[i],
                                            clicktoshow='onoff',
                                            yshift=6
                                            ))

                # resize

            pickle.dump(fig, open(cached_path, 'wb'))
            csv_path = cache.joinpath('lmm.csv')
            results.to_csv(csv_path,index=False)

    return fig


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  COX CALLBACK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
@app.callback(Output('cox_plot', 'figure'),
              [Input('covariates', 'value'),
               Input('time', 'value'),
               Input('event', 'value'),
               Input('cox-plim', 'value'),
               Input('cox-effects', 'value'),
               Input('dataframe', 'data'),
               State('settings', 'data'),])
def univariate_cox(covariates, time, event, plim, effects, df, db):
    fig = {}
    df = load_data()
    if isinstance(df, pd.DataFrame) and time and event and covariates:
        db = get_config()
        #df = pd.read_json(df)
        my_p = 'p-value'
        # generate hash
        all_ingoing = f"{covariates}{time}{event}{plim}{effects}{df}{db}"
        cached_path = cache.joinpath(cache_hash(all_ingoing))
        # someone might get the idea to input survival time or
        # event time as covariates, but this is not allowed
        covariates = [c for c in covariates if c != time and c != event]
        # load pre-generated data
        if cached_path.exists():
            results = pd.read_csv(cached_path)
        else:
            results = {}
            cluster = db['patient_id']
            # DEALING WITH NANS
            df.dropna(subset=[time, event], inplace=True)
            df.dropna(subset=covariates, inplace=True)
            for label in covariates:
                cph = CoxPHFitter()
                try:
                    cph.fit(df,
                            time,
                            event_col=event,
                            cluster_col=cluster,
                            formula=f'{label}'
                            )
                    summary = cph.summary.loc[:, ['exp(coef)', 'p']]
                except Exception as e:
                    debug_log(e)
                # in the current version we are taking e ^ coef
                # which can be turned into a percentual change in hazard
                # if female = 1, and the hazard(female)/hazard(male) =
                # 2 then females are twice as likely to die as males
                # or 2 - 1 = 100% more likely to die at any time t
                try:
                    results[label] = {'exp(coef)': float(summary.values[:, 0]),
                                      'hazard %': (float(summary.values[:, 0]) \
                                                   - 1) * 100,
                                      'p-value': float(summary.values[:, 1])}
                except Exception as e:
                    debug_log(e)
                    # sometimes the model does not converge
                    # setting results to 0 and p to 1
                    results[label] = {'exp(coef)': 0,
                                      'hazard %': 0,
                                      'p-value': 1}
            results = pd.DataFrame(results).transpose()
            results = results.rename_axis("biomarker").reset_index()
            #results['FDR'] = fdr(results['p-value'], alpha=0.05)
            results.to_csv(cached_path)

        fig = px.bar(results.loc[((results['exp(coef)'] >= effects[1]) |
                                  (results['exp(coef)'] <= effects[0])) &
                                 (results[my_p] <= plim)],
                     y='biomarker', x='hazard %',
                     hover_data=['p-value', 'exp(coef)', 'hazard %'],
                     color='p-value',
                     text_auto=True, orientation='h',
                     color_continuous_scale='sunset',
                     color_continuous_midpoint=0.5)

    return fig


if __name__ == '__main__':
    app.run_server(debug=False, host='127.0.0.1', port='8050', proxy=None,
                   dev_tools_ui=None,
                   dev_tools_props_check=None,
                   dev_tools_serve_dev_bundles=None,
                   dev_tools_hot_reload=None,
                   dev_tools_hot_reload_interval=None,
                   dev_tools_hot_reload_watch_interval=None,
                   dev_tools_hot_reload_max_retry=None,
                   dev_tools_silence_routes_logging=None,
                   dev_tools_prune_errors=None
                   )
