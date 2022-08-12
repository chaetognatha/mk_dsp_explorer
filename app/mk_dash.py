# Run this app with `python mk_dash.py` and
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
    Some things will always have to be read in at the start of the program
    we need to be able to update this state and perhaps this is the easiest
    way
    Returns - nothing, this function just restarts the app
    -------

    """
    sys.stdout.flush()
    os.execv(sys.argv[0], sys.argv)

def get_config():
    config_path = Path('conf.ini')
    conf = ConfigParser()
    # read the config file
    conf.read(config_path)
    illegal_chars = re.compile(r'[#@&\s,.\-+]')
    # DEFINE A SET OF FACTORS
    # Should default to a list of all columns in global df
    # has to start as empty list
    # can be replaced with factors from settings popup!
    factors = conf['Settings']['factors'].split(',')
    # also empty list, update after uploading data or in settings
    biomarkers = conf['Settings']['biomarkers'].split(',')
    # another really important global
    pat_id = conf['Settings']['patient_id']
    # controls need to be scrubbed from dataframe
    controls = conf['Settings']['controls'].split(',')
    factors = [re.sub(illegal_chars, '', f) for f in factors]
    biomarkers = [re.sub(illegal_chars, '', b) for b in biomarkers]
    pat_id = re.sub(illegal_chars, '', pat_id)
    controls = [re.sub(illegal_chars, '', c) for c in controls]
    biomarkers = [b for b in biomarkers if b not in controls]
    auto_encode = conf.getboolean('Settings', 'encode_factors')
    auto_transform = conf.getboolean('Settings', 'transform_biomarkers')
    encoding = conf['Settings']['encoding']
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
        print(f'get_config out: {d}')

    return d



def parse_contents(contents, filename, date, enc='utf-8'):
    illegal_chars = re.compile(r'[#@&\s,.\-+]')
    """

    Parameters
    ----------
    contents
    filename
    date
    enc: encoding of the file, default is utf-8 probably best untouched
    it seems dash cannot use anything but utf-8

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
        print('parse_contents: There was an error reading the uploaded file:\n',
              e)

    if dbg:
        print('parse_contents out:')
        print(f'df is a {type(df)}')
        print(df.info())
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
                        ]),
                        dbc.ModalFooter(
                            dbc.Button(
                                "OK", id="close", className="ms-auto",
                                n_clicks=0
                            )
                        ),
                    ],
                    id="modal",
                    is_open=True,
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
    dcc.Store(id='dataframe', storage_type='local'),
    # settings are initialized with conf.ini but are then
    # updated via the settings in the app
    # it is a dictionary
    dcc.Store(id='settings', storage_type='local'),
], style={'textAlign': 'center'})


# CALLBACKS
# MODAL CALLBACK ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
@app.callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    return not is_open if n1 or n2 else is_open

@app.callback(Output('settings', 'data'),
              Input('factors', 'value'),
              Input('biomarkers', 'value'),
              Input('pat-id', 'value'),
              Input('controls', 'value'),
              State('settings', 'data'),
              )
def update_config(f, b, p, c, d):
    """

    Parameters
    ----------
    f : factors
    b : biomarkers
    p : patient id
    c : controls
    d : settings data

    Returns dictionary that can be stored
    -------

    """
    if not d:
        # store preset values in a config file
        # be conservative of what to keep here!
        d = get_config()

    # check which changed
    if f:
        # update locally stored app settings
        d['factors'] = f
        if dbg:
            print(f'update_config: factors: {f}')
    if b:
        d['biomarkers'] = b
        if dbg:
            print(f'update_config: biomarkers: {b}')
    if p:
        d['patient_id'] = p
        if dbg:
            print(f'update_config: patient_id: {p}')
    if c:
        d['controls'] = c
        if dbg:
            print(f'update_config: controls: {c}')
    if dbg:
        print(f'update_config: settings: {PrettyDict(d)}')
    return d


@app.callback(Output('settings-out', 'children'),
              Input('dataframe', 'data'),
              State('settings', 'data'),
              prevent_initial_call=True,)
def init_settings(df, db):
    if df and db:
        if dbg:
            print(f'init_settings: df: {df.shape}')
            print(f'init_settings: db: {PrettyDict(db)}')
        # dataframe is stored as json
        df = pd.read_json(df)
        return html.Div([
            html.Div(['Select biomarkers']),
            dcc.Dropdown(list(df), db['biomarkers'],
                         id='biomarkers',
                         multi=True, persistence=True, persistence_type='local'),
            html.Div(['Select factor columns']),
            dcc.Dropdown(list(df), db['factors'], id='factors',
                         multi=True, persistence=True, persistence_type='local'),
            html.Div(['Select column with patient ID']),
            dcc.Dropdown(list(df), db['patient_id'], id='pat-id',
                         persistence=True, persistence_type='local'),
            html.Div(['Select controls']),
            dcc.Dropdown(db['biomarkers'] + db['controls'], db['controls'],
                         id='controls',
                         persistence=True, persistence_type='local', multi=True),

            # auto encode factors output
            html.Div(id='auto-encode-out'),
            # download data
            html.Button("Download CSV", id="btn_csv"),
            dcc.Download(id="download-dataframe-csv"),
        ])
    else:
        return html.Div(['Upload some data!'])

# xxxxxxxxxxxxxxxxxxxxxxxxxxx DATA UPLOAD xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
@app.callback(
              Output('dataframe', 'data'),
              Input('data_upload', 'contents'),
              State('data_upload', 'filename'),
              State('data_upload', 'last_modified'),
              State('dataframe', 'data'),
              State('settings', 'data'),
              prevent_initial_call=True,
              )
def update_df(list_of_contents, list_of_names, list_of_dates, df, db):
    if list_of_contents:
        # this is the first step in the event chain
        db = get_config()
        # empty cache
        for p in cache.glob('*'):
            p.unlink()
        # parse contents
        df = parse_contents(list_of_contents, list_of_names, list_of_dates,
                            db['encoding'])
        # encode factors
        if db['encode_factors']:
            if dbg:
                print('update_df: encoding factors')
            # just do all columns
            for factor in df.columns:
                if df[factor].dtype == 'object' and factor != db['patient_id']:
                    df[factor] = df[factor].astype('category')
        if db['transform_biomarkers']:
            if dbg:
                print('update_df: transforming biomarkers')
            try:
                df.loc[:, db['biomarkers']] = \
                    df.loc[:, db['biomarkers']].apply(lambda x: np.log2(x))
            except:
                # this just means config was not properly edited
                print('update_df: '
                      'You may have incorrectly '
                      'specified biomarkers in conf.ini')

        return df.to_json()




@app.callback(Output('auto-encode-out', 'children'),
              Input('settings', 'data'),
              Input('dataframe', 'data'),
              prevent_initial_call=True,)
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
    if df and db:
        offenders = []
        df = pd.read_json(df)
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
    df = pd.read_json(df)
    return dcc.send_data_frame(df.to_csv, "download.csv")


# DYNAMIC LAYOUT CALLBACK
@app.callback(Output('main', 'children'),
              Input('dataframe', 'data'),
              Input('settings', 'data'),
              prevent_initial_call=True)
def dynamic_layout(df, db):
    if df:
        return html.Div([
            html.Div(children=[
                # Parallel plot
                html.H2(children='Parallel plot'),
                # dbc layout grid
                dbc.Row([
                    dbc.Col(html.Div([
                                html.Label('Select all factors to show:'),
                                dcc.Dropdown(db['biomarkers'] + db['factors'],
                                db['factors'][0], id='dd-pp',
                                multi=True, persistence=True,
                                             persistence_type='local'),
                    ]), width="auto"),

                    dbc.Col(html.Div([
                            html.Label('Select which factor to use for color:'),
                            dcc.Dropdown(db['biomarkers'] + db['factors'],
                            db['factors'][0], id='color-pp', persistence=True,
                                         persistence_type='local'),
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
                                             persistence_type='local'),


                            ]), width=3
                        ),
                        dbc.Col(
                          html.Div([
                              html.Label('Select factor'),
                              # DROPDOWN 2
                              dcc.Dropdown(db['biomarkers'],
                                           db['biomarkers'][1],
                                           id='2', persistence=True,
                                           persistence_type='local'),
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
                                         persistence_type='local')
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
                                        max=25,
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
                                         persistence_type='local'),
                        ]), width=3
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
                                         persistence_type='local'),
                        ]), width=3
                    ),
                    dbc.Col(
                        html.Div([
                            # EVENT
                            html.Div(['Event']),
                            dcc.Dropdown(db['factors'], id='event',
                                         persistence=True, persistence_type='local'),
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
                                         persistence=True, persistence_type='local'),
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



######################## Paralell plot callback
@app.callback(Output('pp', 'figure'),
              Input('dd-pp', 'value'),
              Input('dataframe', 'data'),
              State('settings', 'data'),
              Input('color-pp', 'value'))
def parallel_plot(dd, df, db, c):
    fig = {}
    if df and len(dd) > 1 and c:
        df = pd.read_json(df)
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
            fig = px.parallel_coordinates(df, dimensions=dd,
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
    if df and value_1 and value_2:
        df = pd.read_json(df)
        all_together = f"{value_1}{value_2}{df}{db}"
        cached_path = cache.joinpath(cache_hash(all_together))
        if cached_path.exists():
            fig = pickle.load(open(cached_path, 'rb'))
        else:
            # the figure/plot created using the data filtered above
            r = stats.spearmanr(df[value_1], df[value_2], nan_policy='omit')
            txt = f"Spearman \u03C1:\t{round(r[0], 4)}\tP-value: {r[1]}"
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
    if df and factor:
        df = pd.read_json(df)
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
               Input('use-fdr', 'value')
               ])
def volcano(fixed_effect, plim, effects, df, db, fdr):
    fig = {}
    if df and fixed_effect:
        my_p = 'FDR' if fdr == 'FDR adj. p' else 'p-values'
        df = pd.read_json(df)
        # to avoid errors caused by getting inf in neglog10
        if plim == 1:
            plim = 0.9999
        all_together = f"{fixed_effect}{plim}{effects}{df}{db}{fdr}"
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
            fig.update_layout(showlegend=False, title='')
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
    if df and time and event and covariates:
        df = pd.read_json(df)
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
                    print(e)
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
                    print(e)
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
    app.run_server(debug=True, host='127.0.0.1', port='8050', proxy=None,
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
