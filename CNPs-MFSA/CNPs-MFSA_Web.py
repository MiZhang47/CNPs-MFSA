"""
Last updated: 2025/05
Complex natural products (CNPs) Structure Annotation Tool constructed based on a modular fragmentation-based structural assembly (MFSA) strategy.
"""

# =========================
#  Import Python libraries
# =========================
# Import libraries for data processing and file handling
import io  # Used for handling file streams
import os  # Provides a way of using operating system-dependent functionalities
import re  # Regular expression operations for string matching
import chardet
import itertools # Generates top-k for neutral loss tab
import base64  # Encoding and decoding operations for binary data
import zipfile  # Handles zip file creation and extraction
import pandas as pd  # Data manipulation and analysis library
import numpy as np  # Numerical computation library
from typing import Optional, Tuple, List, Dict  # Type hinting for function signatures
from itertools import combinations
from collections import defaultdict

# Import external libraries for chemical analysis and handling mass spectrometry data
import sqlite3  # Provides SQLite database operations
from pyteomics import mass, mgf  # Pyteomics library for handling mass spectrometry files and calculations
from matchms import Spectrum, calculate_scores  # MatchMS library for handling mass spectrometry spectra and similarity calculations
from matchms.similarity import CosineGreedy  # Cosine similarity scoring algorithm from MatchMS
from rdkit import Chem  # RDKit library for chemical informatics and computational chemistry
from rdkit.Chem import Draw  # RDKit's drawing utilities for visualizing molecular structures
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Descriptors import MolWt

# Import libraries for building and managing a web-based dashboard using Dash
import dash  # Dash framework for building web applications
from dash import dcc, html, dash_table  # Dash core components, HTML layout, and tables
import dash_bootstrap_components as dbc  # Bootstrap components for Dash
from dash.dependencies import Input, Output, State, MATCH, ALL  # Dash dependencies for interactive components
import plotly.express as px  # Plotly Express for creating visualizations
from flask import send_file  # Flask utility for sending files in response to HTTP requests
from dash import dcc, dash_table
import plotly.graph_objects as go

# =========================
#  Initialize the Dash app
# =========================
app = dash.Dash(__name__,
    title="CNPs-MFSA", # Set the title of the web application
)  # Creates a new Dash application instance

# Define a directory for uploading and storing temporary files
UPLOAD_DIR = 'temp_uploads'  # Directory name for temporary file storage
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR)  # Create the directory if it doesn't exist

# Define global dictionaries to store results from various processing steps
global_results_recognize = {}  # Stores results from 'Target Extract' processing
global_results_composition = {}  # Stores results from 'Composition Calculate' processing
global_results_adduct = {}  # Stores results from 'Adduct Ion Calculate' processing
global_results_neutralloss = {}  # Stores results from 'Neutral Loss Extract' processing
global_results_featureion = {}  # Stores results from 'Feature Ion Extract' processing
global_results_cosine = {}  # Stores results from 'Cosine Score Calculate' processing
global_results_annotate = {}  # Stores results from 'Annotation' processing
global_results_visualization = {}  # Stores results from 'Visualization' processing

# =========================
#  Dash Layout
# =========================
app.layout = html.Div([  # Main container for the application
    # Add a description text at the top of all tabs
    html.Div([
        # Add the logo image on the left
        html.Div([
            # Add the logo image on the left
            html.Img(src='/assets/favicon.ico', style={'height': '60px', 'marginRight': '10px', 'marginLeft': '20px',
                                                       'marginTop': '20px'}),
            # Add the title next to the logo
            html.H3('CNPs-MFSA',
                    style={'textAlign': 'center', 'margin': '0px', 'marginTop': '20px', 'fontFamily': 'Arial',
                           'fontSize': '40px',
                           'fontWeight': 'bold'}),
        ], style={'display': 'flex', 'alignItems': 'center'}),  # Align logo and title horizontally in a single row
        # Add the description paragraph below the title
        html.P(
            'Welcome! Complex natural products (CNPs) Structure Annotation Tool was constructed based on a modular '
            'fragmentation-based structural assembly (MFSA) strategy. This tool provides functionalities to recognize target natural products, calculate molecular compositions, determine possible adduct ions, identify characteristic neutral losses, extract feature ions, calculate cosine similarity scores between standard and query spectra, and provide possible structure annotation. Please select the relevant tab to proceed.',
            style={'textAlign': 'justify', 'margin': '20px', 'marginBottom': '20px', 'fontFamily': 'Arial', 'fontSize':
                '20px'}),
    ], style={
        'backgroundColor': '#f1f1f1',  # background color
        'padding': '5px'
    }),
    dcc.Tabs([  # Create a set of tabs for different functionalities
        # Tab 1: Target Recognition
        dcc.Tab(label='1. Target Recognize', children=[
            html.Div([  # Container for the 'Target Recognize' tab
                html.Div([  # Upload components for MGF and CSV files
                    html.Div([
                        html.H4('Import MGF file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for MGF file upload
                        dcc.Upload(  # Upload component for MGF files
                            id='recognize-upload-mgf',
                            children=html.Div('Drop or Select MGF File'),  # Instructions for file selection
                            style={  # Style settings for the upload component
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple file selection
                        ),
                        html.Div(id='recognize-mgf-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder to show uploaded file paths
                    ], style={'width': '48%'}),  # Container for MGF upload component

                    # Define CSV file upload components
                    html.Div([
                        html.H4('Import CSV file (multiple files):', style={'fontFamily': 'Arial'}),
                        dcc.Upload(
                            id='recognize-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='recognize-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                    ], style={'width': '48%'})
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align MGF and CSV upload sections side-by-side

                # Section for adding new types dynamically
                html.Div(id='target-recognize-div', children=[
                    html.Button('+ Add New Type', id='add-new-type', n_clicks=0, style={'marginTop': '20px'}),
                ]),

                # Button to run extraction process
                html.Button('Run Recognition', id='run-recognition', n_clicks=0,
                            style={
                                'marginTop': '20px',
                                'fontFamily': 'Arial',
                                'padding': '10px 20px',
                                'backgroundColor': '#007bff',
                                'color': 'white',
                                'border': 'none',
                                'cursor': 'pointer',
                                'fontSize': '16px'
                            }),

                # Download link for extraction results
                html.A("Download All Results", id="recognize-download-link", download="Recognition_results.zip",
                       href="",
                       style={'fontFamily': 'Arial', 'marginLeft': '20px', 'fontSize': '16px', 'color': '#007bff'}),

                # Output placeholder for displaying messages or errors
                html.Div(id='recognize-output', style={'marginTop': '20px', 'fontFamily': 'Arial', 'color': '#007bff'}),

                # Dropdown and result table to view specific files and their results
                html.Div([
                    html.Div([
                        html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),
                        dcc.Dropdown(
                            id='recognize-file-dropdown',
                            options=[],
                            placeholder='Select a file',
                            style={'width': '200px', 'fontFamily': 'Arial'}
                        )
                    ], style={'display': 'flex', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),
                    html.Div(id='recognize-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})
                ])
            ])
        ],
            # Styling properties for the 'Target Extract' tab
            style={'border': '2px solid #d6d6d6',
                   'backgroundColor': '#f1f1f1',
                   'padding': '10px',
                   'fontFamily': 'Arial',
                   'fontSize': '16px',
                   'fontWeight': 'bold'},
            selected_style={
                'border': '2px solid #A9A9A9',
                'backgroundColor': '#A9A9A9',
                'color': 'white',
                'padding': '10px',
                'fontFamily': 'Arial',
                'fontSize': '18px',
                'fontWeight': 'bold'
            }),

        # Tab 2: Composition Calculate
        dcc.Tab(label='2. Composition Calculate', children=[
            html.Div([
                # Two main sections: MS spectra upload and formula database upload
                html.Div([
                    # Upload section for MS Spectra (multiple files)
                    html.Div([
                        html.H4('Import query MS Spectra (multiple files):', style={'fontFamily': 'Arial'}),  # Header for MS spectra upload
                        dcc.Upload(
                            id='mf-upload-csv',  # Component ID for uploading MS spectra
                            children=html.Div('Drop or Select CSV Files'),  # Instructions for file upload
                            style={  # Styling for the upload box
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Enable selection of multiple files
                        ),
                        html.Div(id='mf-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder to show uploaded file paths
                    ], style={'width': '48%'}),  # Container width for layout consistency

                    # Upload section for the formula database (POS/NEG polarity)
                    html.Div([
                        html.H4('Import formula database (POS/NEG):', style={'fontFamily': 'Arial'}),  # Header for formula database
                        dcc.Upload(
                            id='mf-upload-db',
                            children=html.Div('Drop or Select CSV Files'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=False  # Single file upload (formula database)
                        ),
                        html.Div(id='mf-db-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for database file path
                    ], style={'width': '48%'})
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align the two upload sections side-by-side

                # Section for input parameters such as Charge and Mass Tolerance (ppm)
                html.Div([
                    html.Label('Charge:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),  # Label for Charge input
                    dcc.Input(id='mf-input-charge', type='number', value=1,
                              style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'}),
                    html.Label('Mass Tolerance (ppm):', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),  # Label for Mass Tolerance
                    dcc.Input(id='mf-input-ppm', type='number', value=5,
                              style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'})
                ], style={'display': 'flex', 'justifyContent': 'start', 'alignItems': 'center', 'fontFamily': 'Arial',
                          'marginTop': '20px'}),  # Align input elements in a row

                # Button to trigger composition calculation
                html.Button('Calculate', id='mf-calculate-button', n_clicks=0, style={
                    'marginTop': '20px', 'fontFamily': 'Arial', 'fontSize': '16px', 'padding': '10px 20px',
                    'backgroundColor': '#007bff', 'color': 'white', 'border': 'none', 'cursor': 'pointer'
                }),

                # Download link for composition results
                html.A("Download All Results", id="mf-download-link", download="composition_results.zip", href="",
                       target="_blank",
                       style={'fontFamily': 'Arial', 'marginLeft': '20px', 'fontSize': '16px', 'color': '#007bff'}),

                # Section to view individual file results using dropdown and tables
                html.Div([
                    html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),
                    dcc.Dropdown(id='mf-file-dropdown', options=[], placeholder='Select a file',
                                 style={'width': '200px', 'fontFamily': 'Arial'})
                ], style={'display': 'flex', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),

                html.Div(id='mf-calculation-output', style={'marginTop': '20px', 'fontFamily': 'Arial'}),  # Placeholder for calculation messages

                html.Div(id='mf-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})  # Table for displaying calculation results
            ])
        ], style={
            'border': '2px solid #d6d6d6',
            'backgroundColor': '#f1f1f1',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '16px',
            'fontWeight': 'bold'},
            selected_style={  # Style settings for the selected tab
            'border': '2px solid #A9A9A9',
            'backgroundColor': '#A9A9A9',
            'color': 'white',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '18px',
            'fontWeight': 'bold'
        }),

        # Tab 3: Adduct Ion Calculate
        dcc.Tab(
            label='3. Adduct Ion Calculate',
            children=[
                html.Div([
                    html.Div([
                        html.Div([
                            html.H4('Import query MS Spectra (multiple files):', style={'fontFamily': 'Arial'}),
                            dcc.Upload(
                                id='adduct-upload-csv',
                                children=html.Div('Drop or Select CSV Files'),
                                style={
                                    'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                    'borderWidth': '1px', 'borderStyle': 'dashed',
                                    'borderRadius': '5px', 'textAlign': 'center',
                                    'margin': '10px', 'fontFamily': 'Arial'
                                },
                                multiple=True
                            ),
                            html.Div(id='adduct-csv-file-path',
                                     style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                        ], style={'width': '48%'}),

                        html.Div([
                            html.H4('Alternative: Import known compound MS database (POS/NEG):', style={'fontFamily':
                                                                                                     'Arial'}),
                            dcc.Upload(
                                id='adduct-upload-db',
                                children=html.Div('Drop or Select DB Files (.csv)'),
                                style={
                                    'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                    'borderWidth': '1px', 'borderStyle': 'dashed',
                                    'borderRadius': '5px', 'textAlign': 'center',
                                    'margin': '10px', 'fontFamily': 'Arial'
                                },
                                multiple=False
                            ),
                            html.Div(id='adduct-db-file-path',
                                     style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                        ], style={'width': '48%'})
                    ], style={'display': 'flex', 'justifyContent': 'space-between'}),

                    html.Div([
                        html.Label('Charge:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                        dcc.Input(
                            id='adduct-input-charge', type='number', value=1,
                            style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'}
                        ),
                        html.Label('Mass Tolerance (ppm):', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                        dcc.Input(
                            id='adduct-input-ppm', type='number', value=5,
                            style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'}
                        )
                    ], style={
                        'display': 'flex', 'justifyContent': 'start', 'alignItems': 'center',
                        'fontFamily': 'Arial', 'marginTop': '20px'
                    }),

                    html.Button(
                        'Calculate', id='adduct-calculate-button', n_clicks=0,
                        style={
                            'marginTop': '20px', 'fontFamily': 'Arial', 'fontSize': '16px', 'padding': '10px 20px',
                            'backgroundColor': '#007bff', 'color': 'white', 'border': 'none', 'cursor': 'pointer'
                        }
                    ),
                    html.A(
                        "Download All Results", id="adduct-download-link",
                        download="adduct_results.zip", href="", target="_blank",
                        style={'fontFamily': 'Arial', 'marginLeft': '20px',
                               'fontSize': '16px', 'color': '#007bff'}
                    ),

                    html.Div([
                        html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),
                        dcc.Dropdown(
                            id='adduct-file-dropdown', options=[], placeholder='Select a file',
                            style={'width': '200px', 'fontFamily': 'Arial'}
                        )
                    ], style={
                        'display': 'flex', 'alignItems': 'center',
                        'fontFamily': 'Arial', 'marginTop': '20px'
                    }),

                    html.Div(id='adduct-calculation-output', style={'marginTop': '20px', 'fontFamily': 'Arial'}),
                    html.Div(id='adduct-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})
                ])
            ],
            style={
                'border': '2px solid #d6d6d6',
                'backgroundColor': '#f1f1f1',
                'padding': '10px',
                'fontFamily': 'Arial',
                'fontSize': '16px',
                'fontWeight': 'bold'
            },
            selected_style={
                'border': '2px solid #A9A9A9',
                'backgroundColor': '#A9A9A9',
                'color': 'white',
                'padding': '10px',
                'fontFamily': 'Arial',
                'fontSize': '18px',
                'fontWeight': 'bold'
            }),

        # Tab 4: Neutral Loss Extract
        dcc.Tab(label='4. Neutral Loss Extract', children=[
            html.Div([
                # Container for MGF, CSV, and database file uploads
                html.Div([
                    # Upload section for MGF files (multiple files)
                    html.Div([
                        html.H4('Import MGF file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for MGF file upload
                        dcc.Upload(
                            id='neutralloss-upload-mgf',
                            children=html.Div('Drop or Select MGF File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple MGF files to be uploaded
                        ),
                        html.Div(id='neutralloss-mgf-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder to display uploaded file paths
                    ], style={'width': '30%'}),  # Container width for layout

                    # Upload section for CSV files (multiple files)
                    html.Div([
                        html.H4('Import CSV file (multiple files):', style={'fontFamily': 'Arial'}),
                        dcc.Upload(
                            id='neutralloss-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files to be uploaded
                        ),
                        html.Div(id='neutralloss-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                    ], style={'width': '30%'}),  # Container width for layout

                    # Upload section for the formula database (POS/NEG)
                    html.Div([
                        html.H4('Import formula database (POS/NEG):', style={'fontFamily': 'Arial'}),
                        dcc.Upload(
                            id='neutralloss-upload-db',
                            children=html.Div('Drop or Select CSV Files'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=False  # Single file for database upload
                        ),
                        html.Div(id='neutralloss-db-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                    ], style={'width': '30%'})  # Container width for layout
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align upload sections side-by-side

                # Charge and Mass Tolerance (ppm) inputs for Neutral Loss calculations
                html.Div([
                    html.Label('Charge:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                    dcc.Input(id='neutralloss-input-charge', type='number', value=1,
                              style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'}),
                    html.Label('Mass Tolerance (ppm):', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                    dcc.Input(id='neutralloss-input-ppm', type='number', value=5,
                              style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'})
                ], style={'display': 'flex', 'justifyContent': 'start', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),

                # Section to add new target loss dynamically
                html.Div(id='neutralloss-div', children=[
                    html.Button('+ Add New Target Loss', id='add-new-targetloss', n_clicks=0, style={'marginTop': '20px'}),
                ]),

                # Button to trigger neutral loss extraction process
                html.Button('Run Neutral Loss Extraction', id='run-neutralloss-extraction', n_clicks=0,
                            style={'marginTop': '20px', 'fontFamily': 'Arial', 'padding': '10px 20px',
                                   'backgroundColor': '#007bff', 'color': 'white', 'border': 'none',
                                   'cursor': 'pointer', 'fontSize': '16px'}),

                # Download link for neutral loss extraction results
                html.A("Download All Results", id="neutralloss-download-link", download="neutralloss_extraction_results.zip", href="",
                       style={'fontFamily': 'Arial', 'marginLeft': '20px', 'fontSize': '16px', 'color': '#007bff'}),

                # Output placeholder to display messages or extraction results
                html.Div(id='neutralloss-extract-output', style={'marginTop': '20px', 'fontFamily': 'Arial', 'color': '#007bff'}),

                # Dropdown and table to view specific file results
                html.Div([
                    html.Div([
                        html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),
                        dcc.Dropdown(
                            id='neutralloss-extract-file-dropdown',
                            options=[],
                            placeholder='Select a file',
                            style={'width': '200px', 'fontFamily': 'Arial'}
                        )
                    ], style={'display': 'flex', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),
                    html.Div(id='neutralloss-extract-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})
                ])
            ])
        ], style={
            'border': '2px solid #d6d6d6',
            'backgroundColor': '#f1f1f1',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '16px',
            'fontWeight': 'bold'},
            selected_style={
            'border': '2px solid #A9A9A9',
            'backgroundColor': '#A9A9A9',
            'color': 'white',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '18px',
            'fontWeight': 'bold'
        }),

        # Tab 5: Feature Ion Extract
        dcc.Tab(label='5. Feature Ion Extract', children=[
            html.Div([  # Main container for 'Feature Ion Extract' tab
                html.Div([  # Section for MGF and CSV file uploads
                    html.Div([
                        html.H4('Import MGF file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for MGF upload
                        dcc.Upload(
                            id='feature-ion-upload-mgf',
                            children=html.Div('Drop or Select MGF File'),
                            style={  # Styling for upload box
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple MGF files
                        ),
                        html.Div(id='feature-ion-mgf-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for MGF file paths
                    ], style={'width': '48%'}),  # Container for MGF upload

                    # Section for CSV file upload
                    html.Div([
                        html.H4('Import CSV file (multiple files):', style={'fontFamily': 'Arial'}),
                        dcc.Upload(
                            id='feature-ion-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='feature-ion-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for CSV file paths
                    ], style={'width': '48%'})  # Container for CSV upload
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align MGF and CSV upload sections side-by-side

                # Section for adding new ion sets dynamically
                html.Div(id='feature-ion-div', children=[
                    html.Button('+ Add New Ion Set', id='add-new-ionset', n_clicks=0, style={'marginTop': '20px'}),
                ]),

                # Button to run feature ion extraction
                html.Button('Run Feature Ion Extraction', id='run-feature-ion-extraction', n_clicks=0,
                            style={'marginTop': '20px', 'fontFamily': 'Arial', 'padding': '10px 20px',
                                   'backgroundColor': '#007bff', 'color': 'white', 'border': 'none',
                                   'cursor': 'pointer', 'fontSize': '16px'}),

                # Download link for feature ion extraction results
                html.A("Download All Results", id="feature-ion-download-link", download="featureion_extraction_results.zip", href="",
                       style={'fontFamily': 'Arial', 'marginLeft': '20px', 'fontSize': '16px', 'color': '#007bff'}),

                # Output section for messages or results
                html.Div(id='feature-ion-extract-output', style={'marginTop': '20px', 'fontFamily': 'Arial', 'color': '#007bff'}),

                # Dropdown and result table for viewing specific file results
                html.Div([
                    html.Div([
                        html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),
                        dcc.Dropdown(
                            id='feature-ion-extract-file-dropdown',
                            options=[],
                            placeholder='Select a file',
                            style={'width': '200px', 'fontFamily': 'Arial'}
                        )
                    ], style={'display': 'flex', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),
                    html.Div(id='feature-ion-extract-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})
                ])
            ])
        ], style={
            'border': '2px solid #d6d6d6',
            'backgroundColor': '#f1f1f1',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '16px',
            'fontWeight': 'bold'},
            selected_style={
            'border': '2px solid #A9A9A9',
            'backgroundColor': '#A9A9A9',
            'color': 'white',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '18px',
            'fontWeight': 'bold'
        }),

        # Tab 6: Cosine Score Calculate
        dcc.Tab(label='6. Cosine Score Calculate', children=[
            html.Div([
                # Container for query and standard MGF and CSV uploads
                html.Div([
                    # Query MGF file upload (multiple files)
                    html.Div([
                        html.H4('Import Query MGF file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for Query MGF upload
                        dcc.Upload(
                            id='cosine-query-upload-mgf',
                            children=html.Div('Drop or Select MGF File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple MGF files
                        ),
                        html.Div(id='cosine-query-mgf-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for uploaded file paths
                    ], style={'width': '25%'}),  # Set container width for layout

                    # Query CSV file upload (multiple files)
                    html.Div([
                        html.H4('Import Query CSV file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for Query CSV upload
                        dcc.Upload(
                            id='cosine-query-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='cosine-query-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                    ], style={'width': '25%'}),  # Container width for layout

                    # Standard MGF file upload (multiple files)
                    html.Div([
                        html.H4('Import Standard MGF file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for Standard MGF upload
                        dcc.Upload(
                            id='cosine-standard-upload-mgf',
                            children=html.Div('Drop or Select MGF File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple MGF files
                        ),
                        html.Div(id='cosine-standard-mgf-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                    ], style={'width': '25%'}),  # Container width for layout

                    # Standard CSV file upload (multiple files)
                    html.Div([
                        html.H4('Import Standard CSV file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for Standard CSV upload
                        dcc.Upload(
                            id='cosine-standard-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='cosine-standard-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                    ], style={'width': '25%'})  # Container width for layout
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align upload sections side-by-side

                # Input fields for specifying the mass range and intensity normalization
                html.Div([
                    html.Label('Mass Range of Product Ions (m/z):', style={'fontFamily': 'Arial', 'fontWeight': 'bold', 'marginRight': '10px'}),  # Label for mass range
                    dcc.Input(id={'type': 'cosine-mass-range-min'}, type='number', placeholder='min',
                              style={'width': '18%', 'marginRight': '10px', 'fontFamily': 'Arial'}),  # Input for minimum mass range
                    html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),  # Separator between min and max range
                    dcc.Input(id={'type': 'cosine-dropdown-and-input'}, type='text',
                              value='Use Precursor Ion', style={'width': '18%', 'fontFamily': 'Arial'})  # Placeholder input for advanced configuration
                ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '20px'}),  # Align inputs in a single row

                # Inputs for specifying intensity normalization and cosine score threshold
                html.Div([
                    html.Label('Intensity Normalization:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),  # Label for Intensity Normalization
                    dcc.Input(id='cosine-input-intensity-normalization', type='number', value=0.1,
                              style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'}),
                    html.Label('Cosine:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),  # Label for Cosine score input
                    dcc.Input(id='cosine-input-cosine-score', type='number', value=0.7,
                              style={'width': '10%', 'marginRight': '20px', 'marginLeft': '10px', 'fontFamily': 'Arial'})
                ], style={'display': 'flex', 'justifyContent': 'start', 'alignItems': 'center', 'fontFamily': 'Arial',
                          'marginTop': '20px'}),  # Align inputs side-by-side

                # Button to trigger cosine score calculation
                html.Button('Run Cosine Score Calculation', id='run-cosine-score-calculation', n_clicks=0,
                            style={'marginTop': '20px', 'fontFamily': 'Arial', 'padding': '10px 20px',
                                   'backgroundColor': '#007bff', 'color': 'white', 'border': 'none',
                                   'cursor': 'pointer', 'fontSize': '16px'}),

                # Download link for cosine score results
                html.A("Download All Results", id="cosine-download-link", download="cosine_results.zip", href="",
                       style={'fontFamily': 'Arial', 'marginLeft': '20px', 'fontSize': '16px', 'color': '#007bff'}),

                # Output placeholder for calculation messages or errors
                html.Div(id='cosine-output', style={'marginTop': '20px', 'fontFamily': 'Arial', 'color': '#007bff'}),

                # Dropdown and table to view results of individual files
                html.Div([
                    html.Div([
                        html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),  # Label for file selection
                        dcc.Dropdown(
                            id='cosine-file-dropdown',
                            options=[],  # Placeholder for dynamically populated file options
                            placeholder='Select a file',
                            style={'width': '200px', 'fontFamily': 'Arial'}
                        )
                    ], style={'display': 'flex', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),  # Container for dropdown
                    html.Div(id='cosine-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})  # Table for displaying file-specific results
                ])
            ])
        ], style={
            'border': '2px solid #d6d6d6',
            'backgroundColor': '#f1f1f1',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '16px',
            'fontWeight': 'bold'},
            selected_style={
            'border': '2px solid #A9A9A9',
            'backgroundColor': '#A9A9A9',
            'color': 'white',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '18px',
            'fontWeight': 'bold'
        }),

        # Tab 7: Annotation
        dcc.Tab(label='7. Annotation', children=[
            html.Div([
                # Container for different file upload sections related to Annotation
                html.Div([
                    # Upload section for the Query CSV file (multiple files)
                    html.Div([
                        html.H4('Import Query CSV file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for Query CSV file upload
                        dcc.Upload(
                            id='annotation-query-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='annotation-query-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for file paths
                    ], style={'width': '33%'}),  # Set container width for layout

                    # Upload section for the pseudo-DB file (single file)
                    html.Div([
                        html.H4('Import pseudo-DB file:', style={'fontFamily': 'Arial'}),  # Header for Skeleton CSV
                        # file upload
                        dcc.Upload(
                            id='annotation-pseudo-DB-upload-db',
                            children=html.Div('Drop or Select DB File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=False
                        ),
                        html.Div(id='annotation-pseudo-DB-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for file path
                    ], style={'width': '33%'}),

                    # Upload section for the Substituent CSV file (multiple files)
                    html.Div([
                        html.H4('Import Substituent CSV file:', style={'fontFamily': 'Arial'}),  # Header for Substituent CSV file upload
                        dcc.Upload(
                            id='annotation-substituent-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='annotation-substituent-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for file paths
                    ], style={'width': '33%'})  # Set container width for layout
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align the three upload sections side-by-side

                # Button to run the annotation process
                html.Button('Run Annotation', id='run-annotation', n_clicks=0,
                            style={'marginTop': '20px', 'fontFamily': 'Arial', 'padding': '10px 20px',
                                   'backgroundColor': '#007bff', 'color': 'white', 'border': 'none',
                                   'cursor': 'pointer', 'fontSize': '16px'}),

                # Download link for annotation results
                html.A("Download All Results", id="annotation-download-link", download="annotation_results.zip", href="",
                       style={'fontFamily': 'Arial', 'marginLeft': '20px', 'fontSize': '16px', 'color': '#007bff'}),

                # Output placeholder for messages or results
                html.Div(id='annotation-output', style={'marginTop': '20px', 'fontFamily': 'Arial', 'color': '#007bff'}),

                # Dropdown and results table to view specific file results
                html.Div([
                    html.Div([
                        html.H4('Select file to view results:', style={'fontFamily': 'Arial', 'marginRight': '20px'}),  # Header for file selection
                        dcc.Dropdown(
                            id='annotation-file-dropdown',
                            options=[],  # Placeholder for dynamically populated file options
                            placeholder='Select a file',
                            style={'width': '200px', 'fontFamily': 'Arial'}
                        )
                    ], style={'display': 'flex', 'alignItems': 'center', 'fontFamily': 'Arial', 'marginTop': '20px'}),  # Align dropdown and label
                    html.Div(id='annotation-results-table', style={'marginTop': '20px', 'fontFamily': 'Arial'})  # Table to display results
                ])
            ])
        ], style={  # Style settings for the 'Annotation' tab
            'border': '2px solid #d6d6d6',
            'backgroundColor': '#f1f1f1',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '16px',
            'fontWeight': 'bold'},
            selected_style={
            'border': '2px solid #A9A9A9',
            'backgroundColor': '#A9A9A9',
            'color': 'white',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '18px',
            'fontWeight': 'bold'
        }),

        # Tab 8: Visualization
        dcc.Tab(label='8. Visualization', children=[
            html.Div([
                html.Div([
                    html.H4('Import Query CSV file (multiple files):', style={'fontFamily': 'Arial'}),
                    dcc.Upload(
                        id='visualization-query-upload-csv',
                        children=html.Div('Drop or Select CSV File'),
                        style={
                            'width': '90%', 'height': '50px', 'lineHeight': '50px',
                            'borderWidth': '1px', 'borderStyle': 'dashed',
                            'borderRadius': '5px', 'textAlign': 'center',
                            'margin': '10px', 'fontFamily': 'Arial'
                        },
                        multiple=True
                    ),
                    html.Div(id='visualization-query-csv-file-path',
                             style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})
                ]),

                html.Button('Run Visualization', id='run-visualization', n_clicks=0,
                            style={'marginTop': '20px', 'fontFamily': 'Arial', 'padding': '10px 20px',
                                   'backgroundColor': '#007bff', 'color': 'white', 'border': 'none',
                                   'cursor': 'pointer', 'fontSize': '16px'}),

                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.Div(
                                dcc.Graph(id='scatter-plot', style={"height": "700px"}),
                                style={'width': '65%', 'paddingRight': '20px'}
                            ),

                            html.Div([
                                html.H5("Annotation Structures", style={"fontSize": "20px", "marginBottom": "10px"}),

                                html.Div(id='structure-gallery', style={
                                    "maxHeight": "650px",
                                    "overflowY": "scroll",
                                    "overflowX": "auto",
                                    "border": "1px solid #ddd",
                                    "padding": "10px",
                                    "backgroundColor": "#fafafa",
                                    "display": "flex",
                                    "flexWrap": "wrap",
                                    "gap": "20px"
                                }),

                                html.Div(id='structure-modal-container')
                            ], style={'width': '35%'})
                        ], style={'display': 'flex', 'marginTop': '20px'})
                    ])
                ], style={
                    'border': '1px solid #ccc',
                    'borderRadius': '5px',
                    'padding': '10px',
                    'marginTop': '50px',
                    'marginBottom': '30px',
                    'marginLeft': '50px',
                    'marginRight': '50px',
                }),

                html.Div([
                    html.Div(id="summary-bar-count", style={'marginBottom': '30px'}),
                    html.Div(id="summary-bar-area")
                ]),

                html.Div(id='filtered-table-container', style={'margin': '50px 50px 100px 50px'})
            ])
        ], style={
            'border': '2px solid #d6d6d6',
            'backgroundColor': '#f1f1f1',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '16px',
            'fontWeight': 'bold'
        }, selected_style={
            'border': '2px solid #A9A9A9',
            'backgroundColor': '#A9A9A9',
            'color': 'white',
            'padding': '10px',
            'fontFamily': 'Arial',
            'fontSize': '18px',
            'fontWeight': 'bold'
        })
    ], style={'fontFamily': 'Arial'})
])


# =======================================
#  Dash callback and function of each tab
# =======================================
"""
Tab 1 Related Functions for CSV and MGF File Processing and Recognition
"""
# Callback to update the CSV file upload status on Tab 1
@app.callback(
    Output('recognize-csv-file-path', 'children'),
    Input('recognize-upload-csv', 'contents'),
    State('recognize-upload-csv', 'filename')
)
def update_recognize_csv_file_path(csv_contents, csv_filenames):
    """
    Update the status message to reflect the uploaded CSV files.

    Parameters:
        csv_contents: Contents of the uploaded CSV files.
        csv_filenames: Filenames of the uploaded CSV files.

    Returns:
        A string message displaying the names of uploaded CSV files or indicating no upload.
    """
    print(f"CSV Contents Received: {csv_contents}")
    print(f"CSV Filenames Received: {csv_filenames}")
    if csv_contents and csv_filenames:
        # Save each uploaded CSV file and generate paths
        csv_paths = [save_tab1_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        print(f"CSV Paths Saved: {csv_paths}")
        # Format the message displaying uploaded file names
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        csv_path_text = "No CSV file uploaded"
    print(f"CSV Path Text: {csv_path_text}")
    return csv_path_text

# Callback to update the MGF file upload status on Tab 1
@app.callback(
    Output('recognize-mgf-file-path', 'children'),
    Input('recognize-upload-mgf', 'contents'),
    State('recognize-upload-mgf', 'filename')
)
def update_recognize_mgf_file_path(mgf_contents, mgf_filenames):
    """
    Update the status message to reflect the uploaded MGF files.

    Parameters:
        mgf_contents: Contents of the uploaded MGF files.
        mgf_filenames: Filenames of the uploaded MGF files.

    Returns:
        A string message displaying the names of uploaded MGF files or indicating no upload.
    """
    print(f"MGF Contents Received: {mgf_contents}")
    print(f"MGF Filenames Received: {mgf_filenames}")
    if mgf_contents and mgf_filenames:
        # Save each uploaded MGF file and generate paths
        mgf_paths = [save_tab1_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
        print(f"MGF Paths Saved: {mgf_paths}")
        # Format the message displaying uploaded file names
        mgf_path_text = f"Uploaded MGF files: {', '.join(mgf_filenames)}"
    else:
        mgf_path_text = "No MGF file uploaded"
    print(f"MGF Path Text: {mgf_path_text}")
    return mgf_path_text

# Callback to validate and update dropdown and input field values in dynamic forms on Tab 1
@app.callback(
    Output({'type': 'dropdown-and-input', 'index': MATCH}, 'value'),
    [Input({'type': 'dropdown-and-input', 'index': MATCH}, 'n_submit')],
    [Input({'type': 'dropdown-and-input', 'index': MATCH}, 'value')]
)
def update_value(n_submit, value):
    """
    Validate and update input values for dropdown and input fields in dynamic forms.

    Parameters:
        n_submit: Number of times the input has been submitted.
        value: The input value, which may be text or numeric.

    Returns:
        The validated input value, which may be converted to a float or retained as a string.
    """
    print(f"Dropdown/Input Value Received: {value}, n_submit: {n_submit}")
    if value == "Use Precursor Ion":
        return value
    try:
        # Attempt to convert the input to a float if it's numeric
        float_value = float(value)
        print(f"Converted Value to Float: {float_value}")
        return float_value
    except ValueError:
        print(f"ValueError encountered. Returning value as is: {value}")
        return value

# Callback to dynamically add or remove type sections in the recognition form on Tab 1
@app.callback(
    Output('target-recognize-div', 'children'),
    [Input('add-new-type', 'n_clicks'),
     Input({'type': 'remove-button-recognize', 'index': ALL}, 'n_clicks')],
    [State('target-recognize-div', 'children')]
)
def modify_types(add_clicks, remove_clicks, children):
    """
    Add new sections for different types of recognition or remove them as needed.

    Parameters:
        add_clicks: Number of clicks on the "Add New Type" button.
        remove_clicks: Number of clicks on any remove button.
        children: Current list of children components in the dynamic form.

    Returns:
        Updated list of children components reflecting the addition or removal of input sections.
    """
    print(f"Add Clicks: {add_clicks}, Remove Clicks: {remove_clicks}")
    print(f"Current Children: {children}")

    # Add a new type section if the "Add New Type" button is clicked
    if add_clicks > 0:
        new_index = len(children) - 1  # Subtract 1 because the last child is the "Add New Type" button
        new_input_area = html.Div(id={'type': 'input-area', 'index': new_index + 1}, children=[
            # Dynamically generated Type Name input field
            html.Div([
                html.Label(f'Type Name {new_index + 1}:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(
                    id={'type': 'input-type-name', 'index': new_index + 1},
                    type='text',
                    placeholder=f'e.g., Daphnane',
                    style={'width': '35%', 'fontFamily': 'Arial', 'marginLeft': '10px'}
                )
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '20px'}),
            html.H4(f'Type {new_index + 1} Feature Formulas:', style={'fontFamily': 'Arial'}),
            dcc.Textarea(
                id={'type': 'input-type-feature-formulas', 'index': new_index + 1},
                placeholder='Input Type feature ion formulas',
                style={'width': '94%', 'height': '100px', 'fontFamily': 'Arial', 'marginLeft': '1%'}
            ),
            html.Div([
                html.Label('Mass Range of Product Ions (m/z):', style={'fontFamily': 'Arial', 'fontWeight': 'bold', 'marginRight': '10px'}),
                dcc.Input(id={'type': 'mass-range-min', 'index': new_index + 1}, type='number', placeholder='min',
                          style={'width': '18%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),  # Ensure "to" is always displayed
                dcc.Input(id={'type': 'dropdown-and-input', 'index': new_index + 1}, type='text',
                          value='Use Precursor Ion', style={'width': '18%', 'fontFamily': 'Arial'}),
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'}),
            html.Div([
                html.Label('Charge:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'charge', 'index': new_index + 1}, type='number', placeholder='e.g., +1/-1',
                          style={'width': '10%', 'marginLeft': '10px', 'marginRight': '20px', 'fontFamily': 'Arial'}),
                html.Label('Mass Tolerance (ppm):', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'tolerance', 'index': new_index + 1}, type='number', placeholder='e.g., 5',
                          style={'width': '10%', 'marginLeft': '10px', 'marginRight': '20px', 'fontFamily': 'Arial'}),
                html.Label('Hit Score:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'hit-score', 'index': new_index + 1}, type='number', placeholder='e.g., 5',
                          style={'width': '10%', 'marginLeft': '10px', 'fontFamily': 'Arial'}),
                html.Button('✖', id={'type': 'remove-button', 'index': new_index + 1}, style={'marginLeft': '10px', 'color': 'red'})
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'})
        ])
        children.insert(-1, new_input_area)  # Insert before the "Add New Type" button

    # Context for determining if a remove button was clicked
    ctx = dash.callback_context
    print(f"Callback Context: {ctx}")
    if ctx.triggered:
        triggered = ctx.triggered[0]['prop_id'].split('.')[0]
        print(f"Triggered Event: {triggered}")
        # Remove the input section corresponding to the clicked remove button
        if 'remove-button-recognize' in triggered:
            button_index = eval(triggered)['index']
            children = [
                child for child in children
                if child['props']['id'].get('index') != button_index
            ]
            print(f"Removed Section with Index: {button_index}")

    print(f"Updated Children: {children}")
    return children

# Callback to handle running the target extraction process
@app.callback(
    Output('recognize-output', 'children'),
    Input('run-recognition', 'n_clicks'),
    State('recognize-upload-csv', 'contents'),
    State('recognize-upload-mgf', 'contents'),
    State('recognize-upload-csv', 'filename'),
    State('recognize-upload-mgf', 'filename'),
    State({'type': 'input-type-name', 'index': ALL}, 'value'),
    State({'type': 'input-type-feature-formulas', 'index': ALL}, 'value'),
    State({'type': 'mass-range-min', 'index': ALL}, 'value'),
    State({'type': 'dropdown-and-input', 'index': ALL}, 'value'),
    State({'type': 'charge', 'index': ALL}, 'value'),
    State({'type': 'tolerance', 'index': ALL}, 'value'),
    State({'type': 'hit-score', 'index': ALL}, 'value')
)
def run_target_recognition(n_clicks, csv_contents, mgf_contents, csv_filenames, mgf_filenames,
                          type_names, type_feature_formulas, mass_ranges_min, mass_ranges_max_options,
                          charges, tolerances, hit_scores):
    """
    Execute the recognition process based on user-defined criteria and input files.

    Parameters:
        n_clicks: Number of clicks on the "Run Recognition" button.
        csv_contents: Contents of the uploaded CSV files.
        mgf_contents: Contents of the uploaded MGF files.
        csv_filenames: Filenames of the uploaded CSV files.
        mgf_filenames: Filenames of the uploaded MGF files.
        type_names: Names of types being processed.
        type_feature_formulas: Formulas associated with each type.
        mass_ranges_min: Minimum mass ranges for filtering.
        mass_ranges_max_options: Maximum mass range options (may be "Use Precursor Ion" or numeric).
        charges: Charges associated with each type for mass calculation.
        tolerances: Tolerance levels for filtering.
        hit_scores: Minimum number of hits required for recognition.

    Returns:
        A message indicating the completion of the recognition process and the location of saved files.
    """
    print("recognition started")
    print(f"Parameters - n_clicks: {n_clicks}, csv_contents: {csv_contents}, mgf_contents: {mgf_contents}")
    print(f"csv_filenames: {csv_filenames}, mgf_filenames: {mgf_filenames}")
    print(f"type_names: {type_names}, type_feature_formulas: {type_feature_formulas}")
    print(f"mass_ranges_min: {mass_ranges_min}, mass_ranges_max_options: {mass_ranges_max_options}")
    print(f"charges: {charges}, tolerances: {tolerances}, hit_scores: {hit_scores}")

    # Check if conditions for extraction are met
    if n_clicks > 0 and csv_contents and mgf_contents:
        all_output_paths = []

        # Process each MGF file uploaded
        for mgf_content, mgf_filename in zip(mgf_contents, mgf_filenames):
            print(f"Processing MGF file: {mgf_filename}")

            mgf_file_path = save_tab1_file(mgf_filename, mgf_content)
            mgf_prefix = mgf_filename.replace('.mgf', '')

            # Match corresponding CSV files based on MGF file name
            matched_csv_files = [csv_file for csv_file in csv_filenames if
                                 csv_file.startswith(mgf_prefix) and csv_file.endswith('_quant.csv')]
            print(f"Matched CSV files: {matched_csv_files}")

            if not matched_csv_files:
                print(f"No matching CSV files found for MGF file: {mgf_filename}")
                continue

            csv_filename = matched_csv_files[0]
            csv_content = csv_contents[csv_filenames.index(csv_filename)]
            csv_file_path = save_tab1_file(csv_filename, csv_content)

            print(f"Processing CSV file: {csv_filename}")

            output_mgf_file = os.path.join(UPLOAD_DIR, f'{mgf_filename}')
            output_csv_file = os.path.join(UPLOAD_DIR, f'{csv_filename}')

            all_spectra = []
            all_ids = []

            # Initialize a dictionary to store scores for each type
            scores_dict = {type_name: {} for type_name in type_names}

            # Process each type name and perform feature extraction
            for i, type_name in enumerate(type_names):
                print(f"Processing type: {type_name}")
                feature_formulas = [x.strip() for x in type_feature_formulas[i].split(',')]
                min_mass = mass_ranges_min[i]
                charge = charges[i]
                tolerance_ppm = tolerances[i]
                hit_score = hit_scores[i]

                # Determine maximum mass based on the option selected
                max_mass_option = mass_ranges_max_options[i]
                if max_mass_option == "Use Precursor Ion":
                    max_mass_value = None
                else:
                    try:
                        max_mass_value = float(max_mass_option)
                    except ValueError:
                        max_mass_value = None

                print(f"Parameters for Type {type_name} - min_mass: {min_mass}, charge: {charge}, tolerance_ppm: {tolerance_ppm}, hit_score: {hit_score}, max_mass_value: {max_mass_value}")

                # Extract spectra from the MGF file matching the criteria
                spectra, type_hit_scores = detect_recognize_spectra(mgf_file_path, feature_formulas, min_mass, charge,
                                                                  tolerance_ppm, hit_score, max_mass_option, max_mass_value)

                print(f"Spectra found for Type {type_name}: {len(spectra)}")
                print(f"Hit scores for Type {type_name}: {type_hit_scores}")

                # Store extracted spectra and their scores
                ids = [int(spectrum['params']['title']) for spectrum in spectra]
                scores_dict[type_name] = {id_: score for id_, score in zip(ids, type_hit_scores)}

                all_spectra.extend(spectra)
                all_ids.extend(ids)

            print(f"Total recognized spectra: {len(all_spectra)}, Unique IDs: {len(set(all_ids))}")

            # Filter and sort unique spectra based on title
            unique_spectra = []
            seen_titles = set()
            for spectrum in all_spectra:
                title = spectrum['params']['title']
                if title not in seen_titles:
                    seen_titles.add(title)
                    unique_spectra.append(spectrum)

            unique_spectra.sort(key=lambda x: int(x['params']['title']))
            print(f"Unique spectra count after filtering: {len(unique_spectra)}")

            # Save unique spectra to an output MGF file
            mgf.write(unique_spectra, output=output_mgf_file)

            # Attempt to load the corresponding CSV file, handling errors gracefully
            try:
                input_csv = pd.read_csv(csv_file_path, index_col='row ID')
                print(f"CSV loaded successfully: {csv_filename}")
            except pd.errors.ParserError:
                input_csv = pd.read_csv(csv_file_path, index_col='row ID', error_bad_lines=False)
                print(f"CSV loaded with errors: {csv_filename}")

            # Filter rows in the CSV file based on extracted spectra IDs
            input_IDs = input_csv.index.tolist()
            drop_IDs = [x for x in input_IDs if x not in all_ids]
            focal_csv = input_csv.drop(drop_IDs, axis=0)
            print(f"Remaining rows after filtering: {focal_csv.shape[0]}")

            if focal_csv.empty:
                print(f"No matching spectra found in the CSV file {csv_filename}.")
                continue

            # Add score columns for each type based on the calculated scores from the detection process
            for i, type_name in enumerate(type_names):
                focal_csv[f'{type_name} Score'] = focal_csv.index.map(lambda x: scores_dict[type_name].get(x, 0))

            # Create a list of score column names for easy access
            score_columns = [f'{type_name} Score' for type_name in type_names]

            # Filter out rows where all the score columns are zero
            focal_csv = focal_csv[(focal_csv[score_columns] != 0).any(axis=1)]
            print(f"Rows after filtering for non-zero scores: {focal_csv.shape[0]}")

            # Determine the classification based on the logic for a two-class scenario:
            # - If Type 1 has a non-zero score and Type 2 has a zero score, classify as Type 1.
            # - If both Type 1 and Type 2 have non-zero scores (or only Type 2 has a score), classify as Type 2.
            if len(score_columns) == 2:
                focal_csv['type'] = focal_csv.apply(
                    lambda row: type_names[0] if (row[score_columns[0]] != 0 and row[score_columns[1]] == 0)
                    else type_names[1],
                    axis=1
                )
            else:
                # For scenarios with more than two types, fallback to the original maximum score logic
                focal_csv['type'] = focal_csv[score_columns].idxmax(axis=1).str.replace(' Score', '')

            # Save the filtered CSV data to an output file
            focal_csv.to_csv(output_csv_file, encoding='UTF-8')
            all_output_paths.append((output_mgf_file, output_csv_file))
            print(f"Output MGF file saved to: {output_mgf_file}")
            print(f"Output CSV file saved to: {output_csv_file}")


        print("Recognition complete")
        return (f"Recognition and classification complete. Processed files saved to: "
                f"{', '.join([f'MGF: {path[0]}, CSV: {path[1]}' for path in all_output_paths])}")

    print("No recognition performed")
    return "No recognition performed. Please upload both CSV and MGF files."


# Callback to update dropdown options after extraction is run
@app.callback(
    Output('recognize-file-dropdown', 'options'),
    Input('run-recognition', 'n_clicks'),
    State('recognize-upload-csv', 'filename'),
    State('recognize-upload-mgf', 'filename')
)
def update_dropdown_options(n_clicks, csv_filenames, mgf_filenames):
    """
    Populate dropdown options with the names of processed files.

    Parameters:
        n_clicks: Number of clicks on the "Run Recognition" button.
        csv_filenames: List of uploaded CSV file names.
        mgf_filenames: List of uploaded MGF file names.

    Returns:
        A list of options for the dropdown containing the names of the processed CSV and MGF files.
    """
    print(f"Update dropdown called with n_clicks: {n_clicks}")
    print(f"CSV Filenames: {csv_filenames}")
    print(f"MGF Filenames: {mgf_filenames}")

    if n_clicks > 0:
        # Combine CSV and MGF filenames into dropdown options
        options = [{'label': file, 'value': file} for file in csv_filenames + mgf_filenames]
        print(f"Dropdown options generated: {options}")
        return options
    print("No options to update in dropdown.")
    return []

# Callback to display the results of extraction in a data table on Tab 1
@app.callback(
    Output('recognize-results-table', 'children'),
    Input('recognize-file-dropdown', 'value'),
    State('recognize-upload-csv', 'filename'),
    State('recognize-upload-mgf', 'filename')
)
def update_results_table(selected_file, csv_filenames, mgf_filenames):
    """
    Display extracted results in a data table for the selected file.

    Parameters:
        selected_file: The filename selected from the dropdown.
        csv_filenames: List of uploaded CSV file names.
        mgf_filenames: List of uploaded MGF file names.

    Returns:
        A Dash DataTable component displaying the extracted data from the selected file.
    """
    print(f"Selected file for results table: {selected_file}")
    print(f"CSV Filenames: {csv_filenames}")
    print(f"MGF Filenames: {mgf_filenames}")

    if selected_file:
        # Check if the selected file is a CSV file and exists in the output directory
        if selected_file in csv_filenames:
            csv_path = os.path.join(UPLOAD_DIR, f'output_{selected_file}')
            print(f"Looking for CSV file at path: {csv_path}")
            if os.path.exists(csv_path):
                print(f"CSV file found. Reading file: {csv_path}")
                # Read the CSV data into a DataFrame
                df = pd.read_csv(csv_path)
                print(f"Data loaded into DataFrame with shape: {df.shape}")
                # Return the data in a formatted Dash DataTable
                return dash_table.DataTable(
                    data=df.to_dict('records'),
                    columns=[{'name': i, 'id': i} for i in df.columns],
                    style_table={'overflowX': 'auto'},
                    style_cell={
                        'fontFamily': 'Arial',
                        'textAlign': 'left',
                        'padding': '10px',
                        'whiteSpace': 'normal',
                        'height': 'auto',
                    },
                    style_header={
                        'fontFamily': 'Arial',
                        'fontWeight': 'bold',
                        'backgroundColor': '#f9f9f9',
                        'border': '1px solid black'
                    },
                    style_data={
                        'fontFamily': 'Arial',
                        'border': '1px solid black',
                        'backgroundColor': 'white',
                    }
                )
            else:
                print(f"CSV file does not exist: {csv_path}")
    print("No file selected or no results available.")
    return "No file selected."

# Callback to generate a download link for the extracted results in a ZIP file
@app.callback(
    Output('recognize-download-link', 'href'),
    Input('run-recognition', 'n_clicks'),
    State('recognize-upload-csv', 'filename'),
    State('recognize-upload-mgf', 'filename')
)
def create_zip_download_link(n_clicks, csv_filenames, mgf_filenames):
    """
    Create a ZIP file containing all processed CSV and MGF results for download.

    Parameters:
        n_clicks: Number of clicks on the "Run Extraction" button.
        csv_filenames: List of uploaded CSV file names.
        mgf_filenames: List of uploaded MGF file names.

    Returns:
        A base64-encoded download link for the ZIP file containing the processed data.
    """
    print(f"Create ZIP download link called with n_clicks: {n_clicks}")
    print(f"CSV Filenames for ZIP: {csv_filenames}")
    print(f"MGF Filenames for ZIP: {mgf_filenames}")

    if n_clicks > 0:
        recognized_folder = 'temp'
        if not os.path.exists(recognized_folder):
            os.makedirs(recognized_folder)

        zip_filename = os.path.join(recognized_folder, "recognition_results.zip")
        print(f"Creating ZIP file at: {zip_filename}")

        # Create a ZIP file for all processed results
        with zipfile.ZipFile(zip_filename, 'w') as zipf:
            # processing CSV
            for csv_file in csv_filenames:
                csv_base = os.path.splitext(csv_file)[0]
                output_csv_path = os.path.join(recognized_folder, f'{csv_base}-re.csv')
                print(f"Checking CSV file for ZIP: {output_csv_path}")
                if os.path.exists(output_csv_path):
                    zipf.write(output_csv_path, os.path.basename(output_csv_path))
                    print(f"Added CSV to ZIP: {output_csv_path}")
                else:
                    print(f"CSV file not found for ZIP: {output_csv_path}")

            # Add processed MGF files to the ZIP
            for mgf_file in mgf_filenames:
                mgf_base = os.path.splitext(mgf_file)[0]
                output_mgf_path = os.path.join(recognized_folder, f'{mgf_base}-re.mgf')
                print(f"Checking MGF file for ZIP: {output_mgf_path}")
                if os.path.exists(output_mgf_path):
                    zipf.write(output_mgf_path, os.path.basename(output_mgf_path))
                    print(f"Added MGF to ZIP: {output_mgf_path}")
                else:
                    print(f"MGF file not found for ZIP: {output_mgf_path}")

        print(f"ZIP file created at: {zip_filename}")
        # return '/download/recognition_results.zip'
        return '/download/recognition_results.zip'

    print("No ZIP link created; no files to download.")
    return ""

# Calculate m/z range based on formula, charge, and tolerance
def recognize_mz_range(formula: str, charge: int, tolerance_ppm: float) -> tuple:
    """
    Calculate the mass-to-charge (m/z) range for a given formula and tolerance.

    Parameters:
        formula: Chemical formula for mass calculation.
        charge: Charge state for m/z calculation.
        tolerance_ppm: Tolerance level in parts per million (PPM).

    Returns:
        A tuple containing the minimum and maximum m/z values.
    """
    print(f"Calculating m/z range for formula: {formula}, charge: {charge}, tolerance: {tolerance_ppm} ppm")

    # Calculate the m/z value and adjust based on the tolerance
    try:
        mz = mass.calculate_mass(formula=formula, charge=charge)
        print(f"Calculated m/z: {mz}")
    except Exception as e:
        print(f"Error calculating m/z for formula: {formula}, error: {e}")
        raise e

    mz_min = mz - mz * tolerance_ppm / 1e6
    mz_max = mz + mz * tolerance_ppm / 1e6
    print(f"Calculated m/z range: min={mz_min}, max={mz_max}")

    return mz_min, mz_max


# Detect and extract spectra from an MGF file based on provided criteria
def detect_recognize_spectra(filename: str, feature_formulas: list, min_mass: float, charge: int, tolerance_ppm: float,
                           hit_score: int, max_mass_option: str, max_mass_value: float = None):
    """
    Recognize relevant spectra from an MGF file based on the provided feature formulas, mass range, and scoring
    criteria.

    Parameters:
        filename: Path to the MGF file.
        feature_formulas: List of chemical formulas to match.
        min_mass: Minimum mass for filtering spectra.
        charge: Charge state for m/z calculation.
        tolerance_ppm: Tolerance level in PPM for mass range.
        hit_score: Minimum number of hits required to consider a spectrum.
        max_mass_option: Option for setting the maximum mass range ("Use Precursor Ion" or a fixed value).
        max_mass_value: Numeric value for the maximum mass range (if not using "Use Precursor Ion").

    Returns:
        A list of matching spectra and their corresponding hit scores.
    """
    print(f"Detecting spectra from file: {filename}")
    print(f"Feature formulas: {feature_formulas}")
    print(f"Min mass: {min_mass}, Charge: {charge}, Tolerance: {tolerance_ppm} ppm")
    print(f"Hit score threshold: {hit_score}, Max mass option: {max_mass_option}, Max mass value: {max_mass_value}")

    spectra_list = []
    hit_scores_list = []

    # Read spectra from the MGF file
    with mgf.read(filename) as spectra:
        for spectrum in spectra:
            # Determine the maximum mass for filtering
            if max_mass_option == "Use Precursor Ion":
                max_mass = spectrum['params']['pepmass'][0]
                if max_mass is None:
                    print(f"Skipping spectrum with missing precursor ion mass.")
                    continue
            else:
                max_mass = max_mass_value

            if max_mass is None:
                max_mass = 2000  # Default to a large value if max mass is not specified

            print(f"Processing spectrum with title: {spectrum['params'].get('title', 'Unknown')}, Max mass: {max_mass}")

            # Filter product ions based on mass range
            product_ions = spectrum['m/z array']
            intensity = spectrum['intensity array']
            ind = np.nonzero((product_ions >= min_mass) & (product_ions <= max_mass))[0]
            print(f"Number of product ions within mass range: {len(ind)}")

            if len(ind) == 0:
                print("No product ions found within specified mass range. Skipping spectrum.")
                continue

            # Normalize intensities for scoring
            intensity_max = np.max(intensity[ind])
            intensity_min = np.min(intensity[ind])
            normalized_intensity = (intensity[ind] - intensity_min) / (intensity_max - intensity_min)
            nor_intensity = np.zeros_like(intensity)
            nor_intensity[ind] = normalized_intensity
            print(f"Normalized intensities for spectrum.")

            # Count hits based on the provided formulas and tolerance
            hits = 0
            for intensity_value, product_ion in zip(nor_intensity, product_ions):
                if intensity_value > 0.05:
                    for formula in feature_formulas:
                        mz_min, mz_max = recognize_mz_range(formula, charge, tolerance_ppm)
                        if mz_min <= product_ion <= mz_max:
                            hits += 1
                            break  # Stop checking further formulas once a hit is found for this product ion
            print(f"Hits found for spectrum: {hits}")

            # Append spectra with sufficient hits to the results list
            if hits >= hit_score:
                spectra_list.append(spectrum)
                hit_scores_list.append(hits)
                print(f"Spectrum added to results with {hits} hits.")

    print(f"Total spectra matching criteria: {len(spectra_list)}")
    return spectra_list, hit_scores_list


# Save an uploaded file to a specified folder for Tab 1
def save_tab1_file(name, content, tab_folder='tab_1_uploads'):
    """
    Save an uploaded file to a specific folder for Tab 1.

    Parameters:
        name: Name of the file being saved.
        content: Base64-encoded content of the file.
        tab_folder: Name of the folder to save the file in (default is 'tab_1_uploads').

    Returns:
        The full path to the saved file.
    """
    print(f"Saving file - Name: {name}, Folder: {tab_folder}")

    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Created folder: {folder_path}")

    # Decode the content and save the file
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)
    print(f"Saving decoded content to: {file_path}")
    with open(file_path, 'wb') as f:
        f.write(decoded)
    print(f"File saved successfully: {file_path}")
    return file_path


"""
Tab 2 Related Functions for CSV and Database File Processing and Elemental Composition Calculation
"""

# Callback to update the status of uploaded CSV and database files on Tab 2
@app.callback(
    [Output('mf-csv-file-path', 'children'), Output('mf-db-file-path', 'children')],
    [Input('mf-upload-csv', 'contents'), Input('mf-upload-db', 'contents')],
    [State('mf-upload-csv', 'filename'), State('mf-upload-db', 'filename')]
)
def update_file_paths(csv_contents, db_content, csv_filenames, db_filename):
    """
    Update and display the paths of uploaded CSV and DB files.

    Parameters:
        csv_contents: Contents of the uploaded CSV files.
        db_content: Content of the uploaded database file.
        csv_filenames: Filenames of the uploaded CSV files.
        db_filename: Filename of the uploaded database file.

    Returns:
        A tuple of strings displaying the names of uploaded CSV and DB files or indicating no upload.
    """
    print(f"Received CSV contents: {csv_contents}")
    print(f"Received DB content: {db_content}")
    print(f"CSV filenames: {csv_filenames}")
    print(f"DB filename: {db_filename}")

    # Process CSV file uploads
    if csv_contents and csv_filenames:
        csv_paths = [save_tab2_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        print(f"Saved CSV file paths: {csv_paths}")
        csv_path_text = f"Uploaded {len(csv_filenames)} CSV file(s): {', '.join(csv_filenames)}"
    else:
        csv_path_text = "No CSV file uploaded"

    # Process database file upload
    if db_content and db_filename:
        db_path = save_tab2_file(db_filename, db_content)
        print(f"Saved DB file path: {db_path}")
        db_path_text = f"Uploaded DB file: {db_filename}"
    else:
        db_path_text = "No DB file uploaded"

    print(f"CSV path text: {csv_path_text}")
    print(f"DB path text: {db_path_text}")
    return csv_path_text, db_path_text


@app.callback(
    [Output('mf-calculation-output', 'children'),
     Output('mf-download-link', 'href'), Output('mf-file-dropdown', 'options')],
    Input('mf-calculate-button', 'n_clicks'),
    [State('mf-upload-csv', 'contents'), State('mf-upload-csv', 'filename')],
    State('mf-upload-db', 'contents'),
    State('mf-input-charge', 'value'),
    State('mf-input-ppm', 'value')
)
def calculate_composition(n_clicks, csv_contents, csv_filenames, db_content, charge, ppm):
    """
    Calculate elemental compositions for the uploaded CSV data based on the database and user input parameters.

    Parameters:
        n_clicks: Number of clicks on the "Calculate" button.
        csv_contents: Contents of the uploaded CSV files.
        csv_filenames: Filenames of the uploaded CSV files.
        db_content: Content of the uploaded database file.
        charge: Charge state for m/z calculation.
        ppm: Tolerance level in parts per million (PPM) for m/z matching.

    Returns:
        A message indicating the completion of calculations, a download link for the ZIP file with results,
        and updated dropdown options for the processed files.
    """
    print(f"Calculation button clicked {n_clicks} times.")
    print(f"Charge state: {charge}, PPM tolerance: {ppm}")
    print(f"CSV contents: {csv_contents}")
    print(f"CSV filenames: {csv_filenames}")
    print(f"DB content: {db_content}")

    if n_clicks > 0 and csv_contents and db_content and charge and ppm:
        zip_filename = f'composition_formula_results.zip'
        zip_filepath = os.path.join('temp', zip_filename)

        print(f"ZIP file path: {zip_filepath}")

        # Ensure the temporary directory for results exists
        if not os.path.exists('temp'):
            os.makedirs('temp')
            print("Created temporary directory for results.")

        # Read formulas from the uploaded database file
        db_path = save_tab2_file('db_file.db', db_content)
        print(f"DB file saved at: {db_path}")
        formula_ranges = read_formulas_from_db(db_path, charge, ppm)
        print(f"Formula ranges loaded from DB: {formula_ranges}")

        # Create a ZIP file to store results
        with zipfile.ZipFile(zip_filepath, 'w') as zf:
            for content, filename in zip(csv_contents, csv_filenames):
                print(f"Processing CSV file: {filename}")

                # Decode the uploaded CSV content
                content_type, content_string = content.split(',')
                decoded_csv = base64.b64decode(content_string)
                df = pd.read_csv(io.StringIO(decoded_csv.decode('utf-8')))
                print(f"CSV DataFrame loaded with shape: {df.shape}")

                # Calculate elemental composition for each row in the CSV
                df['elemental composition'] = df['row m/z'].apply(lambda mz: find_matching_formulas(mz, formula_ranges))
                print(f"Elemental compositions calculated for {filename}.")

                # Store results globally for later use
                global_results_composition[filename] = df

                # Modify output file name: "original_filename-mf.csv"
                output_filename = f"{os.path.splitext(filename)[0]}-mf.csv"
                output_file = os.path.join('temp', output_filename)

                # Save the processed CSV and add it to the ZIP file
                df.to_csv(output_file, index=False)
                print(f"Processed CSV saved at: {output_file}")
                zf.write(output_file, os.path.basename(output_file))
                print(f"Added {output_filename} to ZIP file.")

        print(f"ZIP file created at: {zip_filepath}")
        return (
            "Calculation done. Results saved to the specified path.",
            f'/download/{zip_filename}',
            [{'label': f"{os.path.splitext(name)[0]}-mf", 'value': f"{os.path.splitext(name)[0]}-mf"} for name in csv_filenames]
        )

    print("Conditions not met for calculation. No output generated.")
    return "", "", []


# Callback to display the results of elemental composition calculations in a data table on Tab 2
@app.callback(
    Output('mf-results-table', 'children'),
    Input('mf-file-dropdown', 'value')
)
def update_composition_results_table(selected_file):
    """
    Display the results of elemental composition calculations for the selected file.

    Parameters:
        selected_file: Filename selected from the dropdown.

    Returns:
        A Dash DataTable component displaying the results from the selected file.
    """
    print(f"Selected file for results display: {selected_file}")

    if selected_file and selected_file in global_results_composition:
        df = global_results_composition[selected_file]
        print(f"Displaying results for file: {selected_file}, DataFrame shape: {df.shape}")
        return dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in df.columns],
            style_cell={'fontFamily': 'Arial', 'textAlign': 'left'}
        )

    print("No results available for the selected file.")
    return "No results available"


# Server route to handle file downloads of calculated results in a ZIP file
@app.server.route('/download/<path:filename>')
def download_file(filename):
    """
    Route to download a file from the server.

    Parameters:
        filename: Name of the file to download.

    Returns:
        A response with the specified file to be downloaded as an attachment.
    """
    file_path = os.path.join('temp', filename)
    print(f"Download requested for file: {file_path}")

    if os.path.exists(file_path):
        print(f"File found, preparing download: {filename}")
        return send_file(file_path, as_attachment=True)
    else:
        print(f"File not found: {file_path}")
        return "File not found.", 404


# Function to read formulas from the database and calculate their m/z ranges
def read_formulas_from_db(db_path: str, charge: int, ppm: float) -> Dict[str, Tuple[float, float]]:
    """
    Read formulas from a database and calculate their m/z ranges based on charge and tolerance.

    Parameters:
        db_path: Path to the database file.
        charge: Charge state for m/z calculation.
        ppm: Tolerance level in PPM for mass range.

    Returns:
        A dictionary with formulas as keys and tuples of (min m/z, max m/z) as values.
    """
    print(f"Connecting to database at: {db_path}")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Query the formulas from the "compounds" table
        cursor.execute("SELECT formula FROM compounds")
        formulas = cursor.fetchall()
        print(f"Number of formulas retrieved from database: {len(formulas)}")
    except sqlite3.Error as e:
        print(f"Error querying database: {e}")
        raise e
    finally:
        conn.close()
        print(f"Database connection closed.")

    # Calculate m/z range for each formula
    formula_ranges = {}
    for formula_tuple in formulas:
        formula = formula_tuple[0]
        try:
            mz_range = mf_calculator_mz_range(formula, charge, ppm)
            formula_ranges[formula] = mz_range
            print(f"Calculated m/z range for formula '{formula}': min={mz_range[0]}, max={mz_range[1]}")
        except Exception as e:
            print(f"Error calculating m/z range for formula '{formula}': {e}")
    return formula_ranges


# Function to calculate the m/z range for a given formula, charge, and tolerance
def mf_calculator_mz_range(formula: str, charge: int, ppm: float) -> Tuple[float, float]:
    """
    Calculate the m/z range for a specific formula with the given charge and tolerance.

    Parameters:
        formula: Chemical formula for mass calculation.
        charge: Charge state for m/z calculation.
        ppm: Tolerance level in PPM for the range.

    Returns:
        A tuple containing the absolute minimum and maximum m/z values.
    """
    print(f"Calculating m/z range for formula: {formula}, charge: {charge}, ppm: {ppm}")
    try:
        mz = mass.calculate_mass(formula=formula, charge=charge)
        print(f"Calculated m/z for formula '{formula}': {mz}")
    except Exception as e:
        print(f"Error calculating m/z for formula '{formula}': {e}")
        raise e

    mz_min = mz - mz * ppm / 1e6
    mz_max = mz + mz * ppm / 1e6
    print(f"m/z range for formula '{formula}': min={mz_min}, max={mz_max}")
    return abs(mz_min), abs(mz_max)


# Function to find matching formulas within the m/z range for a given value
def find_matching_formulas(mz_value: float, formula_ranges: Dict[str, Tuple[float, float]]) -> str:
    """
    Find the best matching formula for a given m/z value within specified ranges.

    Parameters:
        mz_value: The m/z value for which to find a matching formula.
        formula_ranges: Dictionary of formulas and their m/z ranges.

    Returns:
        The formula that best matches the m/z value within the given ranges.
    """
    print(f"Finding best matching formula for m/z value: {mz_value}")
    min_ppm = float('inf')
    best_match = ''

    # Iterate through all formulas and their m/z ranges to find the closest match
    for formula, (mz_min, mz_max) in formula_ranges.items():
        if mz_min <= mz_value <= mz_max:
            mz = (mz_max + mz_min) / 2
            ppm = abs((mz_value - mz) / mz * 1e6)
            print(f"Formula '{formula}' matches m/z range: min={mz_min}, max={mz_max}, ppm difference: {ppm}")
            if ppm < min_ppm:
                min_ppm = ppm
                best_match = formula

    if best_match:
        print(f"Best matching formula: {best_match} with ppm difference: {min_ppm}")
    else:
        print("No matching formula found.")

    return best_match


# Utility function to save an uploaded file to a specified folder for Tab 2
def save_tab2_file(name, content, tab_folder='tab_2_uploads'):
    """
    Save an uploaded file to a specific folder for Tab 2.

    Parameters:
        name: Name of the file being saved.
        content: Base64-encoded content of the file.
        tab_folder: Name of the folder to save the file in (default is 'tab_2_uploads').

    Returns:
        The full path to the saved file.
    """
    print(f"Saving file '{name}' to folder '{tab_folder}'")
    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Created folder: {folder_path}")

    # Decode the content and save the file
    try:
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        file_path = os.path.join(folder_path, name)
        with open(file_path, 'wb') as f:
            f.write(decoded)
        print(f"File '{name}' saved successfully at: {file_path}")
    except Exception as e:
        print(f"Error saving file '{name}': {e}")
        raise e

    return file_path


"""
Tab 3 Related Functions for CSV Processing and Adduct Ion Calculation
"""
db_filenames = None  # Will store the uploaded database DataFrame
# Two global dictionaries are used to store known compound configuration information
tab3_typeI = {}
tab3_typeII = {}


@app.callback(
    Output('adduct-db-file-path', 'children'),
    Input('adduct-upload-db', 'contents'),
    State('adduct-upload-db', 'filename')
)
def update_input_db(input_db_content, input_db_filename):
    global db_filenames, tab3_typeI, tab3_typeII

    if input_db_content and input_db_filename:
        content_type, content_string = input_db_content.split(',')
        decoded = base64.b64decode(content_string)
        file_path = os.path.join(UPLOAD_DIR, input_db_filename)
        with open(file_path, 'wb') as f:
            f.write(decoded)

        # read DataFrame
        db_filenames = pd.read_csv(file_path)

        # After each database upload, reset and fill tab3_typeI / tab3_typeII
        tab3_typeI = {}
        tab3_typeII = {}

        # Assume the database CSV contains the following columns: ["Compound Name", "Type", "M", "M+NH3", "M+CH3NH2"]
        for _, row in db_filenames.iterrows():
            compound_name = row.get("Compound Name", "")
            compound_type = row.get("Type", "")
            m_formula = row.get("M", "")
            m_nh3_formula = row.get("M+NH3", "")
            m_ch3nh2_formula = row.get("M+CH3NH2", "")

            if compound_type == "D":
                tab3_typeI[compound_name] = [m_formula, m_nh3_formula, m_ch3nh2_formula]
            elif compound_type == "MD":
                tab3_typeII[compound_name] = [m_formula, m_nh3_formula, m_ch3nh2_formula]

        return f"Uploaded Database File: {input_db_filename}"
    return "No database file uploaded."


@app.callback(
    Output('adduct-csv-file-path', 'children'),
    Input('adduct-upload-csv', 'contents'),
    State('adduct-upload-csv', 'filename')
)
def update_adduct_csv_file_path(csv_contents, csv_filenames):
    """
    Update frontend to show: which CSVs were uploaded
    """
    if csv_contents and csv_filenames:
        return f"Uploaded {len(csv_filenames)} CSV file(s): {', '.join(csv_filenames)}"
    return "No CSV file uploaded"


@app.callback(
    [
        Output('adduct-calculation-output', 'children'),
        Output('adduct-download-link', 'href'),
        Output('adduct-file-dropdown', 'options')
    ],
    Input('adduct-calculate-button', 'n_clicks'),
    State('adduct-upload-csv', 'contents'),
    State('adduct-upload-csv', 'filename'),
    State('adduct-input-charge', 'value'),
    State('adduct-input-ppm', 'value')
)
def calculate_adduct_ions(n_clicks, csv_contents, csv_filenames, charge, ppm):
    """
    Parse the uploaded CSV content
    Perform the following on each CSV:
    - Simple calculation of Adduct Ions (process_adduct_ions)
    - Compound identification/annotation/water loss (NH3/CH3NH2) and other logic (process_annotation_dataframe)
    Write the results to the temp folder and package them into adduct_result.zip
    Return to the download link, drop-down menu options
    """
    if n_clicks > 0 and csv_contents and charge and ppm:
        # The fixed packaged file name is "adduct_result.zip"
        zip_filename = 'adduct_result.zip'
        zip_filepath = os.path.join('temp', zip_filename)
        # Create ZIP and write processed results to CSV
        with zipfile.ZipFile(zip_filepath, 'w') as zf:
            for content, filename in zip(csv_contents, csv_filenames):
                content_type, content_string = content.split(',')
                decoded_csv = base64.b64decode(content_string)
                df = pd.read_csv(io.StringIO(decoded_csv.decode('utf-8')))
                # Re-execute compound identification/annotation logic
                df = process_annotation_dataframe(df)
                # Store the processing results in the global dictionary
                global_results_adduct[filename] = df
                # Construct a new file name with `-ad`
                base_name, extension = os.path.splitext(filename)
                annotated_name = f"{base_name}-ad{extension}"
                # Write the results to a temporary file and then package it
                output_file = os.path.join('temp', annotated_name)
                df.to_csv(output_file, index=False, encoding='utf-8')
                # Write the new file name to the ZIP (internal name is annotated_name)
                zf.write(output_file, annotated_name)
        # The front-end displays that the processing is complete and returns a download link (such as /download/adduct_result.zip)
        return (
            "Calculation done. Results saved.",
            f'/download/{zip_filename}',
            [{'label': name, 'value': name} for name in csv_filenames]
        )
    return "", "", []


@app.callback(
    Output('adduct-results-table', 'children'),
    Input('adduct-file-dropdown', 'value')
)
def update_adduct_results_table(selected_file):
    if selected_file and selected_file in global_results_adduct:
        df = global_results_adduct[selected_file]
        return dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in df.columns],
            style_cell={'fontFamily': 'Arial', 'textAlign': 'left'}
        )
    return "No results available"


@app.server.route('/download/<path:filename>', endpoint='unique_download_file')
def download_file(filename):
    """
    Allow to access zip files in the "temp" folder
    """
    file_path = os.path.join('temp', filename)
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return "File not found.", 404


# Calculate the m/z (mass-to-charge) range for a given molecular formula, charge state, and ppm (parts per million) tolerance.
def calculate_mz_range(formula: str, charge: int = 1, tolerance_ppm: float = 5.0) -> Tuple[float, float]:
    """
    Parameters:
    formula (str): The molecular formula of the compound.
    charge (int): The charge state of the ion. Default is 1, which is typical for charged ions.
    tolerance_ppm (float): The mass accuracy tolerance in parts per million (ppm).

    Returns:
    tuple[float, float]: A tuple containing the lower and upper bounds of the m/z range.
                        The first element is the minimum m/z value, and the second element is the maximum m/z value.
    """
    # Calculate the nominal m/z value for the given formula and charge state using the mass.calculate_mass function from pyteomics.
    mz = mass.calculate_mass(formula=formula, charge=charge)
    # Calculate the m/z tolerance based on the specified ppm. This determines how much the m/z value can vary.
    mz_tolerance = mz * tolerance_ppm / 1e6
    # Return the m/z range as a tuple, where the first element is the lower bound (mz - tolerance)
    # and the second is the upper bound (mz + tolerance).
    return mz - mz_tolerance, mz + mz_tolerance


def identify_compound(row):
    """
    Parameters:
    row (dict): A dictionary representing a row of data from mass spectrometry. It should include 'type' to indicate the
                type of natural compounds (ex.'D' or 'MD') and 'row m/z' for the mass-to-charge ratio to be identified.

    Returns:
    tuple[str, str]: A tuple where the first element indicates if the compound is known or undescribed,
                    and the second element lists the identified compounds or 'none' if no matches are found.
    """
    if row.get('type') == 'D':
        known_compounds = tab3_typeI
    elif row.get('type') == 'MD':
        known_compounds = tab3_typeII
    else:
        known_compounds = {}

    identifications = []
    for compound, formulas in known_compounds.items():
        for formula in formulas:
            mz_min, mz_max = calculate_mz_range(formula, charge=1, tolerance_ppm=5.0)
            if mz_min <= row['row m/z'] <= mz_max:
                identifications.append(compound)

    if identifications:
        return 'possibly known', ', '.join(set(identifications))
    else:
        return 'possibly undescribed', 'none'


def calculate_formula_difference(formula1: str, formula2: str) -> dict:
    """
    Parameters:
    formula1 (str): The chemical formula of the first compound (e.g., "C6H12O6").
    formula2 (str): The chemical formula of the second compound (e.g., "C5H10O5").

    Returns:
    dict: A dictionary where keys are elemental ('C', 'H', 'O', 'N') and values are the differences in
        element counts between the two formulas. Positive values indicate a greater number in formula1;
        negative values indicate a greater number in formula2.
    """
    # Create composition objects for each formula using the Composition class of pyteomics.
    # Provide a dictionary-like interface where elements are keys and their counts are values.
    comp1 = mass.Composition(formula1)
    comp2 = mass.Composition(formula2)
    # Calculate molecular weights for both compositions to determine which is greater.
    mw1 = mass.calculate_mass(composition=comp1)
    mw2 = mass.calculate_mass(composition=comp2)
    # Determine which formula has the greater molecular weight.
    # Calculate the difference in elemental counts for common elements ('C', 'H', 'O', 'N').
    if mw1 > mw2:
        # If the first formula has a greater molecular weight,
        # compute the difference by subtracting elements of the second from the first.
        diff = {elem: comp1.get(elem, 0) - comp2.get(elem, 0) for elem in ('C', 'H', 'O', 'N')}
    else:
        # If the second formula has a greater molecular weight or if they are equal,
        # compute the difference by subtracting elements of the first from the second.
        diff = {elem: comp2.get(elem, 0) - comp1.get(elem, 0) for elem in ('C', 'H', 'O', 'N')}
    # Return the dictionary containing the differences in elemental counts.
    return diff


# Determine if the elemental difference between two formulas corresponds to the molecule H2O (water), NH3 (ammonia) or CH3NH2.
def h2o_nh3_ch3nh2(diff: dict) -> str:
    """
    Parameters:
    diff (dict): A dictionary representing the difference in elemental composition. Keys are elements ('C', 'H', 'O', 'N')
                and values are the differences in the count of each element.

    Returns:
    str: "H2O" if the difference corresponds to a water molecule, "NH3" if it corresponds to an ammonia molecule,
        or None if the difference does not match either.
    """
    # Check if the difference in elemental composition exactly matches that of a water molecule.
    if diff == {'C': 0, 'H': 2, 'O': 1, 'N': 0}:
        return "H2O"  # Return "H2O" if the difference matches a water molecule.
    # Check if the difference in elemental composition exactly matches that of an ammonia molecule.
    elif diff == {'C': 0, 'H': 3, 'O': 0, 'N': 1}:
        return "NH3"  # Return "NH3" if the difference matches an ammonia molecule.
    elif diff == {'C': 1, 'H': 5, 'O': 0, 'N': 1}:
        return "CH3NH2"  # Return "CH3NH2" if the difference matches an CH3NH2 molecule.
    # Return None if the elemental difference does not correspond to either H2O or NH3.
    return None


# Add annotations to the DataFrame, identifying potential adduct ions and water/ammonia loss.
def add_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters:
    df (pd.DataFrame): DataFrame containing MS data with columns for 'elemental composition',
                       'row retention time', and 'row m/z'.

    Returns:
    pd.DataFrame: The updated DataFrame with added annotations for adduct ions.
    """

    # Annotate initial adduct ion types based on presence of nitrogen in the elemental composition.
    def determine_adduct_ion(composition):
        if 'N' not in composition:
            return '[M+H]+'
        else:
            return '[M+NH4]+'

    # Ensure 'elemental composition' column has no NaNs and is of string type
    df['elemental composition'] = df['elemental composition'].fillna('').astype(str)

    # Apply the determine_adduct_ion function safely
    df['adduct ion'] = df['elemental composition'].apply(determine_adduct_ion)

    # Identify duplicated retention times to check for potential molecular relationships.
    duplicated_rts = df['row retention time'].duplicated(keep=False)

    # Loop over each unique duplicated retention time to process potential relationships.
    for rt in df[duplicated_rts]['row retention time'].unique():
        indices = df.index[df['row retention time'] == rt].tolist()
        for i in range(len(indices) - 1):
            for j in range(i + 1, len(indices)):
                diff = calculate_formula_difference(
                    df.at[indices[i], 'elemental composition'],
                    df.at[indices[j], 'elemental composition']
                )
                match_result = h2o_nh3_ch3nh2(diff)

                if match_result == "H2O":

                    if (mass.Composition(df.at[indices[i], 'elemental composition']).get('O', 0)
                            < mass.Composition(df.at[indices[j], 'elemental composition']).get('O', 0)):
                        df.at[indices[i], 'adduct ion'] = '[M+H-H2O]+'
                    else:
                        df.at[indices[j], 'adduct ion'] = '[M+H-H2O]+'

                elif match_result == "NH3":
                    if df.at[indices[i], 'row m/z'] > df.at[indices[j], 'row m/z']:
                        df.at[indices[i], 'adduct ion'] = '[M+NH4]+'
                    else:
                        df.at[indices[j], 'adduct ion'] = '[M+NH4]+'

                elif match_result == "CH3NH2":
                    if df.at[indices[i], 'row m/z'] > df.at[indices[j], 'row m/z']:
                        df.at[indices[i], 'adduct ion'] = '[M+CH3NH2+H]+'
                    else:
                        df.at[indices[j], 'adduct ion'] = '[M+CH3NH2+H]+'

    # Return the annotated DataFrame.
    return df


# Compare pairs of peaks within the same retention time and mark them if composition difference corresponds to H2O or NH3.
def compare_peaks_and_mark(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters:
    df (pandas.DataFrame): DataFrame containing columns 'row retention time', 'type', 'elemental composition', and 'row ID'.

    Returns:
    pandas.DataFrame: Modified DataFrame with two new columns: 'same peak' and 'C part'.
    """
    same_peak_dict = {}  # Dictionary to store related peaks for each row
    grouped = df.groupby('row retention time')      # Group the DataFrame by 'row retention time'

    # Iterate over each group
    for name, group in grouped:
        if len(group) > 1:  # Only process groups with more than one row
            ids = group.index  # Get the indices of the rows in the group
            # Compare each pair of rows within the group
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    # Extract relevant data for the pair of rows
                    formula1 = group.at[ids[i], 'elemental composition']
                    formula2 = group.at[ids[j], 'elemental composition']
                    # Calculate the formula difference and determine if it matches H2O or NH3
                    diff = calculate_formula_difference(formula1, formula2)
                    result = h2o_nh3_ch3nh2(diff)
                    # If the result is H2O or NH3, update the same_peak_dict
                    if result in ["H2O", "NH3", "CH3NH2"]:
                        if ids[i] in same_peak_dict:
                            same_peak_dict[ids[i]].add(group.at[ids[j], 'row ID'])
                        else:
                            same_peak_dict[ids[i]] = {group.at[ids[i], 'row ID'], group.at[ids[j], 'row ID']}
                        if ids[j] in same_peak_dict:
                            same_peak_dict[ids[j]].add(group.at[ids[i], 'row ID'])
                        else:
                            same_peak_dict[ids[j]] = {group.at[ids[i], 'row ID'], group.at[ids[j], 'row ID']}

    # Convert sets in same_peak_dict to sorted, comma-separated strings
    for key, value in same_peak_dict.items():
        same_peak_dict[key] = ', '.join(sorted(str(v) for v in value))
    # Ensure all indices have an entry in same_peak_dict, defaulting to 'none'
    for idx in df.index:
        if idx not in same_peak_dict:
            same_peak_dict[idx] = 'none'
    df['same peak'] = pd.Series(same_peak_dict)
    return df


# Processes a DataFrame with mass spectrometry data by adding identification results, annotations, and organizing the data.
def process_annotation_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters:
    df (pd.DataFrame): The input DataFrame containing mass spectrometry data, expected to include 'row retention time',
    'elemental composition', and 'peak area' columns among others.

    Returns:
    pd.DataFrame: The processed DataFrame with added identifications, annotations, and sorted by peak area.
    """
    # Add identification results for each row. Each row is passed to the identify_compound function,
    # and the results are expanded into two new columns: 'identification' and 'same composition'.
    df['identification'], df['same composition'] = zip(*df.apply(identify_compound, axis=1))
    # Add annotations for Adduct Ions based on elemental compositions and relationships identified by retention times.
    df = add_annotations(df)
    # Round the 'row retention time' values to two decimal places to maintain uniformity and improve readability.
    df = compare_peaks_and_mark(df)
    if 'row retention time' in df.columns:
        df['row retention time'] = df['row retention time'].round(1)
    return df


"""Tab 4 related function"""
# Callback to update the MGF file path upon upload
@app.callback(
    Output('neutralloss-mgf-file-path', 'children'),
    Input('neutralloss-upload-mgf', 'contents'),
    State('neutralloss-upload-mgf', 'filename')
)
def update_neutralloss_mgf_file_path(mgf_contents, mgf_filenames):
    print("MGF upload triggered")
    if mgf_contents and mgf_filenames:
        print(f"MGF file names: {mgf_filenames}")
        # Save each uploaded MGF file
        mgf_paths = [save_tab4_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
        mgf_path_text = f"Uploaded MGF files: {', '.join(mgf_filenames)}"
    else:
        print("No MGF file uploaded")
        mgf_path_text = "No MGF file uploaded"
    return mgf_path_text

# Callback to update the CSV file path upon upload
@app.callback(
    Output('neutralloss-csv-file-path', 'children'),
    Input('neutralloss-upload-csv', 'contents'),
    State('neutralloss-upload-csv', 'filename')
)
def update_neutralloss_csv_file_path(csv_contents, csv_filenames):
    print("CSV upload triggered")
    if csv_contents and csv_filenames:
        print(f"CSV file names: {csv_filenames}")
        # Save each uploaded CSV file
        csv_paths = [save_tab4_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

# Callback to update the database file path upon upload
@app.callback(
    Output('neutralloss-db-file-path', 'children'),
    Input('neutralloss-upload-db', 'contents'),
    State('neutralloss-upload-db', 'filename')
)
def update_neutralloss_db_file_path(db_content, db_filename):
    print("DB upload triggered")
    if db_content and db_filename:
        print(f"DB file name: {db_filename}")
        # Save the uploaded database file
        db_path = save_tab4_file(db_filename, db_content)
        db_path_text = f"Uploaded DB file: {db_filename}"
    else:
        print("No DB file uploaded")
        db_path_text = "No DB file uploaded"
    return [db_path_text]

# Callback to add or remove target loss input areas dynamically
@app.callback(
    Output('neutralloss-div', 'children'),
    [Input('add-new-targetloss', 'n_clicks'),
     Input({'type': 'remove-button-neutralloss', 'index': ALL}, 'n_clicks')],
    [State('neutralloss-div', 'children')]
)
def neutralloss_modify_types(add_clicks, remove_clicks, children):
    print(f"Add new target loss clicked: {add_clicks} times")
    if add_clicks > 0:
        new_index = len(children) - 1 # Update the new index to avoid overwriting existing items
        print(f"Adding new input area at index: {new_index}")
        new_input_area = html.Div(id={'type': 'input-area', 'index': new_index + 1}, children=[
            html.Div([
                html.Label(f'Loss Name {new_index + 1}:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(
                    id={'type': 'input-type-name', 'index': new_index + 1},
                    type='text',
                    placeholder=f'e.g., C3 loss',
                    style={'width': '10%', 'fontFamily': 'Arial', 'marginLeft': '10px'}
                )
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '20px'}),
            html.H4(f'Target Loss {new_index + 1} :', style={'fontFamily': 'Arial'}),
            dcc.Textarea(
                id={'type': 'input-type-targetloss-formulas', 'index': new_index + 1},
                placeholder='Input target loss formulas',
                style={'width': '50%', 'height': '50px', 'fontFamily': 'Arial'}
            ),
            html.Div([
                html.Label('Carbon Range 1 From :',
                           style={'fontFamily': 'Arial', 'fontWeight': 'bold', 'marginRight': '10px'}),
                dcc.Input(id={'type': 'carbon-range-max1', 'index': new_index + 1}, type='text', placeholder='C20',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),
                dcc.Input(id={'type': 'carbon-range-min1', 'index': new_index + 1}, type='text', placeholder='C17',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'}),
            html.Div([
                html.Label('Carbon Range 2 From :',
                           style={'fontFamily': 'Arial', 'fontWeight': 'bold', 'marginRight': '10px'}),
                dcc.Input(id={'type': 'carbon-range-max2', 'index': new_index + 1}, type='text', placeholder='C20',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),
                dcc.Input(id={'type': 'carbon-range-min2', 'index': new_index + 1}, type='text', placeholder='C17',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'}),
            html.Div([
                html.Label('Carbon Range 3 From :',
                           style={'fontFamily': 'Arial', 'fontWeight': 'bold', 'marginRight': '10px'}),
                dcc.Input(id={'type': 'carbon-range-max3', 'index': new_index + 1}, type='text', placeholder='C30',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),
                dcc.Input(id={'type': 'carbon-range-min3', 'index': new_index + 1}, type='text', placeholder='C27',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Button('✖', id={'type': 'remove-button-neutralloss', 'index': new_index + 1},
                            style={'marginLeft': '10px', 'color': 'red'})
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'}),
        ])
        children.insert(-1, new_input_area)

    # Context for determining if a remove button was clicked
    ctx = dash.callback_context
    if ctx.triggered:
        triggered = ctx.triggered[0]['prop_id'].split('.')[0]
        if 'remove-button-neutralloss' in triggered:
            button_index = eval(triggered)['index']
            print(f"Removing input area at index: {button_index}")
            # Filter out the input area that matches the clicked remove button's index
            children = [
                child for child in children
                if isinstance(child, dict) and
                'props' in child and
                isinstance(child['props'].get('id', {}), dict) and
                child['props']['id'].get('index') != button_index and
                child['props']['id'] != 'add-new-targetloss'
            ]
    return children

# Callback to generate a downloadable ZIP file after extraction
@app.callback(
    Output('neutralloss-download-link', 'href'),
    Input('run-neutralloss-extraction', 'n_clicks'),
    State('neutralloss-upload-csv', 'filename'),
    State('neutralloss-upload-csv', 'contents')
)
def neutralloss_generate_zip_for_download(n_clicks, csv_filenames, csv_contents):
    print("Generating ZIP file for download")
    if n_clicks is None or n_clicks == 0:
        print("No clicks to generate ZIP file")
        return ""

    zip_filename = os.path.join(UPLOAD_DIR, 'neutralloss_extraction_results.zip')
    # Create ZIP file for all processed CSV files
    with zipfile.ZipFile(zip_filename, 'w') as zf:
        for csv_filename in csv_filenames:
            processed_file = csv_filename.replace(".csv", "-nl.csv")
            file_path = os.path.join(UPLOAD_DIR, 'tab_4_uploads', processed_file)
            print(f"Adding {processed_file} to ZIP")
            zf.write(file_path, processed_file)

    return f'/download/{os.path.basename(zip_filename)}'

# Callback to update the file dropdown options after extraction
@app.callback(
    Output('neutralloss-extract-file-dropdown', 'options'),
    Input('run-neutralloss-extraction', 'n_clicks'),
    State('neutralloss-upload-csv', 'filename')
)
def neutralloss_update_file_dropdown(n_clicks, csv_filenames):
    print("Updating file dropdown options")
    if n_clicks is None or n_clicks == 0:
        print("No clicks to update dropdown")
        return []

    return [{'label': filename.replace(".csv", "-nl.csv"), 'value': filename.replace(".csv", "-nl.csv")}
            for filename in csv_filenames]

# Callback to display the selected CSV file results in a table
@app.callback(
    Output('neutralloss-extract-results-table', 'children'),
    Input('neutralloss-extract-file-dropdown', 'value')
)
def neutralloss_display_selected_file_results(selected_file):
    if not selected_file:
        print("No file selected for results display")
        return "Please select a file to view results."

    file_path = os.path.join(UPLOAD_DIR, 'tab_4_uploads', selected_file)
    print(f"Loading results from file: {file_path}")
    df = pd.read_csv(file_path)

    return dash_table.DataTable(
        data=df.to_dict('records'),
        columns=[{'name': i, 'id': i} for i in df.columns],
        style_header={'fontFamily': 'Arial', 'fontWeight': 'bold'},
        style_cell={'fontFamily': 'Arial', 'textAlign': 'left'}
    )

# Server route for downloading files
@app.server.route('/download/<path:filename>')
def neutralloss_download(filename):
    print(f"Downloading file: {filename}")
    return send_file(os.path.join(UPLOAD_DIR, filename), as_attachment=True)

# Callback to run the neutral loss extraction process when the button is clicked
@app.callback(
    Output('neutralloss-extract-output', 'children'),
    Input('run-neutralloss-extraction', 'n_clicks'),
    State('neutralloss-upload-csv', 'filename'),
    State('neutralloss-upload-csv', 'contents'),
    State('neutralloss-upload-mgf', 'filename'),
    State('neutralloss-upload-mgf', 'contents'),
    State('neutralloss-upload-db', 'filename'),
    State('neutralloss-upload-db', 'contents'),
    State('neutralloss-input-charge', 'value'),
    State('neutralloss-input-ppm', 'value'),
    State({'type': 'input-type-targetloss-formulas', 'index': ALL}, 'value'),
    State({'type': 'input-type-name', 'index': ALL}, 'value'),
    State({'type': 'carbon-range-max1', 'index': ALL}, 'value'),
    State({'type': 'carbon-range-min1', 'index': ALL}, 'value'),
    State({'type': 'carbon-range-max2', 'index': ALL}, 'value'),
    State({'type': 'carbon-range-min2', 'index': ALL}, 'value'),
    State({'type': 'carbon-range-max3', 'index': ALL}, 'value'),
    State({'type': 'carbon-range-min3', 'index': ALL}, 'value')
)
def run_neutral_loss_extraction(n_clicks, csv_filenames, csv_contents, mgf_filenames, mgf_contents, db_filename,
                                db_content, charge, ppm, target_losses, loss_names,
                                carbon_range_max1, carbon_range_min1, carbon_range_max2, carbon_range_min2,
                                carbon_range_max3, carbon_range_min3):
    """
    Execute neutral loss extraction based on user-provided parameters and uploaded files.

    Parameters:
        n_clicks: Number of clicks on the "Run Extraction" button.
        csv_filenames: List of CSV filenames uploaded.
        csv_contents: Contents of the uploaded CSV files.
        mgf_filenames: List of MGF filenames uploaded.
        mgf_contents: Contents of the uploaded MGF files.
        db_filename: Name of the uploaded database file.
        db_content: Content of the uploaded database file.
        charge: Charge state for mass calculation.
        ppm: Tolerance in parts per million for mass range matching.
        target_losses: List of target loss formulas provided by the user.
        loss_names: List of names corresponding to each target loss.
        carbon_range_max1, carbon_range_min1, carbon_range_max2, carbon_range_min2, carbon_range_max3,
        carbon_range_min3: Carbon range specifications.

    Returns:
        Status message indicating the success or failure of the extraction process.
    """
    print(f"Run button clicked {n_clicks} times")

    # Check if the button has been clicked
    if n_clicks is None or n_clicks == 0:
        print("Button not clicked yet")
        return "Click the button to start the extraction."

    # Ensure that all required files are uploaded
    if not (csv_filenames and mgf_filenames and db_filename):
        print("Missing required files")
        return "Please upload the required files (CSV, MGF, and DB)."

    # Save uploaded files to the appropriate folder
    csv_paths = [save_tab4_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
    mgf_paths = [save_tab4_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
    db_path = save_tab4_file(db_filename, db_content)

    print(f"CSV files saved at: {csv_paths}")
    print(f"MGF files saved at: {mgf_paths}")
    print(f"DB file saved at: {db_path}")

    # Check if target loss formulas are provided
    if not target_losses or all(x is None or x == "" for x in target_losses):
        print("No target losses provided")
        return "Please provide at least one target loss formula."

    # Parse and clean the target loss formulas
    parsed_target_losses = []
    for tl in target_losses:
        if tl:
            formulas = [x.strip() for x in tl.split(',') if x.strip()]
            parsed_target_losses.append(formulas)
        else:
            parsed_target_losses.append([])

    # Set default values for charge and ppm if they are not provided
    charge = charge or 1
    ppm = ppm or 5

    # Process each pair of CSV and MGF files
    for csv_path, mgf_path in zip(csv_paths, mgf_paths):
        output_csv_path = csv_path.replace(".csv", "-nl.csv")
        try:
            neutralloss_process_files(
                input_path_csv=csv_path,
                input_path_mgf=mgf_path,
                output_path_csv=output_csv_path,
                db_path=db_path,
                target_losses=parsed_target_losses,
                charge=charge or 1,
                ppm=ppm or 5,
                loss_names=loss_names,
                carbon_range_max1=carbon_range_max1,
                carbon_range_min1=carbon_range_min1,
                carbon_range_max2=carbon_range_max2,
                carbon_range_min2=carbon_range_min2,
                carbon_range_max3 = carbon_range_max3,
                carbon_range_min3 = carbon_range_min3
            )
        except Exception as e:
            return f"Error processing files: {str(e)}"

    return f"Extraction complete for {len(csv_paths)} CSV file(s)."


def neutralloss_process_files(input_path_csv, input_path_mgf, output_path_csv, db_path,
                              target_losses, charge, ppm,
                              loss_names, carbon_range_max1, carbon_range_min1,
                              carbon_range_max2, carbon_range_min2,
                              carbon_range_max3, carbon_range_min3):
    """
    Process an input CSV and MGF file to extract neutral-loss matches:
      1. Read the CSV and index by 'row ID', cleaning titles.
      2. Load the formula m/z ranges from the SQLite DB (with given charge and ppm).
      3. Read all spectra from the MGF into a dictionary keyed by title.
      4. Initialize output columns for each loss group in the DataFrame.
      5. Parse the user-provided carbon ranges for each loss group.
      6. For each spectrum (peak):
         a. Compute Top-3 matching ion pairs per loss group based on score.
         b. Deduplicate by loss formula (keep only the highest-score pair per formula).
         c. Expand all groups’ Top-K lists via Cartesian product so each combination
            becomes its own output row.
         d. Copy all original columns (including the computed row m/z) into each row.
      7. Write the exploded results to the output CSV.
    """
    # 1) Read CSV and set 'row ID' as index
    df = pd.read_csv(input_path_csv)
    if 'row ID' not in df.columns:
        raise ValueError("'row ID' column not found in input CSV.")
    df.set_index('row ID', inplace=True)
    df.index = df.index.map(clean_title)

    # 2) Load formula-to-m/z ranges from the neutral-loss database
    formula_ranges = read_formulas_from_NLdb(db_path, charge, ppm)

    # 3) Read MGF spectra into a dict for quick lookup
    spectra_dict = {}
    with mgf.read(input_path_mgf) as spectra:
        for spectrum in spectra:
            title = clean_title(spectrum['params']['title'])
            spectra_dict[title] = spectrum

    # 4) Initialize neutral-loss columns to 'none'
    for ln in loss_names:
        df[ln] = 'none'
        df[ln + '_ion_pair'] = 'none'
        df[ln + '_intensity_pair'] = 'none'
        df[ln + '_score'] = 'none'

    # 5) Parse all carbon-range inputs into integer tuples
    parsed_carbon_ranges = []
    for i in range(len(loss_names)):
        try:
            parsed_carbon_ranges.append((
                parse_carbon_range(carbon_range_max1[i]),
                parse_carbon_range(carbon_range_min1[i]),
                parse_carbon_range(carbon_range_max2[i]),
                parse_carbon_range(carbon_range_min2[i]),
                parse_carbon_range(carbon_range_max3[i]),
                parse_carbon_range(carbon_range_min3[i]),
            ))
        except Exception as e:
            raise ValueError(f"Invalid carbon range in group {i+1}: {e}")

    # Keep a deep copy of the original DataFrame for row-by-row expansion
    original_df = df.copy(deep=True)

    output_records = []

    # 6) Process each peak (row) in the original data
    for row_id, orig_row in original_df.iterrows():
        if row_id not in spectra_dict:
            continue

        spectrum = spectra_dict[row_id]
        mz_value = spectrum['params']['pepmass'][0]
        product_ions = spectrum['m/z array']
        intensities  = spectrum['intensity array']

        per_group_matches = []

        # 6a) Collect Top-3 matches per loss group, then dedupe by formula
        for idx, (loss_group, ln) in enumerate(zip(target_losses, loss_names)):
            # If no formulas provided, placeholder a single None entry
            if not loss_group:
                per_group_matches.append([(None, None, None, None, None)])
                continue

            cmax1, cmin1, cmax2, cmin2, cmax3, cmin3 = parsed_carbon_ranges[idx]

            # Build a list of all (index, matched_formula)
            formula_info = [
                (i, find_matching_NLformulas(mz, formula_ranges))
                for i, mz in enumerate(product_ions)
            ]

            # 6a.i) Original matching logic for three carbon-range pairs
            matching_pairs = []
            # First carbon range
            fmax1 = [(i, f) for i, f in formula_info
                     if f != 'none' and mass.Composition(f).get('C',0) == cmax1]
            fmin1 = [(i, f) for i, f in formula_info
                     if f != 'none' and mass.Composition(f).get('C',0) == cmin1]
            for i1, f1 in fmax1:
                for i2, f2 in fmin1:
                    if check_specific_target_loss(specific_formula_difference(f1, f2), loss_group):
                        matching_pairs.append((f1, f2, intensities[i1], intensities[i2]))

            # Second carbon range
            fmax2 = [(i, f) for i, f in formula_info
                     if f != 'none' and mass.Composition(f).get('C',0) == cmax2]
            fmin2 = [(i, f) for i, f in formula_info
                     if f != 'none' and mass.Composition(f).get('C',0) == cmin2]
            for i1, f1 in fmax2:
                for i2, f2 in fmin2:
                    if check_specific_target_loss(specific_formula_difference(f1, f2), loss_group):
                        matching_pairs.append((f1, f2, intensities[i1], intensities[i2]))

            # Third carbon range
            fmax3 = [(i, f) for i, f in formula_info
                     if f != 'none' and mass.Composition(f).get('C',0) == cmax3]
            fmin3 = [(i, f) for i, f in formula_info
                     if f != 'none' and mass.Composition(f).get('C',0) == cmin3]
            for i1, f1 in fmax3:
                for i2, f2 in fmin3:
                    if check_specific_target_loss(specific_formula_difference(f1, f2), loss_group):
                        matching_pairs.append((f1, f2, intensities[i1], intensities[i2]))

            # 6a.ii) Score each pair, sort, take Top-3
            if matching_pairs:
                scored = [(f1, f2, i1, i2, compute_score(i1, i2))
                          for f1, f2, i1, i2 in matching_pairs]
                scored.sort(key=lambda x: x[4], reverse=True)
                top3 = scored[:3]

                # 6a.iii) Deduplicate by loss formula: keep highest-score per unique formula
                deduped = []
                seen = set()
                for f1, f2, i1, i2, score in top3:
                    formula_str = dict_to_formula(specific_formula_difference(f1, f2))
                    if formula_str not in seen:
                        seen.add(formula_str)
                        deduped.append((f1, f2, i1, i2, score))
                per_group_matches.append(deduped)
            else:
                per_group_matches.append([(None, None, None, None, None)])

        # 6b) Perform Cartesian product across all groups to expand rows
        for combo in itertools.product(*per_group_matches):
            rec = orig_row.to_dict()
            rec['row ID']   = row_id
            rec['row m/z']  = mz_value
            for ln, (f1, f2, i1, i2, score) in zip(loss_names, combo):
                if f1 is None:
                    rec[ln]                    = 'none'
                    rec[ln + '_ion_pair']      = 'none'
                    rec[ln + '_intensity_pair']= 'none'
                    rec[ln + '_score']         = 'none'
                else:
                    # Write back selected loss formula, ion pair, intensity pair, and score
                    rec[ln]                    = dict_to_formula(specific_formula_difference(f1, f2))
                    rec[ln + '_ion_pair']      = f"{f1} / {f2}"
                    rec[ln + '_intensity_pair']= f"{i1} / {i2}"
                    rec[ln + '_score']         = score
            output_records.append(rec)

    # 7) Convert to DataFrame
    out_df = pd.DataFrame(output_records)

    # 8) Reorder columns: original columns first, then loss columns
    orig_cols = ['row ID'] + [col for col in original_df.columns if col != 'row ID']
    loss_cols = sum([[ln, ln + '_ion_pair', ln + '_intensity_pair', ln + '_score'] for ln in loss_names], [])

    # Combine and remove any duplicates just in case
    final_cols = list(dict.fromkeys(orig_cols + loss_cols))
    out_df[final_cols].to_csv(output_path_csv, index=False)


def neutralloss_mz_range(formula: str, charge: int, ppm: float) -> Tuple[float, float]:
    mz = mass.calculate_mass(formula=formula, charge=charge)
    mz_min = mz - mz * ppm / 1e6
    mz_max = mz + mz * ppm / 1e6
    return abs(mz_min), abs(mz_max)

def read_formulas_from_NLdb(db_path: str, charge: int, ppm: float) -> Dict[str, Tuple[float, float]]:
    print(f"Reading formulas from DB: {db_path}")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT formula FROM compounds")
    formulas = cursor.fetchall()
    conn.close()
    print(f"Found {len(formulas)} formulas in DB")
    return {formula[0]: neutralloss_mz_range(formula[0], charge, ppm) for formula in formulas}

def find_matching_NLformulas(mz_value: float, formula_ranges: dict) -> str:
    matching_formulas = [
        formula for formula, (mz_min, mz_max) in formula_ranges.items()
        if mz_min <= mz_value <= mz_max
    ]
    print(f"Matching formulas for mz {mz_value}: {matching_formulas}")
    return matching_formulas[0] if len(matching_formulas) == 1 else 'none'

def specific_formula_difference(formula1: str, formula2: str) -> Dict[str, int]:
    comp1 = mass.Composition(formula1)
    comp2 = mass.Composition(formula2)
    differences = {
        'C': abs(comp2.get('C', 0) - comp1.get('C', 0)),
        'H': abs(comp2.get('H', 0) - comp1.get('H', 0)),
        'O': abs(comp2.get('O', 0) - comp1.get('O', 0))
    }
    print(f"Difference between {formula1} and {formula2}: {differences}")
    return differences

def check_specific_target_loss(diff_comp: dict, target_losses: list) -> bool:
    for target_loss in target_losses:
        target_elements = dict(mass.Composition(target_loss))
        print(f"Checking target loss - Differences: {diff_comp}, Target: {target_elements}")
        if all(diff_comp.get(k, 0) == target_elements.get(k, 0) for k in target_elements):
            return True
    return False

def compute_score(intensity1, intensity2):
    IA = (intensity1 + intensity2) / 2
    ID = abs(intensity1 - intensity2)
    IDF = 1 / (1 + (ID / IA)) if IA != 0 else 0
    score = IA * IDF / 100000
    return round(score, 4)

def clean_title(title):
    return str(title).strip()

def parse_carbon_range(carbon_range_str):
    match = re.match(r'C(\d+)', carbon_range_str)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Invalid carbon range input: {carbon_range_str}")

def dict_to_formula(target_elements):
    formula = ''
    for element, count in target_elements.items():
        if count == 1:
            formula += f"{element}"
        else:
            formula += f"{element}{int(count)}"
    return formula

def save_tab4_file(name, content, tab_folder='tab_4_uploads'):
    print(f"Saving file: {name}")
    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)
    print(f"File saved at: {file_path}")
    with open(file_path, 'wb') as f:
        f.write(decoded)
    return file_path


"""Tab 5 related function"""
# Callback for updating the path of the uploaded MGF file
@app.callback(
    Output('feature-ion-mgf-file-path', 'children'),
    Input('feature-ion-upload-mgf', 'contents'),
    State('feature-ion-upload-mgf', 'filename')
)
def update_feature_ion_mgf_file_path(mgf_contents, mgf_filenames):
    print("MGF upload triggered")
    if mgf_contents and mgf_filenames:
        print(f"MGF file names: {mgf_filenames}")
        mgf_paths = [save_tab5_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
        mgf_path_text = f"Uploaded MGF files: {', '.join(mgf_filenames)}"
    else:
        print("No MGF file uploaded")
        mgf_path_text = "No MGF file uploaded"
    return mgf_path_text

# Callback for updating the path of the uploaded CSV file
@app.callback(
    Output('feature-ion-csv-file-path', 'children'),
    Input('feature-ion-upload-csv', 'contents'),
    State('feature-ion-upload-csv', 'filename')
)
def update_feature_ion_csv_file_path(csv_contents, csv_filenames):
    print("CSV upload triggered")
    if csv_contents and csv_filenames:
        print(f"CSV file names: {csv_filenames}")
        csv_paths = [save_tab5_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

# Callback for adding and removing ion sets dynamically
@app.callback(
    Output('feature-ion-div', 'children'),
    [Input('add-new-ionset', 'n_clicks'),
     Input({'type': 'remove-button-feature-ion', 'index': ALL}, 'n_clicks')],
    [State('feature-ion-div', 'children')]
)
def feature_ion_modify_types(add_clicks, remove_clicks, children):
    print(f"Add new ion set clicked: {add_clicks} times")
    if add_clicks > 0:
        new_index = len(children) - 1
        print(f"Adding new input area at index: {new_index}")
        new_input_area = html.Div(id={'type': 'input-area', 'index': new_index + 1}, children=[
            html.Div([
                html.Label(f'Ion Name {new_index + 1}:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(
                    id={'type': 'input-type-name', 'index': new_index + 1},
                    type='text',
                    placeholder=f'e.g., C10 Ion',
                    style={'width': '10%', 'fontFamily': 'Arial', 'marginLeft': '10px'}
                )
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '20px'}),
            html.H4(f'Specific Ion Set {new_index + 1} :', style={'fontFamily': 'Arial'}),
            dcc.Textarea(
                id={'type': 'input-type-ionset-formulas', 'index': new_index + 1},
                placeholder='Input feature ion formulas',
                style={'width': '72%', 'height': '50px', 'fontFamily': 'Arial'}
            ),
            html.Div([
                html.Label('Mass Range of Product Ions (m/z):',
                           style={'fontFamily': 'Arial', 'fontWeight': 'bold', 'marginRight': '10px'}),
                dcc.Input(id={'type': 'mass-range-min', 'index': new_index + 1}, type='number', placeholder='min',
                          style={'width': '18%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),
                dcc.Input(id={'type': 'dropdown-and-input', 'index': new_index + 1}, type='text',
                          value='Use Precursor Ion', style={'width': '18%', 'fontFamily': 'Arial'}),
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'}),
            html.Div([
                html.Label('Charge:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'charge', 'index': new_index + 1}, type='number', placeholder='e.g., +1/-1',
                          style={'width': '5%', 'marginLeft': '10px', 'marginRight': '20px', 'fontFamily': 'Arial'}),

                html.Label('Mass Tolerance (ppm):', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'tolerance', 'index': new_index + 1}, type='number', placeholder='e.g., 5',
                          style={'width': '5%', 'marginLeft': '10px', 'marginRight': '20px', 'fontFamily': 'Arial'}),

                html.Label('Intensity Normalization:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'normalization', 'index': new_index + 1}, type='number', placeholder='e.g., 0.1',
                          style={'width': '5%', 'marginLeft': '10px', 'marginRight': '20px', 'fontFamily': 'Arial'}),

                html.Label('Top Number:', style={'fontFamily': 'Arial', 'fontWeight': 'bold'}),
                dcc.Input(id={'type': 'top', 'index': new_index + 1}, type='number', placeholder='e.g., 3',
                          style={'width': '5%', 'marginLeft': '10px', 'marginRight': '20px', 'fontFamily': 'Arial'}),

                html.Button('✖', id={'type': 'remove-button-feature-ion', 'index': new_index + 1},
                            style={'marginLeft': '10px', 'color': 'red'})
            ], style={'display': 'flex', 'alignItems': 'center', 'marginTop': '10px'})
        ])
        children.insert(-1, new_input_area)

    # Handling remove button functionality
    ctx = dash.callback_context
    if ctx.triggered:
        triggered = ctx.triggered[0]['prop_id'].split('.')[0]
        if 'remove-button-feature-ion' in triggered:
            button_index = eval(triggered)['index']
            children = [
                child for child in children
                if isinstance(child, dict) and
                'props' in child and
                isinstance(child['props'].get('id', {}), dict) and
                child['props']['id'].get('index') != button_index and
                child['props']['id'] != 'add-new-ionset'
            ]

    return children

# Callback for generating download link for extracted feature ions
@app.callback(
    Output('feature-ion-download-link', 'href'),
    Input('run-feature-ion-extraction', 'n_clicks'),
    State('feature-ion-upload-csv', 'filename'),
    State('feature-ion-upload-csv', 'contents'),
    State('feature-ion-upload-mgf', 'filename'),
    State('feature-ion-upload-mgf', 'contents'),
    State({'type': 'input-type-name', 'index': ALL}, 'value'),
    State({'type': 'input-type-ionset-formulas', 'index': ALL}, 'value'),
    State({'type': 'mass-range-min', 'index': ALL}, 'value'),
    State({'type': 'dropdown-and-input', 'index': ALL}, 'value'),
    State({'type': 'charge', 'index': ALL}, 'value'),
    State({'type': 'tolerance', 'index': ALL}, 'value'),
    State({'type': 'normalization', 'index': ALL}, 'value'),
    State({'type': 'top', 'index': ALL}, 'value')
)
def feature_ion_generate_zip_for_download(n_clicks, csv_filenames, csv_contents, mgf_filenames, mgf_contents,
                                          ion_names, ion_formulas, mass_min, mass_max_option, charges, tolerances,
                                          normalizations, tops):
    if n_clicks is None or n_clicks == 0:
        return ""

    zip_filename = os.path.join(UPLOAD_DIR, 'feature_ion_extraction_results.zip')

    min_len = min(len(ion_names), len(ion_formulas), len(mass_min), len(mass_max_option),
                  len(charges), len(tolerances), len(normalizations), len(tops))

    ion_sets_params = []
    for i in range(min_len):
        ion_set_params = {
            'name': ion_names[i],
            'formulas': ion_formulas[i],
            'mass_min': mass_min[i],
            'mass_max_option': mass_max_option[i],
            'charge': charges[i],
            'tolerance': tolerances[i],
            'intensity_normalization': normalizations[i],
            'top_number': tops[i]
        }
        ion_sets_params.append(ion_set_params)

    # Save MGF and CSV files
    mgf_paths = [save_tab5_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
    csv_paths = [save_tab5_file(name, content) for name, content in zip(csv_filenames, csv_contents)]

    # Process each pair of MGF and CSV files
    for i, (mgf_path, csv_path) in enumerate(zip(mgf_paths, csv_paths)):
        processed_file = csv_filenames[i].replace(".csv", "-fi.csv")
        output_path = os.path.join(UPLOAD_DIR, 'tab_5_uploads', processed_file)

        print(f"Processing MGF: {mgf_path} with CSV: {csv_path} -> Saving to {output_path}")

        # Process files based on the given parameters
        process_FeatureIonExtract_parallel(mgf_path, csv_path, ion_sets_params, output_path)

    # Create a ZIP file of processed files
    with zipfile.ZipFile(zip_filename, 'w') as zf:
        for csv_filename in csv_filenames:
            processed_file = csv_filename.replace(".csv", "-fi.csv")
            file_path = os.path.join(UPLOAD_DIR, 'tab_5_uploads', processed_file)
            if os.path.exists(file_path):
                zf.write(file_path, processed_file)
                print(f"Added {file_path} to ZIP.")
            else:
                print(f"Warning: {file_path} not found.")

    return f'/download/{os.path.basename(zip_filename)}'

# Callback for updating file dropdown after extraction
@app.callback(
    Output('feature-ion-extract-file-dropdown', 'options'),
    Input('run-feature-ion-extraction', 'n_clicks'),
    State('feature-ion-upload-csv', 'filename')
)
def featureion_update_file_dropdown(n_clicks, csv_filenames):
    if n_clicks is None or n_clicks == 0:
        return []

    return [{'label': filename.replace(".csv", "-fi.csv"), 'value': filename.replace(".csv", "-fi.csv")}
            for filename in csv_filenames]

# Callback for displaying the results of the extraction in a table
@app.callback(
    Output('feature-ion-extract-results-table', 'children'),
    Input('feature-ion-extract-file-dropdown', 'value')
)
def feature_ion_display_selected_file_results(selected_file):
    if not selected_file:
        return "Please select a file to view results."

    file_path = os.path.join(UPLOAD_DIR, 'tab_5_uploads', selected_file)
    df = pd.read_csv(file_path)

    return dash_table.DataTable(
        data=df.to_dict('records'),
        columns=[{'name': i, 'id': i} for i in df.columns],
        style_header={'fontFamily': 'Arial', 'fontWeight': 'bold'},
        style_cell={'fontFamily': 'Arial', 'textAlign': 'left'}
    )

# Route for downloading the generated ZIP file
@app.server.route('/download/<path:filename>')
def feature_ion_download(filename):
    return send_file(os.path.join(UPLOAD_DIR, filename), as_attachment=True)

def safe_get_field(mr_dict, key, field):
    """
    Secure the fields safely from matched_results.
    If the mr_dictkeyis dictates, it returns the mr_dictkeyfield,
 otherwise it returns none.
    """
    val = mr_dict.get(key)
    if isinstance(val, dict):
        return val.get(field, 'none')
    return 'none'

def detect_formulas_in_spectra(spectrum, specific_formulas, params):
    """
    Detect specific formulas in a given spectrum based on mass and intensity parameters.

    Parameters:
        spectrum (dict): The mass spectrum data containing 'params', 'm/z array', and 'intensity array'.
        specific_formulas (list): List of chemical formulas to match against the spectrum.
        params (dict): A dictionary containing parameters for detection including mass range, charge, tolerance, etc.

    Returns:
        matched_formulas (dict): A dictionary where keys are spectrum titles and values are matched formulas.
    """
    matched_formulas = {}
    try:
        # Extract parameters for formula detection
        min_mass = params['mass_min']
        charge = params['charge']
        tolerance_ppm = params['tolerance']
        norScore = params['intensity_normalization']
        top_number = params['top_number']
        max_mass_option = params.get('max_mass_option')

        # Obtain the spectrum's title and mass range
        title = spectrum['params']['title']
        if max_mass_option == "Specify Value":
            max_mass = params['mass_max']
        else:
            max_mass = spectrum['params']['pepmass'][0]

        # Filter product ions within the specified mass range and above the intensity threshold
        product_ions = spectrum['m/z array']
        intensity = spectrum['intensity array']
        ind = np.nonzero((product_ions > min_mass) & (product_ions < max_mass) & (intensity > norScore))[0]

        # Check if any ions meet the conditions
        if len(ind) > 0:
            intensity_max = np.max(intensity[ind])
            intensity_min = np.min(intensity[ind])
            normalized_intensity = np.zeros_like(intensity)
            normalized_intensity[ind] = (intensity[ind] - intensity_min) / (intensity_max - intensity_min)
            formula_intensity_pairs = []

            # Compare each formula's m/z range against the product ions in the spectrum
            for formula in specific_formulas:
                mz = mass.calculate_mass(formula=formula, charge=charge)
                mz_min = mz - mz * tolerance_ppm / 1e6
                mz_max = mz + mz * tolerance_ppm / 1e6

                # Identify matches within the specified m/z range and add to result
                for i in ind:
                    if mz_min <= product_ions[i] <= mz_max and normalized_intensity[i] > norScore:
                        formula_intensity_pairs.append((formula, normalized_intensity[i]))

            # Sort and select the top formulas based on intensity
            formula_intensity_pairs.sort(key=lambda x: x[1], reverse=True)
            top_formulas = [pair[0] for pair in formula_intensity_pairs[:top_number]]
            if top_formulas:
                matched_formulas[int(title)] = ', '.join(top_formulas)

    except Exception as e:
        print("Error during formula detection:", str(e))

    return matched_formulas

def process_FeatureIonExtract_parallel(mgf_path, csv_path, ion_sets_params, full_output_path):
    """
    For each peak in the CSV/MGF pair, extract feature-ion matches per ion-set,
    keep the top N by absolute intensity, then expand into one row per combination
    of those Top-N lists across all ion-sets.
    """
    # 1) Read the input CSV file and set 'row ID' as the DataFrame index if present
    df = pd.read_csv(csv_path, encoding='utf-8')
    if 'row ID' in df.columns:
        df.set_index('row ID', inplace=True)
    try:
        df.index = df.index.astype(int)
    except:
        # If conversion to int fails, leave the index as-is
        pass

    # 2) Load all spectra from the MGF file into a dictionary keyed by title
    spectra_dict = {}
    with mgf.read(mgf_path) as spectra:
        for spectrum in spectra:
            title = int(spectrum['params']['title'])
            spectra_dict[title] = spectrum

    # 3) Prepare a dictionary of formulas for each ion-set
    specific_formulas = {
        params['name']: [f.strip() for f in params['formulas'].split(',')]
        for params in ion_sets_params
    }

    # 4) Make a deep copy of the original DataFrame for row expansion
    original_df = df.copy(deep=True)
    output_records = []

    # 5) Process each peak (row) in the original DataFrame
    for row_id, orig_row in original_df.iterrows():
        if row_id not in spectra_dict:
            continue  # Skip if no matching spectrum

        spectrum = spectra_dict[row_id]
        pepmz = spectrum['params']['pepmass'][0]
        mzs   = spectrum['m/z array']
        ints  = spectrum['intensity array']

        per_group_topN = []  # Will hold Top-N lists for each ion-set

        # 5a) For each ion set, find all formula matches, sort, and keep Top-N
        for params in ion_sets_params:
            name       = params['name']
            formulas   = specific_formulas[name]
            min_mass   = params['mass_min']
            tol_ppm    = params['tolerance']
            charge     = params['charge']
            norm_thr   = params['intensity_normalization']
            top_n      = params['top_number'] or 3
            max_opt    = params.get('mass_max_option', 'Use Precursor Ion')

            # Determine upper mass limit: either user-specified or the precursor m/z
            if max_opt != "Use Precursor Ion":
                try:
                    max_mass = float(max_opt)
                except:
                    max_mass = pepmz
            else:
                max_mass = pepmz

            # Filter indices by m/z range and intensity threshold
            idxs = [
                i for i,(mz_val,inten) in enumerate(zip(mzs,ints))
                if mz_val>min_mass and mz_val<max_mass and inten>norm_thr
            ]

            # Collect all matching formula-intensity pairs
            candidates = []
            for formula in formulas:
                target_mz = mass.calculate_mass(formula=formula, charge=charge)
                lo = target_mz*(1 - tol_ppm/1e6)
                hi = target_mz*(1 + tol_ppm/1e6)
                for i in idxs:
                    if lo <= mzs[i] <= hi:
                        candidates.append((formula, ints[i]))

            # Sort candidates by absolute intensity, descending, then take Top-N
            candidates.sort(key=lambda x: x[1], reverse=True)
            top_list = candidates[:top_n] if candidates else [(None,None)]

            # Attach a normalized score (intensity/100000) to each top hit
            scored = []
            for formula,inten in top_list:
                if formula is None:
                    scored.append((None, None, None))
                else:
                    score = round(inten/100000, 4)
                    scored.append((formula, inten, score))
            per_group_topN.append(scored)

        # 6) Expand into one row per combination of Top-N hits across ion-sets
        for combo in itertools.product(*per_group_topN):
            rec = orig_row.to_dict()                # Copy original row data
            rec['row ID']  = row_id                # Preserve original row ID
            rec['row m/z'] = pepmz                 # Add precursor m/z value
            # Assign each ion-set's selected formula/intensity/score
            for params,(formula,inten,score) in zip(ion_sets_params,combo):
                name = params['name']
                if formula is None:
                    rec[name]             = 'none'
                    rec[f"{name}_intensity"] = 'none'
                    rec[f"{name}_score"]     = 'none'
                else:
                    rec[name]             = formula
                    rec[f"{name}_intensity"] = f"{inten:.2f}"
                    rec[f"{name}_score"]     = f"{score:.4f}"
            output_records.append(rec)

    # 7) Build a DataFrame from the exploded records and write to CSV
    out_df = pd.DataFrame(output_records)
    # Ensure 'row ID' is the first column
    cols = ['row ID'] + [c for c in out_df.columns if c!='row ID']
    out_df.reset_index(inplace=True)
    out_df = out_df[cols]
    # Create output directory if needed, then write CSV without index
    os.makedirs(os.path.dirname(full_output_path), exist_ok=True)
    out_df.to_csv(full_output_path, index=False)
    print(f"Saved processed file to: {full_output_path}")

def save_tab5_file(name, content, tab_folder='tab_5_uploads'):
    """
    Save uploaded MGF or CSV content to a specific folder for further processing.

    Parameters:
        name (str): Filename of the uploaded content.
        content (str): Base64-encoded content of the file.
        tab_folder (str): Folder path to save the file.

    Returns:
        file_path (str): Full path where the file is saved.
    """
    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)
    print(f"Saving file {name} to {file_path}")
    with open(file_path, 'wb') as f:
        f.write(decoded)
    return file_path


""" Tab 6 related Functions """
@app.callback(
    Output('cosine-download-link', 'href'),
    Input('run-cosine-score-calculation', 'n_clicks'),
    State('cosine-query-upload-csv', 'filename'),
    State('cosine-query-upload-csv', 'contents'),
    State('cosine-query-upload-mgf', 'filename'),
    State('cosine-query-upload-mgf', 'contents'),
    State('cosine-standard-upload-csv', 'filename'),
    State('cosine-standard-upload-csv', 'contents'),
    State('cosine-standard-upload-mgf', 'filename'),
    State('cosine-standard-upload-mgf', 'contents'),
    State('cosine-input-intensity-normalization', 'value'),
    State('cosine-input-cosine-score', 'value'),
    State({'type': 'cosine-mass-range-min'}, 'value'),
    State({'type': 'cosine-dropdown-and-input'}, 'value')
)
def cosine_generate_zip_for_download(n_clicks, query_csv_filenames, query_csv_contents, query_mgf_filenames, query_mgf_contents,
                                     std_csv_filenames, std_csv_contents, std_mgf_filenames, std_mgf_contents,
                                     intensity_normalization, cosine_score, min_mass, max_mass_option):
    if n_clicks is None or n_clicks == 0:
        return ""

    if min_mass is None or intensity_normalization is None or cosine_score is None:
        return ""

    zip_filename = os.path.join(UPLOAD_DIR, 'cosine_results.zip')

    mgf_paths_query = [save_tab6_file(name, content) for name, content in zip(query_mgf_filenames, query_mgf_contents)]
    csv_paths_query = [save_tab6_file(name, content) for name, content in zip(query_csv_filenames, query_csv_contents)]
    mgf_paths_std = [save_tab6_file(name, content) for name, content in zip(std_mgf_filenames, std_mgf_contents)]
    csv_paths_std = [save_tab6_file(name, content) for name, content in zip(std_csv_filenames, std_csv_contents)]

    for query_mgf, query_csv, std_mgf, std_csv in zip(mgf_paths_query, csv_paths_query, mgf_paths_std, csv_paths_std):
        processed_file = query_csv_filenames[0].replace(".csv", "-cs.csv")
        output_path = os.path.join(UPLOAD_DIR, 'tab_6_uploads', processed_file)

        process_cosine_score_calculate(query_mgf, std_mgf, query_csv, std_csv, output_path, min_mass, intensity_normalization, max_mass_option, cosine_score)

    with zipfile.ZipFile(zip_filename, 'w') as zf:
        for csv_filename in query_csv_filenames:
            processed_file = csv_filename.replace(".csv", "-cs.csv")
            file_path = os.path.join(UPLOAD_DIR, 'tab_6_uploads', processed_file)
            if os.path.exists(file_path):
                zf.write(file_path, processed_file)
            else:
                print(f"Warning: {file_path} not found.")

    return f'/download/{os.path.basename(zip_filename)}'

@app.callback(
    Output({'type': 'cosine-dropdown-and-input', 'index': MATCH}, 'value'),
    [Input({'type': 'cosine-dropdown-and-input', 'index': MATCH}, 'n_submit')],
    [Input({'type': 'cosine-dropdown-and-input', 'index': MATCH}, 'value')]
)
def update_value(n_submit, value):
    if value == "Use Precursor Ion":
        return value
    try:
        float_value = float(value)
        return float_value
    except ValueError:
        return value

@app.callback(
    Output('cosine-query-mgf-file-path', 'children'),
    Input('cosine-query-upload-mgf', 'contents'),
    State('cosine-query-upload-mgf', 'filename')
)
def update_cosine_query_mgf_file_path(mgf_contents, mgf_filenames):
    print("MGF upload triggered")
    if mgf_contents and mgf_filenames:
        print(f"MGF file names: {mgf_filenames}")
        mgf_paths = [save_tab6_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
        mgf_path_text = f"Uploaded MGF files: {', '.join(mgf_filenames)}"
    else:
        print("No MGF file uploaded")
        mgf_path_text = "No MGF file uploaded"
    return mgf_path_text

@app.callback(
    Output('cosine-query-csv-file-path', 'children'),
    Input('cosine-query-upload-csv', 'contents'),
    State('cosine-query-upload-csv', 'filename')
)
def update_cosine_query_csv_file_path(csv_contents, csv_filenames):
    print("CSV upload triggered")
    if csv_contents and csv_filenames:
        print(f"CSV file names: {csv_filenames}")
        csv_paths = [save_tab6_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

@app.callback(
    Output('cosine-standard-mgf-file-path', 'children'),
    Input('cosine-standard-upload-mgf', 'contents'),
    State('cosine-standard-upload-mgf', 'filename')
)
def update_cosine_standard_mgf_file_path(mgf_contents, mgf_filenames):
    print("MGF upload triggered")
    if mgf_contents and mgf_filenames:
        print(f"MGF file names: {mgf_filenames}")
        mgf_paths = [save_tab6_file(name, content) for name, content in zip(mgf_filenames, mgf_contents)]
        mgf_path_text = f"Uploaded MGF files: {', '.join(mgf_filenames)}"
    else:
        print("No MGF file uploaded")
        mgf_path_text = "No MGF file uploaded"
    return mgf_path_text

@app.callback(
    Output('cosine-standard-csv-file-path', 'children'),
    Input('cosine-standard-upload-csv', 'contents'),
    State('cosine-standard-upload-csv', 'filename')
)
def update_cosine_standard_csv_file_path(csv_contents, csv_filenames):
    print("CSV upload triggered")
    if csv_contents and csv_filenames:
        print(f"CSV file names: {csv_filenames}")
        csv_paths = [save_tab6_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

@app.callback(
    Output('cosine-file-dropdown', 'options'),
    Input('run-cosine-score-calculation', 'n_clicks'),
    State('cosine-query-upload-csv', 'filename')
)
def cosine_update_file_dropdown(n_clicks, csv_filenames):
    if n_clicks is None or n_clicks == 0:
        return []

    return [{'label': filename.replace(".csv", "-cs.csv"), 'value': filename.replace(".csv", "-cs.csv")}
            for filename in csv_filenames]

@app.callback(
    Output('cosine-results-table', 'children'),
    Input('cosine-file-dropdown', 'value')
)
def cosine_display_selected_file_results(selected_file):
    if not selected_file:
        return "Please select a file to view results."

    file_path = os.path.join(UPLOAD_DIR, 'tab_6_uploads', selected_file)
    df = pd.read_csv(file_path)

    return dash_table.DataTable(
        data=df.to_dict('records'),
        columns=[{'name': i, 'id': i} for i in df.columns],
        style_header={'fontFamily': 'Arial', 'fontWeight': 'bold'},
        style_cell={'fontFamily': 'Arial', 'textAlign': 'left'}
    )

@app.server.route('/download/<path:filename>')
def cosine_download(filename):
    return send_file(os.path.join(UPLOAD_DIR, filename), as_attachment=True)

def filter_spectra_with_normalization(filename, norScore, mass_ranges_min, mass_ranges_max_options):
    if mass_ranges_min is None or norScore is None:
        raise ValueError("min_mass and norScore must be provided")

    filtered_spectra = []
    with mgf.read(filename) as spectra:
        for spectrum in spectra:
            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']
            max_mass_option = mass_ranges_max_options

            if max_mass_option == "Use Precursor Ion":
                max_mass_value = spectrum['params']['pepmass'][0]
            else:
                try:
                    max_mass_value = float(max_mass_option)
                except ValueError:
                    max_mass_value = None

            ind = (mz_array >= mass_ranges_min) & (mz_array <= max_mass_value) & (intensity_array > norScore)
            if np.any(ind):
                mz_array = mz_array[ind]
                intensity_array = intensity_array[ind]
                intensity_max = np.max(intensity_array)
                intensity_min = np.min(intensity_array)
                normalized_intensities = (intensity_array - intensity_min) / (intensity_max - intensity_min)
                filtered_spectra.append({
                    'm/z array': mz_array,
                    'intensity array': normalized_intensities,
                    'params': spectrum['params']
                })
    return filtered_spectra

def convert_to_matchms_spectra(filtered_spectra):
    matchms_spectra = []
    for spec in filtered_spectra:
        spectrum = Spectrum(mz=spec['m/z array'],
                            intensities=spec['intensity array'],
                            metadata=spec['params'])
        matchms_spectra.append(spectrum)
    return matchms_spectra

def process_cosine_score_calculate(query_path_mgf, std_path_mgf, query_path_csv, std_path_csv, full_output_path, mass_ranges_min, norScore, mass_ranges_max_options, CosineScore):
    # Step 1: Normalize and filter the spectra
    query_spectra_filtered = convert_to_matchms_spectra(filter_spectra_with_normalization(query_path_mgf, mass_ranges_min, norScore, mass_ranges_max_options))
    std_spectra_filtered = convert_to_matchms_spectra(filter_spectra_with_normalization(std_path_mgf, mass_ranges_min, norScore, mass_ranges_max_options))

    # Step 2: Calculate similarity scores using matchms
    similarity_scores = calculate_scores(query_spectra_filtered, std_spectra_filtered, CosineGreedy(tolerance=0.005), is_symmetric=False)

    # Step 3: Read query and standard CSV files
    query_csv = pd.read_csv(query_path_csv)
    std_csv = pd.read_csv(std_path_csv)
    query_csv['row ID'] = query_csv['row ID'].astype(str)
    std_csv['row ID'] = std_csv['row ID'].astype(str)

    best_matches = {}
    cosine_scores = {}

    # Step 4: Find the best matches based on cosine score
    for score in similarity_scores:
        query_spectrum, std_spectrum, score_values = score
        cosine_score, matched_peaks = score_values
        if cosine_score > CosineScore:
            query_title = query_spectrum.metadata['title']
            std_title = std_spectrum.metadata['title']
            if std_title in std_csv['row ID'].values:
                substructure = std_csv.loc[std_csv['row ID'] == std_title, 'substructure type'].iloc[0]
                best_matches[query_title] = substructure
                cosine_scores[query_title] = cosine_score

    # Step 5: Update query CSV with the results, appending columns instead of inserting
    query_csv['cosine score'] = query_csv['row ID'].map(cosine_scores).round(2).fillna('value < 0.7')
    query_csv['substructure type'] = query_csv['row ID'].map(best_matches).fillna('others')

    # Just append the new columns at the end
    query_csv.to_csv(full_output_path, index=False)

def save_tab6_file(name, content, tab_folder='tab_6_uploads'):
    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)
    with open(file_path, 'wb') as f:
        f.write(decoded)
    return file_path


""" Tab 7 related Functions """
@app.callback(
    Output('annotation-query-csv-file-path', 'children'),
    Input('annotation-query-upload-csv', 'contents'),
    State('annotation-query-upload-csv', 'filename')
)
def update_annotation_query_csv_file_path(csv_contents, csv_filenames):
    print("CSV upload triggered")
    # Check if both content and filenames are provided
    if csv_contents and csv_filenames:
        print(f"CSV file names: {csv_filenames}")
        # Save each file and store the paths
        csv_paths = [save_tab7_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        print(f"CSV paths saved: {csv_paths}")
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

@app.callback(
    Output('annotation-pseudo-DB-file-path', 'children'),
    Input('annotation-pseudo-DB-upload-db', 'contents'),
    State('annotation-pseudo-DB-upload-db', 'filename')
)
def update_annotation_pseudo_DB_file_path(content, filename):
    print("pseudo-DB upload triggered")
    if content and filename:
        print(f"pseudo-DB file name: {filename}")
        saved_file = save_uploaded_file(filename, content)
        if saved_file:
            print(f"File saved at: {saved_file}")
            return f"Uploaded pseudo-DB file: {filename}"
        else:
            print("Failed to upload the file.")
            return "Failed to upload the file."
    print("No DB file uploaded")
    return "No DB file uploaded"

@app.callback(
    Output('annotation-substituent-csv-file-path', 'children'),
    Input('annotation-substituent-upload-csv', 'contents'),
    State('annotation-substituent-upload-csv', 'filename')
)
def update_annotation_substituent_csv_file_path(csv_contents, csv_filenames):
    print("Substituent CSV upload triggered")
    if csv_contents and csv_filenames:
        print(f"Substituent CSV file names: {csv_filenames}")
        csv_paths = [save_tab7_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        print(f"CSV paths saved: {csv_paths}")
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

@app.callback(
    Output('annotation-file-dropdown', 'options'),
    Input('run-annotation', 'n_clicks'),
    State('annotation-query-upload-csv', 'filename')
)
def annotation_file_dropdown(n_clicks, csv_filenames):
    print(f"Annotation run clicked: {n_clicks} times")
    if n_clicks is None or n_clicks == 0:
        return []
    print(f"Files available for dropdown: {csv_filenames}")
    return [{'label': filename.replace(".csv", "-an.csv"), 'value': filename.replace(".csv", "-an.csv")}
            for filename in csv_filenames]

@app.callback(
    Output('annotation-results-table', 'children'),
    Input('annotation-file-dropdown', 'value')
)
def annotation_display_selected_file_results(selected_file):
    if not selected_file:
        print("No file selected for results display.")
        return "Please select a file to view results."

    file_path = os.path.join(UPLOAD_DIR, 'tab_7_uploads', selected_file)
    print(f"Trying to open file: {file_path}")

    try:
        df = pd.read_csv(file_path)
        print(f"File {selected_file} loaded successfully.")
        return dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns],
            style_header={'fontFamily': 'Arial', 'fontWeight': 'bold'},
            style_cell={'fontFamily': 'Arial', 'textAlign': 'left'}
        )
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return f"Error: {str(e)}"

@app.server.route('/download/<path:filename>')
def annotation_download(filename):
    file_path = os.path.join(UPLOAD_DIR, filename)
    print(f"Initiating download for file: {file_path}")
    return send_file(file_path, as_attachment=True)

@app.callback(
    Output('annotation-output', 'children'),
    Input('run-annotation', 'n_clicks'),
    State('annotation-query-upload-csv', 'filename'),
    State('annotation-pseudo-DB-upload-db', 'filename'),
    State('annotation-substituent-upload-csv', 'filename'),
)
def run_annotation(n_clicks, query_fns, pseudo_db_fn, acyl_fns):
    if not n_clicks:
        return "Click Run Annotation to start processing."
    query_file = query_fns[0]          if isinstance(query_fns, list) else query_fns
    acyl_file  = acyl_fns[0]            if isinstance(acyl_fns, list)  else acyl_fns
    query_path = os.path.join(UPLOAD_DIR, 'tab_7_uploads', query_file)
    acyl_path  = os.path.join(UPLOAD_DIR, 'tab_7_uploads', acyl_file)
    pseudo_path= os.path.join(UPLOAD_DIR, 'tab_7_uploads', pseudo_db_fn)

    base, ext = os.path.splitext(query_file)
    out_fname  = f"{base}-an{ext}"
    out_path   = os.path.join(UPLOAD_DIR, 'tab_7_uploads', out_fname)

    annotation_process_data(
        query_path,
        acyl_path,
        pseudo_path,
        out_path
    )
    return f"Annotation completed. Output saved to {out_path}"

def save_uploaded_file(name, content, folder='tab_7_uploads'):
    folder_path = os.path.join(UPLOAD_DIR, folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    print(f"Saving file {name} in folder {folder}")

    try:
        content_type, content_string = content.split(',')
        print(f"Content type: {content_type}")
    except ValueError as e:
        print(f"Error splitting content for file {name}: {e}")
        return None

    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)

    with open(file_path, 'wb') as f:
        f.write(decoded)
    print(f"File saved successfully at {file_path}")
    return file_path

def save_tab7_file(name, content, tab_folder='tab_7_uploads'):
    print(f"Saving file: {name} to folder: {tab_folder}")
    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)
    with open(file_path, 'wb') as f:
        f.write(decoded)
    print(f"File saved at: {file_path}")
    return file_path

# --- Utility function: Count the number of a specific element in a formula ---
def count_element(formula: str, element: str) -> int:
    pattern = re.compile(rf"{element}(\d*)")
    match = pattern.search(formula)
    if match:
        return int(match.group(1) or 1)
    return 0

# --- Compute the Index of Hydrogen Deficiency (IHD) ---
def determine_ihd(row) -> float:
    comp = row.get("Elemental Composition","")
    C = count_element(comp,"C")
    H = count_element(comp,"H")
    N = count_element(comp,"N")
    return (2*C + 2 + N - H)/2

# --- Normalize formula (C, H first, then others alphabetically) ---
def normalize_formula(formula: str) -> str:
    comp = mass.Composition(formula or "")
    parts = []
    for e in ("C","H"):
        if comp.get(e):
            parts.append(f"{e}{comp[e]}")
    for e in sorted(set(comp)-{"C","H"}):
        parts.append(f"{e}{comp[e]}")
    return "".join(parts)

# --- Solver: Find (x, y) such that 2*x + 7*y = remain (with limits) ---
def solve_2_7(remain: int) -> Optional[Tuple[int, int]]:
    max_x, max_y = 3, 3
    for y in range(0, min(remain // 7, max_y) + 1):
        rest = remain - 7 * y
        if rest >= 0 and rest % 2 == 0:
            x = rest // 2
            if x <= max_x:
                return x, y
    return None

# --- Solver: 2*x + 7*y + 4*z + 5*m = remain, minimize total parts, all vars ≤ 5 ---
def solve_2_9_10_7_11(remain:int)->Optional[Tuple[int,int,int,int,int]]:
    best=None
    for a in range(6):
        for b in range(3):
            for c in range(3):
                for d in range(4):
                    for e in range(3):
                        if 2*a+9*b+10*c+7*d+11*e==remain:
                            tot=a+b+c+d+e
                            if best is None or tot<best[0]:
                                best=(tot,a,b,c,d,e)
    return best[1:] if best else None

def solve_2_7_4_5(remain:int)->List[Tuple[int,int,int,int]]:
    out=[]
    for x in range(6):
        for y in range(6):
            for z in range(6):
                if 2*x+7*y+4*z==remain:
                    out.append((x,y,z,0))
    for x in range(6):
        for m in range(6):
            if 2*x+5*m==remain:
                out.append((x,0,0,m))
    return out

# --- Determine number and type of substituent fragments ---
def determine_substitution(row):
    C_total=count_element(row.get("elemental composition",""),"C")
    t=row.get("type","")
    rows=[]
    if t=="D":
        remain=C_total-20
        sol=solve_2_9_10_7_11(remain)
        if sol:
            names=["C2","C9","C10","C7","C11"]
            terms=[f"{n}*{c}" for n,c in zip(names,sol) if c>0]
            subs=sum(([n]*c for n,c in zip(names,sol)),[])
            new=row.copy()
            new["Sub No."]=len(subs)
            new["Sub Carbon"]=", ".join(sorted(set(subs),key=lambda x:int(x[1:])))
            new["Sub_An"]="+".join(terms)
            rows.append(new)
    elif t=="MD":
        core_map={"C9":29,"C10":30,"C14":34}
        for frag,core in core_map.items():
            if frag in row.get("M part",""):
                remain=C_total-core
                for sol in solve_2_7_4_5(remain):
                    names=["C2","C7","C4","C5"]
                    terms=[f"{n}*{c}" for n,c in zip(names,sol) if c>0]
                    subs=sum(([n]*c for n,c in zip(names,sol)),[])
                    new=row.copy()
                    new["Sub No."]=len(subs)
                    new["Sub Carbon"]=", ".join(sorted(set(subs),key=lambda x:int(x[1:])))
                    new["Sub_An"]="+".join(terms)
                    rows.append(new)
                break
    if rows: return rows
    row["Sub No."]=0; row["Sub Carbon"]=""; row["Sub_An"]=""
    return [row]

# ------ Part A ------
def determine_a_parts(row):
    if row['type']=="D":
        o = count_element(row['C17 ion'],"O")
        h = count_element(row['C17 ion'],"H")
        o18 = count_element(row.get('C18 ion',""),"O")
        if o in (2,3):
            if h==16: return "1,2-en-3-one"
            if h==20: return "1,2-dihydro-3-ol"
            if h==18:
                if o18==1: return "1,2-en-3-ol"
                if o18==2: return "1,2-dihydro-3-one"
    elif row['type']=="MD":
        try: c16 = float(row.get("C16 ion_score",0))
        except: c16=0
        c15 = str(row.get("C15 ion","")).strip().lower()
        if c15=="none":
            if row.get("C16 ion")== "C16H16O3" and row.get("C20 ion")== "C20H16O2":
                return "1-alkyl-2-ol-3-one"
            if row.get("C16 ion")== "C16H18O3" and row.get("C20 ion")== "C20H18O2":
                return "1-alkyl-2-ol-3-ol"
        if row.get("C15 ion")== "C15H22O2" and row.get("C16 ion")== "C16H16O3":
            return "bicyclo [2.2.1] heptane ring"
        if row.get("C14 ion")== "C14H14O3" and row.get("C19 ion") in ("C19H18O4","C19H20O4"):
            return "bicyclo [2.2.1] heptane ring"
        if row.get("C14 loss")== "C14H24O2":
            if row.get("C14 ion")== "C14H14O3" and row.get("C19 ion") in ("C19H18O4","C19H20O4"):
                return "bicyclo [2.2.1] heptane ring"
            elif row.get("C19 ion")== "C19H18O4":
                return "bicyclo [2.2.1] heptane ring"
        if row.get("C14 ion")== "C14H14O3" and row.get("C15 ion") in ("C15H22O2","C15H22O3") and row.get("C19 ion")== "C19H18O3":
            return "bicyclo [2.2.1] heptane ring"
        if row.get("C15 ion")== "C15H22O3" and row.get("C14 loss")== "C14H24O4" and (row.get("C16 ion","")=="none" or c16<1):
            return "3,4-seco"
        if row.get("C20 ion") in ("C20H22O3","C20H20O3","C20H18O2") and row["IHD"]>11:
            return "1-alkyl-3-ol"
        return "1-alkyl-3-one"
    return "others"

# ------ Part B ------
def determine_b_parts(row, df):
    if row['type'] == "D":
        try:
            s = float(row.get("C19 ion_score", 0))
        except:
            s = 0
        if row.get("C19 ion") == "C19H20O1" and s >= 1:
            try:
                c3 = float(row.get("C3 loss_score", 0))
            except:
                c3 = 0
            if c3 <= 100:
                return "5-ol-6,7-diol"

        if row.get("C3 loss") == "C3H4O2" and count_element(row.get("C17 ion", ""), "O") in (2, 3):
            return "5-ol-6,7-epoxy"

        if row.get("C15 ion") == "C15H18O2":
            if row.get("adduct ion") == "[M+NH4]+":
                rt, ec = row["row retention time"], row["elemental composition"]
                cand = df[
                    (abs(df["row retention time"] - rt) <= 0.1) &
                    (df["elemental composition"] == ec) &
                    (df["adduct ion"] == "[M+H]+")
                    ]
                if not cand.empty:
                    return "5-ol-4,6-epoxy"
            return "5-ol-4,7-epoxy"

        return "5-ol-6,7-ene"

    elif row['type'] == "MD":
        if any(row.get(c, "none") != "none" for c in ("C3 loss", "C17 ion", "C27 ion")):
            return "5-ol-6,7-epoxy"

    return "others"

# ------ Part C ------
def determine_c_parts(df: pd.DataFrame) -> pd.Series:
    triol_ids=set()
    for _,r in df.iterrows():
        if r["type"]=="D" and r["adduct ion"]=="[M+H-H2O]+":
            sp=r.get("same peak","")
            if isinstance(sp,str):
                for i in sp.split(","):
                    if i.isdigit(): triol_ids.add(int(i))
    def f(row):
        if int(row["row ID"]) in triol_ids: return "9,13,14-triol"
        return "9,13,14-orthoester"
    return df.apply(f,axis=1)

# ------ Part C-12 ------
def determine_c12(row):
    if row["type"]=="D":
        if count_element(row["C17 ion"],"O")==2: return "no C-12"
        return "with C-12"
    elif row["type"]=="MD":
        if count_element(row.get("C14 loss",""),"O")==2:
            return "with C-12"
    return "no C-12"

# ------ Part C-18 ------
def determine_c18(row,df):
    if row["type"]!="MD": return "no C-18"
    ad,rt,ec=row["adduct ion"],row["row retention time"],row["elemental composition"]
    if ad=="[M+H]+": targets=("[M+H]+","[M+NH4]+")
    elif ad=="[M+NH4]+": targets=("[M+H]+",)
    else: return "no C-18"
    cands=df[(abs(df["row retention time"]-rt)<=0.05)&(df["adduct ion"].isin(targets))&(df.index!=row.name)]
    allowed = {
        "[M+H]+":[{"C":7,"H":6,"O":2,"N":0},{"C":2,"H":4,"O":2,"N":0},{"C":5,"H":10,"O":2,"N":0},{"C":4,"H":8,"O":2,"N":0}],
        "[M+NH4]+":[{"C":7,"H":10,"O":2,"N":1},{"C":2,"H":8,"O":2,"N":1},{"C":5,"H":14,"O":2,"N":1},{"C":4,"H":12,"O":2,"N":1}]
    }
    for _,m in cands.iterrows():
        diff=calculate_formula_difference(ec,m["elemental composition"])
        if diff in allowed.get(ad,[]): return "with C-18"
    return "no C-18"

# ------ Olefinic ------
def determine_olefinic_part(row):
    if row.get("C5 ion")== "C5H4O" or row.get("C6 ion")== "C6H6O":
        return "2,4-olefinic"
    return "others"

# ------ Part M ------
def determine_m_parts(row):
    if row["type"]!="MD": return ["others"]
    parts=[]
    ion=row.get("C10 ion","")
    try: score=float(row.get("C10 ion_score",0))
    except: score=0
    if ion=="C10H16O3": parts.append("C10 ring with one substituents and a hydroxyl group")
    elif ion=="C10H14O3": parts.append("C10 ring with two substituents and a hydroxyl group")
    elif ion=="C10H12O3" and score>0.015: parts.append("C10 ring with three substituents and a hydroxyl group")
    # loss-based
    loss_map={
        "C9":("C9 loss","C9 loss_score"),
        "C10":("C10 loss","C10 loss_score"),
        "C14":("C14 loss","C14 loss_score")
    }
    for k,(col,_) in loss_map.items():
        lf=row.get(col,"")
        if lf!="none":
            h,o=count_element(lf,"H"),count_element(lf,"O")
            if k=="C10":
                if h==18 and o==2: parts.append("C10 ring with no substituents")
                elif h==16 and o==2: parts.append("C10 ring with one substituent")
                elif h==16 and o==3: parts.append("C10 ring with one substituents and a hydroxyl group")
                elif h==14 and o==2: parts.append("C10 ring with two substituents")
            elif k=="C9" and h==16: parts.append("C9 ring with no substituents")
            elif k=="C14" and o==2: parts.append("C14 ring with an olenfinic")
    return parts or ["others"]

# ------ Acyl ------
def determine_acyl_parts(row,acyl_df):
    # loss-based
    parts=[]
    loss_map={
        "C2 loss":{"C2H4O2":"acyl"},
        "C7 loss":{"C7H6O2":"benzoyl"},
        "C9 loss":{"C9H8O2":"cinnamoyl","C9H8O3":"coumaroyl","C9H18O2":"nonanoyl"},
        "C10 loss":{"C10H10O4":"feruloyl","C10H20O2":"decanoyl"},
        "C11 loss":{"C11H22O2":"undecanoyl"}
    }
    for col,d in loss_map.items():
        f=row.get(col,"")
        if f in d:
            try: sc=float(row.get(f"{col}_score",0))
            except: sc=None
            parts.append(f"{d[f]}({sc:.2f})" if sc is not None else d[f])
    # ion-based
    norm_map=acyl_df.groupby(acyl_df["Acyl Formula"].map(normalize_formula))["Identification"].apply(list).to_dict()
    for ion_col in ("Acyl A ion","Acyl B ion"):
        ions=str(row.get(ion_col,"")).split(",")
        scores=str(row.get(f"{ion_col}_score","")).split(",")
        for ion,sc in zip(ions+[""]*10,scores+[""]*10):
            ion=ion.strip(); sc=sc.strip()
            if not ion or ion.lower()=="none": break
            for ident in norm_map.get(normalize_formula(ion),[]):
                try: s=float(sc)
                except: s=None
                parts.append(f"{ident}({s:.2f})" if s is not None else ident)
    return ",".join(dict.fromkeys(parts)) or "others"

def refine_acyl_parts_by_olefinic(row,acyl_df):
    if row.get("olefinic position")!="2,4-olefinic": return row.get("Acyl part","")
    ents=[e.strip() for e in str(row.get("Acyl part","")).split(",") if e.strip()]
    name2f={}
    for e in ents:
        nm=re.match(r"(.+?)(?:\(|$)",e).group(1)
        sel=acyl_df[acyl_df["Identification"]==nm]
        if not sel.empty: name2f[nm]=sel.iloc[0]["Acyl Formula"]
    rm=set()
    for grp in defaultdict(list, **{f:[] for f in set(name2f.values())}).values():
        ien=[n for n in name2f if "ienoyl" in n]
        ox=[n for n in name2f if "oxo" in n]
        if ien and ox: rm.update(ox)
    return ",".join([e for e in ents if re.match(r"(.+?)(?:\(|$)",e).group(1) not in rm])

def determine_select_acyl(row,acyl_df):
    try: sub_no=int(row.get("Sub No.",0))
    except: sub_no=0
    if sub_no==0: return "none"
    subs=[c.strip() for c in str(row.get("Sub Carbon","")).split(",") if c.strip().startswith("C")]
    raw=[p.strip() for p in str(row.get("Acyl part","")).split(",") if p.strip()]
    parsed=[(re.match(r"(.+?)\(([\d\.]+)\)",p).group(1),float(re.match(r".+?\(([\d\.]+)\)",p).group(1))) if "(" in p else (p,None) for p in raw]
    corr={2:"acyl",4:"isobutyryl",5:"2-methylbutyryl",7:"benzoyl"}
    def cands(parts):
        out=[]
        for name,score in parts:
            sel=acyl_df[acyl_df["Identification"]==name]
            if not sel.empty:
                out.append((score or 0,name,count_element(sel.iloc[0]["Acyl Formula"],"C")))
        return out
    # two distinct
    if sub_no==2 and len(subs)==2:
        n1,n2=int(subs[0][1:]),int(subs[1][1:]); total=n1+n2
        allc=cands(parsed)
        pairs=[((s1,n1),(s2,n2)) for s1,n1,c1 in allc for s2,n2,c2 in allc if n1!=n2 and c1+c2==total]
        if pairs:
            best=max(pairs,key=lambda pr: pr[0][0]+pr[1][0])
            return ",".join(dict.fromkeys([best[0][1],best[1][1]]))
    # identical double
    if sub_no==2 and len(subs)==1:
        n=int(subs[0][1:]); total=2*n
        allc=cands(parsed)
        pairs=[((s1,n1),(s2,n2)) for (s1,n1,c1),(s2,n2,c2) in combinations(allc,2) if c1+c2==total]
        if pairs:
            best=max(pairs,key=lambda pr: pr[0][0]+pr[1][0])
            return ",".join(dict.fromkeys([best[0][1],best[1][1]]))
    # per-carbon
    valid=[]
    for sc in subs:
        n=int(sc[1:]) if sc[1:].isdigit() else None
        if n:
            cands_i=[(score,name) for name,score in parsed for _,row in acyl_df[acyl_df["Identification"]==name].iterrows() if count_element(row["Acyl Formula"],"C")==n]
            if cands_i:
                valid.append(max(cands_i,key=lambda x:x[0])[1])
            elif corr.get(n):
                valid.append(corr[n])
    return ",".join(dict.fromkeys(valid)) or "none"

def double_check_acyl_parts(row,acyl_df):
    corr_map={2:"acyl",4:"isobutyryl",5:"2-methylbutyryl",7:"benzoyl"}
    orig=[p.strip() for p in str(row.get("Correct Acyl part","")).split(",") if p.strip()]
    pat=re.findall(r"C(\d+)\*(\d+)",row.get("Sub_An",""))
    if not pat: return ",".join(orig)
    carbon_to_ident={}
    for ident in orig:
        sel=acyl_df[acyl_df["Identification"]==ident]
        if not sel.empty:
            c=count_element(sel.iloc[0]["Acyl Formula"],"C")
            carbon_to_ident[c]=ident
    res=[]
    for c_str,count_str in pat:
        c_num,cnt=int(c_str),int(count_str)
        ident=carbon_to_ident.get(c_num, corr_map.get(c_num))
        if ident:
            res+= [ident]*cnt
    return ",".join(res) or ",".join(orig)

def _sort_key(s: str):
    m = re.match(r'^([A-Za-z]+)(\d+)$', s)
    if m:
        return (m.group(1), int(m.group(2)))
    else:
        return (s, 0)

def generate_composite_key_list(row: pd.Series, id_cols: list) -> list:
    parts = []
    for col in id_cols:
        val = row.get(col)
        if isinstance(val, str):
            v = val.strip()
            if v and v.lower() != 'none':
                parts.append(v)
    return sorted(parts, key=_sort_key)

def build_db_annotation_mapping(db_path: str) -> dict:
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [r[0] for r in cursor.fetchall()]
    mapping = {}
    for tbl in tables:
        df_tbl = pd.read_sql_query(f"SELECT * FROM '{tbl}'", conn)
        part_cols = [c for c in df_tbl.columns if c.endswith(' parts')]
        if not part_cols:
            continue
        for _, db_row in df_tbl.iterrows():
            parts = []
            for col in part_cols:
                val = str(db_row.get(col, '')).strip()
                if not val or val.lower() in ('none', 'others'):
                    continue
                parts.append(val)
            key = tuple(sorted(parts, key=_sort_key))
            mapping.setdefault(key, []).append(
                (db_row.get('New SMILES'), db_row.get('Number'))
            )
    conn.close()
    return mapping

def annotate_df_with_db(df: pd.DataFrame, db_path: str) -> pd.DataFrame:
    mapping = build_db_annotation_mapping(db_path)
    id_cols = [c for c in df.columns if c.endswith('_id')]

    expanded_rows = []
    for _, row in df.iterrows():
        key = tuple(generate_composite_key_list(row, id_cols))
        matches = mapping.get(key)
        if matches:
            for smiles, num in matches:
                new_row = row.copy()
                new_row['Annotation SMILES'] = smiles
                new_row['Annotation ID']      = num
                expanded_rows.append(new_row)
        else:
            new_row = row.copy()
            new_row['Annotation SMILES'] = None
            new_row['Annotation ID']      = None
            expanded_rows.append(new_row)

    df_out = pd.DataFrame(expanded_rows).reset_index(drop=True)
    return df_out

def add_aglycone_type(df: pd.DataFrame) -> pd.DataFrame:
    def f(r):
        for c in ("A_id","B_id","C_id"):
            if str(r[c]).strip().lower() in ("","none"): return "none"
        parts=[r["A_id"],r["B_id"],r["C_id"]]
        m=str(r["M_id"]).strip()
        if m.lower() not in ("","none"): parts.append(m)
        return "".join(parts)
    df["Aglycone Type"]=df.apply(f,axis=1)
    cols=list(df.columns)
    ai=cols.index("Aglycone Type"); si=cols.index("Annotation SMILES")
    cols.insert(si,cols.pop(ai))
    return df[cols]

# --- Scoring function for each row ---
def get_scores_by_parts(row):
    score_cols = []
    type = row.get('type', '').strip()

    if type == 'MD':
        # M part scoring
        m_part = row.get('M part', '')
        if 'C10 ring' in m_part:
            if 'substituent' in m_part or 'substituents' in m_part:
                score_cols.append('C10 ion_score')
            else:
                score_cols.append('C10 loss_score')
        elif 'C9 ring' in m_part:
            score_cols.append('C9 loss_score')
        elif 'C14 ring' in m_part:
            score_cols.append('C14 loss_score')

        # A part scoring
        a_part = row.get('A part', '')
        if a_part == 'bicyclo [2.2.1] heptane ring':
            score_cols.extend(['C14 ion_score', 'C19 ion_score'])
        elif a_part == '3,4-seco':
            score_cols.extend(['C15 ion_score', 'C16 ion_score'])
        elif a_part in ('1-alkyl-3-ol', '1-alkyl-3-one'):
            score_cols.append('C20 ion_score')
        elif a_part in ('1-alkyl-2-ol-3-one', '1-alkyl-2-ol-3-ol'):
            score_cols.extend(['C16 ion_score', 'C20 ion_score'])

        # B part scoring
        b_part = row.get('B part', '')
        if b_part == '5-ol-6,7-epoxy':
            score_cols.extend(['C17 ion_score', 'C27 ion_score', 'C3 loss_score'])

        # C-12 and C-18 bonus scores
        if row.get('C-12') == 'with C-12':
            score_cols.append('C19 ion_score')
        if row.get('C-18') == 'with C-18':
            score_cols.extend(['C18 ion_score', 'C18 loss_score'])

    elif type == 'D':
        # A part
        a_part = row.get('A part', '')
        if a_part in ('1,2-en-3-one', '1,2-dihydro-3-ol'):
            score_cols.append('C17 ion_score')
        elif a_part in ('1,2-en-3-ol', '1,2-dihydro-3-one'):
            score_cols.extend(['C17 ion_score', 'C18 ion_score'])

        # B part
        b_part = row.get('B part', '')
        if b_part == '5-ol-6,7-diol':
            score_cols.append('C19 ion_score')
        elif b_part == '5-ol-6,7-epoxy':
            score_cols.append('C17 ion_score')
        elif b_part in ('5-ol-4,6-epoxy', '5-ol-4,7-epoxy'):
            score_cols.append('C15 ion_score')

    # Olefinic scoring
    if row.get('olefinic position') == '2,4-olefinic':
        score_cols.extend(['C5 ion_score', 'C6 ion_score'])

    # Acyl part scoring
    acyl_part = row.get('Acyl part', '')
    loss_score_map = {
        'acyl': 'C2 loss_score',
        'benzoyl': 'C7 loss_score',
        'cinnamoyl': 'C9 loss_score',
        'coumaroyl': 'C9 loss_score',
        'nonanoyl': 'C9 loss_score',
        'feruloyl': 'C10 loss_score',
        'decanoyl': 'C10 loss_score',
        'undecanoyl': 'C11 loss_score',
    }
    if any(keyword in acyl_part for keyword in loss_score_map):
        score_cols.extend(['Acyl A ion_score', 'Acyl B ion_score'])

    for keyword, loss_col in loss_score_map.items():
        if keyword in acyl_part:
            score_cols.append(loss_col)

    return score_cols

# --- Compute total score by summing relevant ion/loss intensities ---
def precise_total_score(row):
    score_cols = get_scores_by_parts(row)
    total = 0
    for col in score_cols:
        if not isinstance(col, str):
            continue
        val = pd.to_numeric(row.get(col, 0), errors='coerce')
        if pd.notna(val):
            total += val
    return total

def _calc_formula(smiles:str)->str:
    if not smiles: return ""
    m=Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcMolFormula(m) if m else ""

def _calc_mw(smiles:str)->float:
    if not smiles: return 0.0
    m=Chem.MolFromSmiles(smiles)
    return round(MolWt(m),4) if m else 0.0

# ------ Mapping tables ------
df_a = pd.DataFrame({
    "A number":["A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11"],
    "A Struc":["1,2-en-3-one","1,2-dihydro-3-one","1,2-en-3-ol","1,2-dihydro-3-ol",
               "1,10-en-3-one","1-alkyl-3-one","1-alkyl-2-ol-3-one","1-alkyl-3-ol",
               "1-alkyl-2-ol-3-ol","3,4-seco","bicyclo [2.2.1] heptane"]
})
df_b = pd.DataFrame({
    "B number":["B1","B2","B3","B4","B5","B6","B7","B8"],
    "B Struc":["5-ol-6,7-epoxy","5-ol-6,7-diol","5-ol-6,7-ene",
               "5-dehydro-6,7-epoxy","5-dehydro-6,7-diol","5-dehydro-6,7-ene",
               "5-ol-4,7-epoxy","5-ol-4,6-epoxy"]
})
df_c = pd.DataFrame({
    "C number":["C1","C2","C3","C4","C5","C6"],
    "C Struc":["9,13,14-orthoester","9,13,14-orthoester","9,13,14-triol",
               "9,13,14-triol","9,13,14-orthoester","9,13,14-orthoester"],
    "C-12Struc":["no C-12","with C-12","no C-12","with C-12","no C-12","with C-12"],
    "C-18Struc":["no C-18","no C-18","no C-18","no C-18","with C-18","with C-18"]
})
df_m = pd.DataFrame({
    "M number":["M1","M2","M3","M4","M5","M6","M7","M8","M9"],
    "M Struc":[["C9 ring with no substituents"],["C10 ring with no substituents"],
               ["C10 ring with one substituents"],["C10 ring with two substituents","C10 ring with one substituents and a hydroxyl group"],
               ["C10 ring with three substituents","C10 ring with two substituents and a hydroxyl group"],
               ["C10 ring with four substituents","C10 ring with three substituents and a hydroxyl group"],
               ["C14 ring with an olenfinic"],["C16 ring with an olenfinic"],["C16 ring with one substituents"]]
})

def map_a_number(row):
    s = str(row['A part']).strip()
    m = df_a[df_a['A Struc'].str.strip() == s]
    return m['A number'].iat[0] if not m.empty else 'none'

def map_b_number(row):
    s = str(row['B part']).strip()
    m = df_b[df_b['B Struc'].str.strip() == s]
    return m['B number'].iat[0] if not m.empty else 'none'

def map_c_number(row):
    s = str(row['C part']).strip()
    c12 = str(row['C-12']).strip()
    c18 = str(row['C-18']).strip()
    m = df_c[
        (df_c['C Struc'].str.strip()   == s) &
        (df_c['C-12Struc'].str.strip() == c12) &
        (df_c['C-18Struc'].str.strip() == c18)
    ]
    return m['C number'].iat[0] if not m.empty else 'none'

def map_m_number(row):
    s = str(row['M part']).strip()
    m = df_m[df_m['M Struc'].apply(lambda lst: s in lst)]
    return m['M number'].iat[0] if not m.empty else 'none'

def _read_csv_with_fallback(path):
    raw = open(path, 'rb').read(4096)
    guess = chardet.detect(raw)['encoding'] or 'utf-8'
    for enc in [guess, 'utf-8', 'latin-1']:
        try:
            return pd.read_csv(path, encoding=enc)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, encoding='utf-8', errors='ignore')

def annotation_process_data(query_csv, acyl_csv, pseudo_db, output_csv):
    df      = _read_csv_with_fallback(query_csv)
    acyl_df = _read_csv_with_fallback(acyl_csv)

    df["IHD"]=df.apply(determine_ihd,axis=1)
    df["A part"]=df.apply(determine_a_parts,axis=1)
    df["B part"]=df.apply(lambda r:determine_b_parts(r,df),axis=1)
    df["C part"]=determine_c_parts(df)
    df["C-12"]=df.apply(determine_c12,axis=1)
    df["C-18"]=df.apply(lambda r:determine_c18(r,df),axis=1)
    df["olefinic position"]=df.apply(determine_olefinic_part,axis=1)

    expanded=[]
    for _,r in df.iterrows():
        for m in determine_m_parts(r):
            rr=r.copy(); rr["M part"]=m
            expanded.append(rr)
    df_exp=pd.DataFrame(expanded)

    sub_rows=[]
    for _,r in df_exp.iterrows():
        sub_rows.extend(determine_substitution(r))
    df_exp=pd.DataFrame(sub_rows)

    df_exp["Acyl part"]=df_exp.apply(lambda r:determine_acyl_parts(r,acyl_df),axis=1)
    df_exp["Acyl part"]=df_exp.apply(lambda r:refine_acyl_parts_by_olefinic(r,acyl_df),axis=1)
    df_exp["Correct Acyl part"]=df_exp.apply(lambda r:determine_select_acyl(r,acyl_df),axis=1)
    df_exp["Correct Acyl part"]=df_exp.apply(lambda r:double_check_acyl_parts(r,acyl_df),axis=1)

    df_exp["A_id"] = df_exp.apply(map_a_number, axis=1)
    df_exp["B_id"] = df_exp.apply(map_b_number, axis=1)
    df_exp["C_id"] = df_exp.apply(map_c_number, axis=1)
    df_exp["M_id"] = df_exp.apply(map_m_number, axis=1)

    r_map={}
    for _,r in acyl_df.iterrows():
        for ident in re.split(r"[,/]",str(r["Identification"])):
            if ident.strip(): r_map[ident.strip()]=r.get("R parts","")
    r_series=df_exp["Correct Acyl part"].apply(lambda s:[r_map.get(i,"none") for i in str(s).split(",")])
    maxr=r_series.map(len).max()
    for i in range(maxr):
        df_exp[f"R{i+1}_id"]=r_series.map(lambda lst: lst[i] if len(lst)>i else "none")

    df_exp["Total Score"]=df_exp.apply(precise_total_score,axis=1)
    df_exp=df_exp.sort_values("Total Score",ascending=False)
    key_cols=["A part","B part","C part","C-12","C-18","M part","Sub Carbon","row ID"]
    df_exp=df_exp.drop_duplicates(subset=key_cols)

    df_exp=annotate_df_with_db(df_exp,pseudo_db)
    df_exp=add_aglycone_type(df_exp)
    df_exp["Annotation Formula"]=df_exp["Annotation SMILES"].apply(_calc_formula)
    df_exp["Annotation MW"]=df_exp["Annotation SMILES"].apply(_calc_mw)

    df_exp["TopK_num"]=df_exp.groupby("row ID")["Total Score"].rank(method="dense",ascending=False).astype(int)
    df_exp["Top K"]=df_exp["TopK_num"].map(lambda x:f"Top{x}")
    df_exp["confidence score"]=df_exp.groupby("row ID")["Total Score"].transform(lambda s:(s/s.max()).round(3))

    df_exp.to_csv(output_csv,index=False)
    base,ext=os.path.splitext(output_csv)
    top1=f"{base}_top1{ext}"
    df_exp[df_exp["Top K"]=="Top1"].drop_duplicates().to_csv(top1,index=False)
    zip_path=f"{base}.zip"
    with zipfile.ZipFile(zip_path,"w") as z:
        z.write(output_csv,os.path.basename(output_csv))
        z.write(top1,os.path.basename(top1))
    return zip_path

""" Tab 8 related Functions """
@app.callback(
    Output('visualization-query-csv-file-path', 'children'),
    Input('visualization-query-upload-csv', 'contents'),
    State('visualization-query-upload-csv', 'filename')
)
def update_visualization_query_csv_file_path(csv_contents, csv_filenames):
    if csv_contents and csv_filenames:
        print(f"Received CSV filenames: {csv_filenames}")
        csv_paths = [save_tab8_file(name, content) for name, content in zip(csv_filenames, csv_contents)]
        print(f"CSV paths saved: {csv_paths}")
        csv_path_text = f"Uploaded CSV files: {', '.join(csv_filenames)}"
    else:
        print("No CSV file uploaded")
        csv_path_text = "No CSV file uploaded"
    return csv_path_text

@app.callback(
    Output('scatter-plot', 'figure'),
    Input('run-visualization', 'n_clicks'),
    State('visualization-query-upload-csv', 'filename'),
    Input('scatter-plot', 'clickData')
)
def generate_scatter_plot(n_clicks, csv_filenames, clickData):
    if n_clicks == 0 or not csv_filenames:
        return {}

    file_path = os.path.join(UPLOAD_DIR, 'tab_8_uploads', csv_filenames[0])
    df = pd.read_csv(file_path)

    if 'index' not in df.columns:
        unique_row_ids = df["row ID"].dropna().unique()
        row_id_to_index = {rid: idx + 1 for idx, rid in enumerate(sorted(unique_row_ids))}
        df["index"] = df["row ID"].map(row_id_to_index)

    peak_area_cols = [col for col in df.columns if col.endswith("Peak area")]

    long_df = df.melt(
        id_vars=["index", "row ID", "row m/z", "row retention time", "elemental composition", "adduct ion", "identification"],
        value_vars=peak_area_cols,
        var_name="Plant Sample",
        value_name="Concentration"
    )

    long_df["Concentration"] = pd.to_numeric(long_df["Concentration"], errors="coerce")
    long_df = long_df[long_df["Concentration"] > 0]

    def base_color(ident):
        return "possibly undescribed" if ident == "possibly undescribed" else "possibly known"

    long_df["Color"] = long_df["identification"].map(base_color)

    if clickData:
        selected_index = clickData["points"][0]["x"]
        selected_sample = clickData["points"][0]["y"]
        mask = (long_df["index"] == selected_index) & (long_df["Plant Sample"] == selected_sample)
        long_df.loc[mask, "Color"] = "highlight"

    color_map = {
        "possibly undescribed": "red",
        "possibly known": "grey",
        "highlight": "yellow"
    }

    fig = px.scatter_3d(
        long_df,
        x="index",
        y="Plant Sample",
        z="row retention time",
        size="Concentration",
        color="Color",
        color_discrete_map=color_map,
        hover_data={
            "row ID": True,
            "row m/z": True,
            "row retention time": True,
            "elemental composition": True,
            "adduct ion": True
        },
        size_max=40
    )

    fig.update_layout(
        scene=dict(
            xaxis_title="Peak Index",
            yaxis_title="Plant Sample",
            zaxis_title="Retention Time (min)"
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )

    return fig

@app.callback(
    Output('filtered-table-container', 'children'),
    [Input('scatter-plot', 'clickData'),
     Input('visualization-query-upload-csv', 'filename')]
)
def update_filtered_table(clickData, csv_filenames):
    if not csv_filenames:
        return html.Div("No CSV uploaded.")

    file_path = os.path.join(UPLOAD_DIR, 'tab_8_uploads', csv_filenames[0])
    df = pd.read_csv(file_path)

    if 'index' not in df.columns:
        unique_row_ids = df['row ID'].dropna().unique()
        row_id_to_index = {rid: idx + 1 for idx, rid in enumerate(sorted(unique_row_ids))}
        df['index'] = df['row ID'].map(row_id_to_index)

    if clickData:
        selected_index = clickData['points'][0]['x']
        df = df[df['index'] == selected_index]

    columns_to_show = [
        "row ID", "row m/z", "row retention time", "type", "elemental composition", "identification",
        "same composition", "adduct ion", "same peak",
        "A part", "B part", "C part", "C-12", "C-18", "M part",
        "Correct Acyl part", "Sub No.", "Sub_An", "Aglycone Type",
        "Top K", "confidence score"
    ]
    available_columns = [col for col in columns_to_show if col in df.columns]

    return dash_table.DataTable(
        data=df[available_columns].to_dict("records"),
        columns=[{"name": col, "id": col} for col in available_columns],
        style_table={
            'overflowX': 'auto',
            'overflowY': 'auto',
            'maxHeight': '500px',
            'border': '1px solid #ccc'
        },
        style_cell={
            'minWidth': '100px',
            'maxWidth': '200px',
            'whiteSpace': 'normal',
            'textAlign': 'left',
            'fontFamily': 'Arial',
            'fontSize': '14px'
        },
        page_size=20,
        filter_action="native",
        sort_action="native",
        style_header={
            'backgroundColor': '#f2f2f2',
            'fontWeight': 'bold'
        }
    )


@app.callback(
    [Output('structure-gallery', 'children'),
     Output('structure-modal-container', 'children')],
    Input('scatter-plot', 'clickData'),
    State('visualization-query-upload-csv', 'filename')
)
def update_structure_gallery(clickData, csv_filenames):
    if not csv_filenames:
        return html.Div("No structure to display."), []

    if not clickData:
        return html.Div("Click on a data point to view annotated structures."), []

    file_path = os.path.join(UPLOAD_DIR, 'tab_8_uploads', csv_filenames[0])
    df = pd.read_csv(file_path)

    if 'index' not in df.columns:
        unique_row_ids = df['row ID'].dropna().unique()
        row_id_to_index = {rid: idx + 1 for idx, rid in enumerate(sorted(unique_row_ids))}
        df['index'] = df['row ID'].map(row_id_to_index)

    selected_index = clickData["points"][0]["x"]
    filtered_df = df[df["index"] == selected_index]
    filtered_df = filtered_df[filtered_df["Annotation SMILES"].notna()]

    cards = []
    modals = []

    for i, row in filtered_df.iterrows():
        smiles = row["Annotation SMILES"]
        annotation_id = row.get("Annotation ID", "N/A")
        formula = row.get("Annotation Formula", "N/A")
        mw = row.get("Annotation MW", "N/A")

        # 处理 confidence score 保留三位小数
        score_raw = row.get("confidence score")
        try:
            confidence = f"{float(score_raw):.3f}"
        except (ValueError, TypeError):
            confidence = "N/A"

        mol_img = smiles_to_image(smiles)
        if not mol_img:
            continue

        # 卡片区域
        cards.append(html.Div([
            html.Div([
                html.Div(f"{annotation_id}", style={"fontWeight": "bold", "fontSize": "13px"}),
                html.Div(f"{formula} | MW={mw}", style={"fontSize": "12px"}),
                html.Div(f"Confidence: {confidence}", style={"fontSize": "12px", "color": "#444"})
            ], style={"textAlign": "center", "marginBottom": "8px"}),

            html.Img(
                src=mol_img,
                id={"type": "img", "index": i},
                n_clicks=0,
                style={
                    "maxWidth": "200px", "height": "200px",
                    "display": "block", "margin": "auto",
                    "cursor": "pointer"
                }
            )
        ], style={
            "border": "1px solid #ccc",
            "borderRadius": "8px",
            "padding": "10px",
            "width": "220px",
            "backgroundColor": "#ffffff",
            "boxShadow": "0 2px 4px rgba(0,0,0,0.1)"
        }))

        # Modal 放大图
        modals.append(
            dbc.Modal([
                dbc.ModalHeader(f"{annotation_id} | {formula} | MW={mw}"),
                dbc.ModalBody(html.Img(src=mol_img, style={"width": "100%"})),
            ],
                id={"type": "modal", "index": i},
                is_open=False,
                size="xl")
        )

    return cards, modals

@app.callback(
    Output({'type': 'modal', 'index': MATCH}, 'is_open'),
    Input({'type': 'img', 'index': MATCH}, 'n_clicks'),
    State({'type': 'modal', 'index': MATCH}, 'is_open')
)
def toggle_modal(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

@app.callback(
    Output("summary-bar-count", "children"),
    Output("summary-bar-area", "children"),
    Input('run-visualization', 'n_clicks'),
    State('visualization-query-upload-csv', 'filename')
)
def update_summary_charts(n_clicks, csv_filenames):
    if not csv_filenames:
        return None, None
    file_path = os.path.join(UPLOAD_DIR, 'tab_8_uploads', csv_filenames[0])
    df = pd.read_csv(file_path)

    if "index" not in df.columns:
        unique_row_ids = df["row ID"].dropna().unique()
        df["index"] = df["row ID"].map({rid: idx + 1 for idx, rid in enumerate(sorted(unique_row_ids))})

    count_fig, area_fig = generate_summary_bar_charts_plotly(df)
    return count_fig, area_fig

def generate_summary_bar_charts_plotly(df):
    peak_cols = [c for c in df.columns if c.endswith("Peak area")]
    if not peak_cols:
        return None, None

    df_long = df.melt(
        id_vars=["index", "identification"],
        value_vars=peak_cols,
        var_name="Sample",
        value_name="Area"
    )

    df_long["Area"] = pd.to_numeric(df_long["Area"], errors="coerce")
    df_long = df_long[df_long["Area"] > 0]
    df_long["Sample"] = df_long["Sample"].str.replace(" Peak area", "", regex=False)
    df_long["count"] = 1

    count_summary = df_long.groupby(["Sample", "identification"])["count"].sum().reset_index()
    count_summary = count_summary.pivot(index="Sample", columns="identification", values="count").fillna(0)

    area_summary = df_long.groupby(["Sample", "identification"])["Area"].sum().reset_index()
    area_summary = area_summary.pivot(index="Sample", columns="identification", values="Area").fillna(0)

    for col in ["possibly undescribed", "possibly known"]:
        if col not in count_summary.columns:
            count_summary[col] = 0
        if col not in area_summary.columns:
            area_summary[col] = 0

    samples = count_summary.index.tolist()
    known_counts = count_summary.get("possibly known").values
    new_counts = count_summary.get("possibly undescribed").values
    known_areas = area_summary.get("possibly known").values
    new_areas = area_summary.get("possibly undescribed").values

    fig_count = go.Figure(data=[
        go.Bar(name="possibly known", x=samples, y=known_counts, marker=dict(color='lightgray')),
        go.Bar(name="possibly undescribed", x=samples, y=new_counts, marker=dict(color='red'))
    ])
    fig_count.update_layout(
        barmode='stack',
        title="Possibly Undescribed Compound Count",
        yaxis_title="Count",
        xaxis_tickangle=60,
        height=400
    )

    fig_area = go.Figure(data=[
        go.Bar(name="possibly known", x=samples, y=known_areas, marker=dict(color='lightgray')),
        go.Bar(name="possibly undescribed", x=samples, y=new_areas, marker=dict(color='red'))
    ])

    fig_area.update_layout(
        barmode='stack',
        title="New Compound Area",
        yaxis_title="Total Area",
        xaxis_tickangle=60,
        height=400
    )

    return dcc.Graph(figure=fig_count), dcc.Graph(figure=fig_area)

def smiles_to_image(smiles, size=(250, 250)):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=size)
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"
    else:
        return None

def save_tab8_file(name, content, tab_folder='tab_8_uploads'):
    print(f"Saving file: {name} to folder: {tab_folder}")
    folder_path = os.path.join(UPLOAD_DIR, tab_folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    file_path = os.path.join(folder_path, name)
    with open(file_path, 'wb') as f:
        f.write(decoded)
    print(f"File saved at: {file_path}")
    return file_path

# =========================
#  Main Entrance
# =========================
# if __name__ == '__main__':
    # app.run_server(debug=True)
