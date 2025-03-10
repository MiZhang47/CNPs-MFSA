"""
Last updated: 2025/02
Complex natural products (CNPs) Structure Annotation Tool constructed based on a modular fragmentation-based structural assembly (MFSA) strategy.
"""

# =========================
#  Import Python libraries
# =========================
# Import libraries for data processing and file handling
import io  # Used for handling file streams
import os  # Provides a way of using operating system-dependent functionalities
import re  # Regular expression operations for string matching
import uuid  # Generates unique identifiers for various objects
import base64  # Encoding and decoding operations for binary data
import zipfile  # Handles zip file creation and extraction
import pandas as pd  # Data manipulation and analysis library
import numpy as np  # Numerical computation library
from typing import Dict, Tuple  # Type hinting for function signatures
from concurrent.futures import ThreadPoolExecutor  # Multithreading support to handle concurrent operations

# Import external libraries for chemical analysis and handling mass spectrometry data
import sqlite3  # Provides SQLite database operations
from pyteomics import mass, mgf  # Pyteomics library for handling mass spectrometry files and calculations
from matchms import Spectrum, calculate_scores  # MatchMS library for handling mass spectrometry spectra and similarity calculations
from matchms.similarity import CosineGreedy  # Cosine similarity scoring algorithm from MatchMS
from rdkit import Chem  # RDKit library for chemical informatics and computational chemistry
from rdkit.Chem import Draw  # RDKit's drawing utilities for visualizing molecular structures

# Import libraries for building and managing a web-based dashboard using Dash
import dash  # Dash framework for building web applications
from dash import dcc, html, dash_table  # Dash core components, HTML layout, and tables
import dash_bootstrap_components as dbc  # Bootstrap components for Dash
from dash.dependencies import Input, Output, State, MATCH, ALL  # Dash dependencies for interactive components
import plotly.express as px  # Plotly Express for creating visualizations
from flask import send_file  # Flask utility for sending files in response to HTTP requests


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

                    # Upload section for the Skeleton CSV file (single file)
                    html.Div([
                        html.H4('Import Skeleton CSV file:', style={'fontFamily': 'Arial'}),  # Header for Skeleton CSV file upload
                        dcc.Upload(
                            id='annotation-skeleton-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=False  # Only one Skeleton CSV file allowed
                        ),
                        html.Div(id='annotation-skeleton-csv-file-path',
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
                # Container for uploading query CSV files
                html.Div([
                    html.Div([
                        html.H4('Import Query CSV file (multiple files):', style={'fontFamily': 'Arial'}),  # Header for CSV file upload
                        dcc.Upload(
                            id='visualization-query-upload-csv',
                            children=html.Div('Drop or Select CSV File'),
                            style={
                                'width': '90%', 'height': '50px', 'lineHeight': '50px',
                                'borderWidth': '1px', 'borderStyle': 'dashed',
                                'borderRadius': '5px', 'textAlign': 'center',
                                'margin': '10px', 'fontFamily': 'Arial'
                            },
                            multiple=True  # Allow multiple CSV files
                        ),
                        html.Div(id='visualization-query-csv-file-path',
                                 style={'fontFamily': 'Arial', 'marginTop': '10px', 'color': '#007bff'})  # Placeholder for file paths
                    ], style={'width': '48%'}),
                ], style={'display': 'flex', 'justifyContent': 'space-between'}),  # Align upload section

                # Button to trigger visualization
                html.Button('Run Visualization', id='run-visualization', n_clicks=0,
                            style={'marginTop': '20px', 'fontFamily': 'Arial', 'padding': '10px 20px',
                                   'backgroundColor': '#007bff', 'color': 'white', 'border': 'none',
                                   'cursor': 'pointer', 'fontSize': '16px'}),
                dbc.Card(
                    [
                        dbc.CardBody(
                            [
                                dbc.Row([
                                    dbc.Col(
                                        dcc.Graph(
                                            id='scatter-plot',
                                            style={"height": "700px"}
                                        ),
                                        width=8
                                    ),
                                    dbc.Col([
                                        html.H3("Skeleton Structure",
                                                style={"font-size": "18px"}),
                                        html.Img(id='skeleton-image',
                                                 src='',
                                                 style={"max-width": "200px"}),
                                        html.Br(),
                                        html.H3("Acyl Structure R1",
                                                style={"font-size": "18px"}),
                                        html.Img(id='acyl-image-1',
                                                 src='',
                                                 style={"max-width": "200px"}),
                                        html.Br(),
                                        html.H3("Acyl Structure R2",
                                                style={"font-size": "18px"}),
                                        html.Img(id='acyl-image-2',
                                                 src='',
                                                 style={"max-width": "200px"})
                                    ], width=3)
                                ])
                            ]
                        )
                    ],
                    style={
                        'border': '1px solid #ccc',
                        'borderRadius': '5px',
                        'padding': '10px',
                        'marginTop': '50px',
                        'marginBottom': '150px',
                        'marginLeft': '150px',
                        'marginRight': '150px',
                    }
                )
            ])
        ], style={  # Style settings for the 'Visualization' tab
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

            # Add score columns for each type and filter non-zero score rows
            for i, type_name in enumerate(type_names):
                focal_csv[f'{type_name} Score'] = focal_csv.index.map(lambda x: scores_dict[type_name].get(x, 0))

            score_columns = [f'{type_name} Score' for type_name in type_names]
            focal_csv = focal_csv[(focal_csv[score_columns] != 0).any(axis=1)]
            print(f"Rows after filtering for non-zero scores: {focal_csv.shape[0]}")

            # Determine the type associated with the highest score
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
                if intensity_value > 0.1:
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
                dcc.Input(id={'type': 'carbon-range-max2', 'index': new_index + 1}, type='text', placeholder='C30',
                          style={'width': '10%', 'marginRight': '10px', 'fontFamily': 'Arial'}),
                html.Span('to', style={'marginRight': '10px', 'fontFamily': 'Arial'}),
                dcc.Input(id={'type': 'carbon-range-min2', 'index': new_index + 1}, type='text', placeholder='C27',
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
    State({'type': 'carbon-range-min2', 'index': ALL}, 'value')
)
def run_neutral_loss_extraction(n_clicks, csv_filenames, csv_contents, mgf_filenames, mgf_contents, db_filename,
                                db_content, charge, ppm, target_losses, loss_names,
                                carbon_range_max1, carbon_range_min1, carbon_range_max2, carbon_range_min2):
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
        carbon_range_max1, carbon_range_min1, carbon_range_max2, carbon_range_min2: Carbon range specifications.

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
                carbon_range_min2=carbon_range_min2
            )
        except Exception as e:
            return f"Error processing files: {str(e)}"

    return f"Extraction complete for {len(csv_paths)} CSV file(s)."


def neutralloss_process_files(input_path_csv, input_path_mgf, output_path_csv, db_path, target_losses, charge, ppm,
                              loss_names,
                              carbon_range_max1, carbon_range_min1, carbon_range_max2, carbon_range_min2):
    """
    Process CSV and MGF files to find neutral loss matches based on target losses.

    Parameters:
        input_path_csv: Path to the input CSV file.
        input_path_mgf: Path to the input MGF file.
        output_path_csv: Path where the processed CSV file will be saved.
        db_path: Path to the database file containing compounds.
        target_losses: List of target loss formulas.
        charge: Charge state for calculating m/z.
        ppm: Tolerance in PPM for m/z matching.
        loss_names: Names of the loss columns in the output CSV.
        carbon_range_max1, carbon_range_min1, carbon_range_max2, carbon_range_min2: Carbon ranges for filtering.

    Returns:
        None
    """
    print("Starting file processing")

    # Read the input CSV
    df = pd.read_csv(input_path_csv)
    print(f"Read CSV with {len(df)} rows")

    # Check if 'row ID' column exists
    if 'row ID' not in df.columns:
        raise ValueError("'row ID' column not found in the CSV file.")

    df.set_index('row ID', inplace=True)
    df.index = df.index.map(clean_title)

    # Compute m/z ranges for formulas from the DB
    formula_ranges = read_formulas_from_NLdb(db_path, charge, ppm)
    print(f"Formula ranges computed for DB")

    # Read the MGF file and create a dictionary of spectra
    spectra_dict = {}
    with mgf.read(input_path_mgf) as spectra:
        for spectrum in spectra:
            title = clean_title(spectrum['params']['title'])
            spectra_dict[title] = spectrum
        print(f"Processed {len(spectra_dict)} spectra from MGF")

    # Initialize columns for loss names in the dataframe
    for ln in loss_names:
        df[ln] = 'none'

    # Parse all carbon ranges first and save them in a list to avoid repeated parsing
    # Assume that loss_names, target_losses, carbon_range_xxx have the same length
    parsed_carbon_ranges = []
    for i in range(len(loss_names)):
        try:
            cmax1_val = parse_carbon_range(carbon_range_max1[i])
            cmin1_val = parse_carbon_range(carbon_range_min1[i])
            cmax2_val = parse_carbon_range(carbon_range_max2[i])
            cmin2_val = parse_carbon_range(carbon_range_min2[i])
        except Exception as e:
            raise ValueError(f"Invalid carbon range input in group {i+1}: {e}")
        parsed_carbon_ranges.append((cmax1_val, cmin1_val, cmax2_val, cmin2_val))

    # Process each row in the CSV
    for row_id in df.index:
        if row_id not in spectra_dict:
            continue
        spectrum = spectra_dict[row_id]
        precursor_ion = spectrum['params']['pepmass'][0]
        df.at[row_id, 'row m/z'] = precursor_ion

        product_ions = spectrum['m/z array']
        intensities = spectrum['intensity array']

        # Find possible formulas based on m/z and db range
        # find_matching_NLformulas here returns a single match or 'none'
        # But if want to support multiple matches, need to rewrite
        formulas = [find_matching_NLformulas(mz, formula_ranges) for mz in product_ions]

        # Test each group of Target Loss
        for i, (loss_group, ln) in enumerate(zip(target_losses, loss_names)):
            if not loss_group:
                # If no formula is entered in this group, skip
                continue
            cmax1_val, cmin1_val, cmax2_val, cmin2_val = parsed_carbon_ranges[i]

            # Get the formula of carbon number = cmax1_val / cmin1_val
            formulas_max_carbon_1 = [
                f for f in formulas
                if f != 'none' and mass.Composition(f).get('C', 0) == cmax1_val
            ]
            formulas_min_carbon_1 = [
                f for f in formulas
                if f != 'none' and mass.Composition(f).get('C', 0) == cmin1_val
            ]

            # Get the formula of carbon number = cmax2_val / cmin2_val
            formulas_max_carbon_2 = [
                f for f in formulas
                if f != 'none' and mass.Composition(f).get('C', 0) == cmax2_val
            ]
            formulas_min_carbon_2 = [
                f for f in formulas
                if f != 'none' and mass.Composition(f).get('C', 0) == cmin2_val
            ]

            matching_pairs = []

            # The first pair of carbon ranges
            for f1 in formulas_max_carbon_1:
                for f2 in formulas_min_carbon_1:
                    diff = specific_formula_difference(f1, f2)
                    if check_specific_target_loss(diff, loss_group):
                        idx1 = formulas.index(f1)
                        idx2 = formulas.index(f2)
                        intensity_avg = (intensities[idx1] + intensities[idx2]) / 2
                        matching_pairs.append((f1, f2, intensity_avg))

            # The second pair of carbon ranges
            for f1 in formulas_max_carbon_2:
                for f2 in formulas_min_carbon_2:
                    diff = specific_formula_difference(f1, f2)
                    if check_specific_target_loss(diff, loss_group):
                        idx1 = formulas.index(f1)
                        idx2 = formulas.index(f2)
                        intensity_avg = (intensities[idx1] + intensities[idx2]) / 2
                        matching_pairs.append((f1, f2, intensity_avg))

            # Select the highest peak
            if matching_pairs:
                highest_intensity_pair = max(matching_pairs, key=lambda x: x[2])
                best_formula_f1, best_formula_f2 = highest_intensity_pair[0], highest_intensity_pair[1]
                best_loss = None
                for target_loss in loss_group:
                    target_elements = dict(mass.Composition(target_loss))
                    if specific_formula_difference(best_formula_f1, best_formula_f2) == target_elements:
                        best_loss = target_elements
                        break

                if best_loss:
                    best_loss_formula = dict_to_formula(best_loss)
                    df.at[row_id, ln] = best_loss_formula

    df.reset_index().to_csv(output_path_csv, index=False)
    print(f"Processing complete for {len(df)} rows")

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

    ion_sets_params = []
    for i in range(len(ion_names)):
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
    Process MGF and CSV files to detect feature ions and extract matching formulas in parallel.

    Parameters:
        mgf_path (str): Path to the MGF file containing mass spectra.
        csv_path (str): Path to the CSV file containing metadata or row IDs.
        ion_sets_params (list): List of parameters for each ion set including formulas, charge, tolerance, etc.
        full_output_path (str): Path to save the processed CSV output.

    Returns:
        None: Saves the processed CSV data to the specified output path.
    """
    try:
        # Read the CSV file
        df = pd.read_csv(csv_path, encoding='utf-8')
        print(f"CSV Columns: {df.columns.tolist()}")

        # Check and set 'row ID' as index if present
        if 'row ID' in df.columns:
            df.set_index('row ID', inplace=True)
            print(f"'row ID' column found and set as index.")
        else:
            print(f"'row ID' column not found in {csv_path}. Proceeding without setting index.")

        # Read the MGF file and store spectra by title
        spectra_dict = {}
        with mgf.read(mgf_path) as spectra:
            for spectrum in spectra:
                title = int(spectrum['params']['title'])
                spectra_dict[title] = spectrum
        print(f"Total spectra read from MGF: {len(spectra_dict)}")

        # Prepare formulas for each ion set
        specific_formulas = {}
        for ion_set_params in ion_sets_params:
            ion_set_name = ion_set_params['name']
            specific_formulas[ion_set_name] = [x.strip() for x in ion_set_params['formulas'].split(',')]
        print(f"Specific formulas prepared for ion sets: {specific_formulas.keys()}")

        # Function to process individual ion sets
        def process_ion_set(ion_set_params, ion_set_name):
            matched_formulas = {}
            for title, spectrum in spectra_dict.items():
                result = detect_formulas_in_spectra(spectrum, specific_formulas[ion_set_name], ion_set_params)
                if result:
                    matched_formulas.update(result)
            return ion_set_name, matched_formulas

        # Process ion sets in parallel
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_ion_set, ion_set_params, ion_set_params['name']) for ion_set_params in ion_sets_params]

            # Collect results and update the DataFrame
            for future in futures:
                ion_set_name, matched_formulas = future.result()
                df[ion_set_name] = df.index.map(lambda x: matched_formulas.get(x, 'none'))

        # Save the processed DataFrame to a CSV file
        df.reset_index(inplace=True)
        df.to_csv(full_output_path, index=False)
        print(f"Data processed and saved to {full_output_path}")

    except Exception as e:
        print(f"Error during processing: {e}")
        raise

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
    Output('annotation-skeleton-csv-file-path', 'children'),
    Input('annotation-skeleton-upload-csv', 'contents'),
    State('annotation-skeleton-upload-csv', 'filename')
)
def update_annotation_skeleton_csv_file_path(content, filename):
    print("Skeleton CSV upload triggered")
    if content and filename:
        print(f"Skeleton CSV file name: {filename}")
        saved_file = save_uploaded_file(filename, content)
        if saved_file:
            print(f"File saved at: {saved_file}")
            return f"Uploaded Skeleton CSV file: {filename}"
        else:
            print("Failed to upload the file.")
            return "Failed to upload the file."
    print("No CSV file uploaded")
    return "No CSV file uploaded"

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
    State('annotation-skeleton-upload-csv', 'filename'),
    State('annotation-substituent-upload-csv', 'filename'),
)
def run_annotation(n_clicks, query_csv, skeleton_csv, substituent_csv):
    print(f"Run annotation clicked {n_clicks} times")
    if n_clicks > 0:
        # File paths for query, skeleton, and substituent files
        query_file_path = os.path.join(UPLOAD_DIR, 'tab_7_uploads', query_csv[0])
        skeleton_file_path = os.path.join(UPLOAD_DIR, 'tab_7_uploads', skeleton_csv)
        substituent_file_path = os.path.join(UPLOAD_DIR, 'tab_7_uploads', substituent_csv[0])

        print(
            f"Processing annotation with Query: {query_file_path}, Skeleton: {skeleton_file_path}, Substituent: {substituent_file_path}")

        # Define output file path
        original_filename = query_csv[0]
        base_name, ext = os.path.splitext(original_filename)
        output_filename = f"{base_name}-an{ext}"
        output_file_path = os.path.join(UPLOAD_DIR, 'tab_7_uploads', output_filename)

        # Execute annotation process
        annotation_process_data(query_file_path, substituent_file_path, skeleton_file_path, output_file_path)
        print(f"Annotation completed. Output saved to {output_file_path}")

        return f"Annotation completed. Output saved to {output_file_path}"
    return "Click Run Annotation to start processing."

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

def count_element(molecular_formula, element):
    pattern = re.compile(rf'{element}(\d*)')
    match = pattern.search(molecular_formula)
    if match:
        return int(match.group(1) or 1)
    return 0

def annotation_process_data(path_csv, acyl_path_csv, skeleton_db_path_csv, output_path_csv):
    print("Starting annotation process")
    # Load input dataframes
    df = pd.read_csv(path_csv)
    acyl_df = pd.read_csv(acyl_path_csv)
    skeleton_db = pd.read_csv(skeleton_db_path_csv)
    print(f"Loaded query CSV with {len(df)} rows")
    print(f"Loaded acyl CSV with {len(acyl_df)} rows")
    print(f"Loaded skeleton DB CSV with {len(skeleton_db)} rows")

    def determine_a_parts(row):
        if row['type'] == 'D':
            o_count = count_element(row['C17 ion'], 'O')
            h_count = count_element(row['C17 ion'], 'H')
            o_count_c18 = count_element(row.get('C18 ion', ''), 'O')
            if o_count == 2:
                if h_count == 16:
                    return '1,2-en-3-one'
                elif h_count == 20:
                    return '1,2-dihydro-3-ol'
                elif h_count == 18:
                    if o_count_c18 == 1:
                        return '1,2-en-3-ol'
                    elif o_count_c18 == 2:
                        return '1,2-dihydro-3-one'
            elif o_count == 3:
                if h_count == 16:
                    return '1,2-en-3-one'
                elif h_count == 20:
                    return '1,2-dihydro-3-ol'
                elif h_count == 18:
                    if o_count_c18 == 1:
                        return '1,2-en-3-ol'
                    elif o_count_c18 == 2:
                        return '1,2-dihydro-3-one'
        elif row['type'] == 'MD':
            if row['C15 ion'] != 'none' and row.get('C3 loss', '') != 'none':
                return 'bicyclo [2.2.1] heptane ring'
            elif row['C15 ion'] == 'none':
                return '1-alkyl-3-ol'
        return 'others'

    df['A part'] = df.apply(determine_a_parts, axis=1)

    def determine_b_parts(row):
        if row['type'] == 'D':
            return row['B part']
        elif row['type'] == 'MD':
            c_count_c27 = count_element(row.get('C27 ion', ''), 'C')
            if c_count_c27 == 27 or row.get('C3 loss', '') != 'none':
                return '6,7-epoxy'
            elif row.get('C27 ion', '') == 'none' or row.get('C3 loss', '') == 'none':
                return '6,7-epoxy'
    df['B part'] = df.apply(determine_b_parts, axis=1)

    def determine_c_parts(row):
        if row['type'] == 'D' and row['adduct ion'] == '[M+H-H2O]+':
            return '9,13,14-triol'
        elif row['type'] == 'MD':
            return '9,13,14-orthoester'
        else:
            return '9,13,14-orthoester'
    df['C part'] = df.apply(determine_c_parts, axis=1)

    def determine_c12(row):
        if row['type'] == 'D':
            o_count = count_element(row['C17 ion'], 'O')
            if o_count == 2:
                return 'no C-12'
            elif o_count == 3:
                return 'with C-12'
        return 'others'
    df['C-12'] = df.apply(determine_c12, axis=1)

    def determine_m_parts(row):
        if row['type'] == 'MD':
            C10_h_count = count_element(row.get('C10 loss', ''), 'H')
            if C10_h_count == 18:
                return 'C10 ring with no substituents'
            elif C10_h_count == 16:
                return 'C10 ring with one substituents'
            elif C10_h_count == 14:
                return 'C10 ring with two substituents'
            elif C10_h_count == 12:
                return 'C10 ring with three substituents'
            elif C10_h_count == 10:
                return 'C10 ring with four substituents'
            C9_h_count = count_element(row.get('C9 loss', ''), 'H')
            if C9_h_count == 16:
                return 'C9 ring with no substituents'
            C14_h_count = count_element(row.get('C14 loss', ''), 'H')
            if C14_h_count == 24:
                return 'C14 ring with an olenfinic'
        return 'others'
    df['M part'] = df.apply(determine_m_parts, axis=1)

    def determine_acyl_parts(row):
        acyl_parts = []
        if row['C2 loss'] == 'C2H4O2':
            acyl_parts.append('acyl')

        if row['Acyl ion'] != 'none':
            acyl_ions = row['Acyl ion'].split(',')
            acyl_formulas = acyl_df['Acyl Formula'].tolist()
            identifications = acyl_df['Identification'].tolist()

            for ion in acyl_ions:
                matches = [identifications[i] for i, formula in enumerate(acyl_formulas) if
                           formula.strip() == ion.strip()]
                if matches:
                    acyl_parts.extend(matches)

        if not acyl_parts:
            return 'others'
        return ','.join(acyl_parts)
    df['Acyl part'] = df.apply(determine_acyl_parts, axis=1)

    def determine_acyl_numb(row):
        acyl_numb = []
        if row['C2 loss'] == 'C2H4O2':
            acyl_numb.append('R1')

        if row['Acyl ion'] != 'none':
            acyl_ions = row['Acyl ion'].split(',')
            acyl_formulas = acyl_df['Acyl Formula'].tolist()
            r_parts = acyl_df['R parts'].tolist()

            for ion in acyl_ions:
                matches = [r_parts[i] for i, formula in enumerate(acyl_formulas) if formula.strip() == ion.strip()]
                if matches:
                    acyl_numb.extend(matches)

        if not acyl_numb:
            return 'others'
        return ','.join(acyl_numb)
    df['Acyl Number'] = df.apply(determine_acyl_numb, axis=1)

    def determine_acyl_simile(row):
        acyl_simile = []
        if row['C2 loss'] == 'C2H4O2':
            acyl_simile.append('CC(=O)*')

        if row['Acyl ion'] != 'none':
            acyl_ions = row['Acyl ion'].split(',')
            acyl_formulas = acyl_df['Acyl Formula'].tolist()
            acyl_SMILES = acyl_df['Acyl SMILES'].tolist()

            for ion in acyl_ions:
                matches = [acyl_SMILES[i] for i, formula in enumerate(acyl_formulas) if formula.strip() == ion.strip()]
                if matches:
                    acyl_simile.extend(matches)
        if not acyl_simile:
            return 'others'
        return ','.join(acyl_simile)

    df['Acyl SMILES'] = df.apply(determine_acyl_simile, axis=1)

    data_a = {
        'A number': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11'],
        'A Struc': [
            '1,2-en-3-one', '1,2-dihydro-3-one', '1,2-en-3-ol', '1,2-dihydro-3-ol', '1,10-en-3-one',
            '1-alkyl-3-one', '1-alkyl-2-ol-3-one', '1-alkyl-3-ol', '1-alkyl-2-ol-3-ol',
            '3,4-seco ring', 'bicyclo [2.2.1] heptane ring'
        ]
    }
    data_b = {
        'B number': ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8'],
        'B Struc': [
            '6,7-epoxy', '6,7-diol', '6,7-ene',
            '5-dehydro-6,7-epoxy', '5-dehydro-6,7-diol', '5-dehydro-6,7-ene',
            '4,7-epoxy', '4,6-epoxy'
        ]
    }
    data_c = {
        'C number': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
        'C Struc': [
            '9,13,14-orthoester', '9,13,14-orthoester', '9,13,14-triol', '9,13,14-triol',
            '9,13,14-orthoester and 18-ol', '9,13,14-orthoester and 18-ol'
        ],
        'C-12Struc': [
            'no C-12', 'with C-12', 'no C-12', 'with C-12', 'no C-12', 'with C-12'
        ]
    }
    data_m = {
        'M number': ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9'],
        'M Struc': [
            'C9 ring with no substituents', 'C10 ring with no substituents',
            'C10 ring with one substituents', 'C10 ring with two substituents',
            'C10 ring with three substituents', 'C10 ring with four substituents',
            'C14 ring with an olenfinic',
            'C16 ring with an olenfinic', 'C16 ring with one substituents'
        ]
    }

    df_a = pd.DataFrame(data_a)
    df_b = pd.DataFrame(data_b)
    df_c = pd.DataFrame(data_c)
    df_m = pd.DataFrame(data_m)

    def map_a_number(row):
        match = df_a[df_a['A Struc'] == row['A part']]
        if not match.empty:
            return match['A number'].values[0]
        return 'others'
    df['A number'] = df.apply(map_a_number, axis=1)

    def map_b_number(row):
        match = df_b[df_b['B Struc'] == row['B part']]
        if not match.empty:
            return match['B number'].values[0]
        return 'others'
    df['B number'] = df.apply(map_b_number, axis=1)

    def map_c_number(row):
        if row['type'] == 'D':
            match = df_c[(df_c['C Struc'] == row['C part']) & (df_c['C-12Struc'] == row['C-12'])]
            if not match.empty:
                return match['C number'].values[0]
        elif row['type'] == 'MD':
            match = df_c[df_c['C Struc'] == row['C part']]
            if not match.empty:
                return match['C number'].values[0]
        return 'others'
    df['C number'] = df.apply(map_c_number, axis=1)

    def map_m_number(row):
        match = df_m[df_m['M Struc'] == row['M part']]
        if not match.empty:
            return match['M number'].values[0]
        return 'others'
    df['M number'] = df.apply(map_m_number, axis=1)

    def find_skeleton_smiles(row, db):
        match = pd.DataFrame()
        if row['type'] == 'D':
            match = db[
                (db['A parts'] == row['A number']) &
                (db['B parts'] == row['B number']) &
                (db['C parts'] == row['C number'])
                ]
        elif row['type'] == 'MD':
            match = db[
                (db['A parts'] == row['A number']) &
                (db['B parts'] == row['B number']) &
                (db['C parts'] == row['C number']) &
                (db['M parts'] == row['M number'])
                ]
        if not match.empty:
            return match['New SMILES'].values[0]
        return 'others'

    df['A number'] = df.apply(map_a_number, axis=1)
    df['B number'] = df.apply(map_b_number, axis=1)
    df['C number'] = df.apply(map_c_number, axis=1)
    df['M number'] = df.apply(map_m_number, axis=1)

    df['Skeleton SMILES'] = df.apply(find_skeleton_smiles, axis=1, db=skeleton_db)

    insert_index = df.columns.get_loc('Acyl part') + 1
    for col in ['A number', 'B number', 'C number', 'M number', 'Acyl Number', 'Acyl SMILES', 'Skeleton SMILES']:
        cols = df.columns.tolist()
        cols.insert(insert_index, cols.pop(cols.index(col)))
        df = df[cols]
        insert_index += 1

    df.to_csv(output_path_csv, index=False)
    print(f"Annotation process complete. Data saved to {output_path_csv}")


""" Tab 8 related Functions """
@app.callback(
    Output('visualization-query-csv-file-path', 'children'),
    Input('visualization-query-upload-csv', 'contents'),
    State('visualization-query-upload-csv', 'filename')
)
def update_visualization_query_csv_file_path(csv_contents, csv_filenames):
    if csv_contents and csv_filenames:
        print(f"Received CSV filenames: {csv_filenames}")
        # Save the CSV files to the specified path
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
    State('visualization-query-upload-csv', 'filename')
)
def generate_scatter_plot(n_clicks, csv_filenames):
    if n_clicks == 0 or not csv_filenames:
        print("Scatter plot generation not triggered or no CSV file provided")
        return {}

    file_path = os.path.join(UPLOAD_DIR, 'tab_8_uploads', csv_filenames[0])
    print(f"Reading file for scatter plot: {file_path}")

    df = pd.read_csv(file_path)

    # Reshape the DataFrame to a long format for plotting
    long_df = df.melt(
        id_vars=[
            'index', 'row ID', 'row m/z', 'elemental composition',
            'adduct ion', 'row retention time',
            'Skeleton SMILES', 'Acyl SMILES', 'identification'
        ],
        value_vars=df.columns.difference([
            'index', 'row m/z', 'elemental composition', 'adduct ion',
            'row retention time', 'Skeleton SMILES',
            'Acyl SMILES', 'identification'
        ]),
        var_name='Plant Sample',
        value_name='Concentration'
    )

    long_df['Concentration'] = pd.to_numeric(long_df['Concentration'], errors='coerce')
    long_df = long_df[long_df['Concentration'] > 0]

    print(f"Data filtered for visualization: {len(long_df)} rows")

    # Define color mapping based on 'identification'
    color_map = {'possibly undescribed': 'red', 'possibly known': 'grey'}
    long_df['Color'] = long_df['identification'].apply(
        lambda x: 'possibly undescribed' if x == 'possibly undescribed' else 'possibly known'
    )

    fig = px.scatter_3d(
        long_df,
        x='index',
        y='Plant Sample',
        z='row retention time',
        size='Concentration',

        size_max=40,
        color='Color',
        color_discrete_map=color_map,
        hover_data={
            'row m/z': True,
            'elemental composition': True,
            'adduct ion': True,
            'row retention time': True
        }
    )
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title=dict(text="Peak Number",
                           font=dict(size=14))
            ),
            yaxis=dict(
                title=dict(text="\n\n\n\nPlant Sample",
                           font=dict(size=14))
            ),
            zaxis=dict(
                title=dict(text="Retention Time (min)",
                           font=dict(size=14))
            )
        )
    )

    return fig


@app.callback(
    [Output('skeleton-image', 'src'),
     Output('acyl-image-1', 'src'),
     Output('acyl-image-2', 'src')],
    [Input('scatter-plot', 'clickData')],
    State('visualization-query-upload-csv', 'filename')
)
def update_images(clickData, csv_filenames):
    if not clickData or not csv_filenames:
        print("No click data or no CSV filenames provided")
        return '', '', ''

    file_path = os.path.join(UPLOAD_DIR, 'tab_8_uploads', csv_filenames[0])
    print(f"Reading file for image update: {file_path}")

    df = pd.read_csv(file_path)

    point = clickData['points'][0]
    compound_number = point['x']
    plant_sample = str(point['y'])

    print(f"Selected compound_number: {compound_number}, plant_sample: {plant_sample}")
    print("Columns in df:", df.columns)

    matching_columns = [col for col in df.columns if plant_sample in col]
    if not matching_columns:
        print(f"Error: No matching column found for plant_sample value '{plant_sample}' in DataFrame columns.")
        return '', '', ''

    sample_column = matching_columns[0]
    print(f"Using column '{sample_column}' for plant sample '{plant_sample}'.")

    try:
        compound_number = float(compound_number)
    except ValueError:
        print(f"Warning: compound_number '{compound_number}' could not be converted to float.")

    filtered_df = df[(df['index'] == compound_number) & (df[sample_column] > 0)]
    if filtered_df.empty:
        print(f"No matching data found for compound_number: {compound_number} and plant_sample column: {sample_column}")
        return '', '', ''

    row = filtered_df.iloc[0]
    skeleton_smiles = row.get('Skeleton SMILES', '')
    acyl_smiles = row.get('Acyl SMILES', '')

    print(f"Skeleton SMILES: {skeleton_smiles}, Acyl SMILES: {acyl_smiles}")

    skeleton_img = smiles_to_image(skeleton_smiles) if skeleton_smiles and skeleton_smiles != 'others' else ''
    acyl_img_1, acyl_img_2 = '', ''
    if acyl_smiles and acyl_smiles != 'others':
        acyl_parts = acyl_smiles.split(',')
        acyl_img_1 = smiles_to_image(acyl_parts[0].strip()) if len(acyl_parts) > 0 else ''
        acyl_img_2 = smiles_to_image(acyl_parts[1].strip()) if len(acyl_parts) > 1 else ''

    return skeleton_img, acyl_img_1, acyl_img_2


def smiles_to_image(smiles, size=(200, 200)):
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
if __name__ == '__main__':
    app.run_server(debug=True)