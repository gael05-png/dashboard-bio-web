import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import matplotlib.pyplot as plt
from Bio.PDB import PDBList, PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from stmol import showmol
import py3Dmol
import os

# --- Configuración Visual Profesional ---
st.set_page_config(page_title="BioDashboard | C.G.L.S.", layout="wide", page_icon="🧬")

# CSS Personalizado: Estilo Minimalista y Transparente
st.markdown("""
<style>
    /* Hace que los contenedores sean transparentes y sutiles */
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        height: 45px;
        white-space: pre-wrap;
        background-color: transparent;
        border-radius: 5px;
        border: 1px solid rgba(150, 150, 150, 0.2);
        color: inherit;
    }
    .stTabs [aria-selected="true"] {
        background-color: rgba(100, 100, 100, 0.1);
        border-bottom: 2px solid #FF4B4B;
        font-weight: bold;
    }
    /* Estilo para tarjetas de métricas sin fondo blanco duro */
    div[data-testid="stMetric"] {
        background-color: rgba(128, 128, 128, 0.05);
        border: 1px solid rgba(128, 128, 128, 0.1);
        padding: 10px;
        border-radius: 8px;
    }
</style>
""", unsafe_allow_html=True)

# --- Encabezado ---
st.title("🧬 BioSuite: Análisis Estructural Avanzado")
st.markdown(f"**Dashboard Científico de Bioinformática** | Análisis en tiempo real de estructuras PDB.")
st.divider()

# --- Barra Lateral (Tu Firma) ---
with st.sidebar:
    st.header('🎛️ Panel de Control')
    pdb_id = st.text_input("Ingresa ID PDB:", "1A2C").upper()
    
    st.divider()
    st.subheader("Visualización")
    style_3d = st.selectbox("Estilo 3D", ["Cartoon", "Stick", "Sphere", "Line"])
    color_3d = st.selectbox("Color", ["spectrum", "secondary structure", "chain", "residue"])
    
    st.divider()
    # TU NOMBRE DESTACADO AQUÍ
    st.markdown("### 👨‍💻 Autor Principal")
    st.markdown("**Cristo Gael Lopezportillo Sánchez**")
    st.caption("Bioinformática & Desarrollo Web")
    st.info("Versión Final v3.5 Pro")

# --- Funciones (Backend) ---
def get_structure_data(id_pdb):
    if not os.path.exists('pdb'): os.makedirs('pdb')
    pdbl = PDBList(verbose=False)
    file = pdbl.retrieve_pdb_file(id_pdb, pdir='pdb', file_format='pdb')
    parser = PDBParser(QUIET=True)
    return parser.get_structure(id_pdb, file)

def calculate_contact_map(structure):
    model = structure[0]
    atoms = [res['CA'] for chain in model for res in chain if 'CA' in res]
    size = len(atoms)
    contact_map = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            contact_map[i, j] = atoms[i] - atoms[j]
    return contact_map

# Mapeo de colores simple para gráficos
aa_properties = {'A':'Hidrofóbico','V':'Hidrofóbico','L':'Hidrofóbico','I':'Hidrofóbico','M':'Hidrofóbico','F':'Hidrofóbico','W':'Hidrofóbico','P':'Hidrofóbico','G':'Polar','S':'Polar','T':'Polar','C':'Polar','Y':'Polar','N':'Polar','Q':'Polar','D':'Ácido','E':'Ácido','K':'Básico','R':'Básico','H':'Básico'}

# --- Lógica Principal ---
if pdb_id:
    try:
        structure = get_structure_data(pdb_id)
        # Limpieza de secuencia
        residues = [r for m in structure for c in m for r in c if r.get_id()[0]==' ']
        seq = "".join([seq1(r.get_resname()) for r in residues if seq1(r.get_resname()) != 'X'])
        analysed_seq = ProteinAnalysis(seq)
        
        # PESTAÑAS
        t1, t2, t3, t4 = st.tabs(["🧊 Visor 3D", "📊 Reporte Ejecutivo", "🔥 Mapa de Contactos", "⚗️ Simulación"])

        # 1. VISOR 3D
        with t1:
            col_viewer, col_info = st.columns([3, 1])
            with col_viewer:
                view = py3Dmol.view(query='pdb:'+pdb_id)
                if style_3d == "Cartoon": view.setStyle({'cartoon':{'color':color_3d}})
                elif style_3d == "Stick": view.setStyle({'stick':{'color':color_3d}})
                elif style_3d == "Sphere": view.setStyle({'sphere':{'color':color_3d}})
                elif style_3d == "Line": view.setStyle({'line':{'color':color_3d}})
                view.zoomTo()
                showmol(view, height=500, width=800)
            with col_info:
                st.markdown("### Datos Estructurales")
                st.
