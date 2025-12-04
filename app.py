import streamlit as st
import pandas as pd
import altair as alt
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from stmol import showmol
import py3Dmol
import os

st.set_page_config(page_title="BioDashboard Web", layout="wide")
st.title("🧬 BioDashboard: Análisis de Proteínas")
st.markdown("Entrega de Proyecto Final - Bioinformática")

# Sidebar
pdb_id = st.sidebar.text_input("ID PDB:", "6LU7").upper()
st.sidebar.write("Equipo: [Cristo Gael Lopezportillo Sánchez]")

# Funciones
def get_seq(structure):
    residues = [r.get_resname() for m in structure for c in m for r in c if r.get_id()[0]==' '][:200] # Limitado para rapidez
    seq = ""
    for r in residues:
        try: seq += seq1(r)
        except: pass
    return seq

if pdb_id:
    try:
        if not os.path.exists('pdb'): os.makedirs('pdb')
        pdbl = PDBList(verbose=False)
        file = pdbl.retrieve_pdb_file(pdb_id, pdir='pdb', file_format='pdb')
        
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(pdb_id, file)
        seq = get_seq(struct)
        
        # Tabs
        t1, t2 = st.tabs(["🧊 3D", "📊 Datos"])
        
        with t1:
            view = py3Dmol.view(query='pdb:'+pdb_id)
            view.setStyle({'cartoon':{'color':'spectrum'}})
            view.zoomTo()
            showmol(view, height=500, width=800)
            
        with t2:
            st.write("Secuencia detectada:", seq)
            # Grafico simple
            df = pd.DataFrame({'Amino': list(seq)}).value_counts().reset_index()
            df.columns = ['Amino', 'Count']
            st.altair_chart(alt.Chart(df).mark_bar().encode(x='Amino', y='Count'), use_container_width=True)
            
            # Quimica
            pa = ProteinAnalysis(seq)
            st.metric("Peso Molecular", f"{pa.molecular_weight():.1f}")
            
    except Exception as e:
        st.error(f"Error o ID inválido: {e}")
