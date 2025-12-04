import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
from Bio.PDB import PDBList, PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from stmol import showmol
import py3Dmol
import os

# --- Configuración Visual ---
st.set_page_config(page_title="BioSuite Ultimate", layout="wide", page_icon="🧬")

st.markdown("""
<style>
    .metric-card {
        background-color: #f0f2f6;
        border-left: 5px solid #4B4B4B;
        padding: 15px;
        border-radius: 5px;
        margin-bottom: 10px;
    }
    h3 { color: #2c3e50; }
</style>
""", unsafe_allow_html=True)

st.title("🧬 BioSuite: Análisis Estructural y Mutagénesis")
st.markdown("**Proyecto Final de Bioinformática** | Análisis de estructuras PDB en tiempo real.")
st.divider()

# --- Barra Lateral ---
with st.sidebar:
    st.header('🎛️ Configuración')
    pdb_id = st.text_input("ID PDB:", "1A2C").upper()
    
    st.divider()
    st.subheader("Visualización 3D")
    style_3d = st.selectbox("Estilo", ["Cartoon", "Stick", "Sphere", "Line"])
    color_3d = st.selectbox("Color", ["spectrum", "chain", "secondary structure", "residue"])
    
    st.divider()
    st.info("Desarrollado por el Equipo BioInfo 2025")

# --- Funciones ---
def get_structure_data(id_pdb):
    if not os.path.exists('pdb'): os.makedirs('pdb')
    pdbl = PDBList(verbose=False)
    file = pdbl.retrieve_pdb_file(id_pdb, pdir='pdb', file_format='pdb')
    parser = PDBParser(QUIET=True)
    return parser.get_structure(id_pdb, file)

def calculate_ramachandran(structure):
    phi_psi = []
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        for phi, psi in pp.get_phi_psi_list():
            if phi and psi:
                phi_psi.append([phi, psi])
    return np.array(phi_psi)

# Diccionario de propiedades químicas
aa_properties = {
    'A': 'Hidrofóbico', 'V': 'Hidrofóbico', 'L': 'Hidrofóbico', 'I': 'Hidrofóbico', 
    'M': 'Hidrofóbico', 'F': 'Hidrofóbico', 'W': 'Hidrofóbico', 'P': 'Hidrofóbico',
    'G': 'Polar', 'S': 'Polar', 'T': 'Polar', 'C': 'Polar', 'Y': 'Polar', 'N': 'Polar', 'Q': 'Polar',
    'D': 'Ácido (-)', 'E': 'Ácido (-)',
    'K': 'Básico (+)', 'R': 'Básico (+)', 'H': 'Básico (+)'
}

# --- Lógica Principal ---
if pdb_id:
    try:
        structure = get_structure_data(pdb_id)
        residues = [r for m in structure for c in m for r in c if r.get_id()[0]==' ']
        seq = "".join([seq1(r.get_resname()) for r in residues if seq1(r.get_resname()) != 'X'])
        analysed_seq = ProteinAnalysis(seq)
        
        # Tabs
        tab1, tab2, tab3, tab4 = st.tabs(["🧊 Visor 3D", "📐 Ramachandran", "⚗️ Mutaciones", "📊 Reporte Ejecutivo"])

        # TAB 1: 3D
        with tab1:
            c1, c2 = st.columns([3, 1])
            with c1:
                view = py3Dmol.view(query='pdb:'+pdb_id)
                if style_3d == "Cartoon": view.setStyle({'cartoon':{'color':color_3d}})
                elif style_3d == "Stick": view.setStyle({'stick':{'color':color_3d}})
                elif style_3d == "Sphere": view.setStyle({'sphere':{'color':color_3d}})
                elif style_3d == "Line": view.setStyle({'line':{'color':color_3d}})
                view.zoomTo()
                showmol(view, height=450, width=700)
            with c2:
                st.markdown(f"**Detalles:**")
                st.write(f"ID: `{pdb_id}`")
                st.write(f"Longitud: `{len(seq)} resid`")

        # TAB 2: Ramachandran
        with tab2:
            st.subheader("Validación Geométrica")
            angles = calculate_ramachandran(structure)
            if len(angles) > 0:
                angles = angles * 180 / np.pi
                df_rama = pd.DataFrame(angles, columns=['Phi', 'Psi'])
                chart = alt.Chart(df_rama).mark_circle(size=60, opacity=0.5).encode(
                    x=alt.X('Phi', scale=alt.Scale(domain=[-180, 180])),
                    y=alt.Y('Psi', scale=alt.Scale(domain=[-180, 180])),
                    color=alt.value('#6c5ce7'),
                    tooltip=['Phi', 'Psi']
                ).properties(height=400).interactive()
                st.altair_chart(chart, use_container_width=True)

        # TAB 3: Mutaciones
        with tab3:
            st.subheader("Simulador de Mutagénesis")
            c_input, c_result = st.columns(2)
            with c_input:
                pos = st.number_input("Posición", 1, len(seq), 1)
                new_aa = st.selectbox("Nuevo Aminoácido", list(aa_properties.keys()))
                orig = seq[pos-1]
                st.write(f"Mutación: **{orig}{pos}{new_aa}**")
            with c_result:
                mut_seq = seq[:pos-1] + new_aa + seq[pos:]
                p_wt = ProteinAnalysis(seq)
                p_mut = ProteinAnalysis(mut_seq)
                d_mw = p_mut.molecular_weight() - p_wt.molecular_weight()
                st.metric("Cambio Peso Molecular", f"{d_mw:.2f} Da", delta=d_mw)

        # TAB 4: REPORTE MEJORADO (Aquí está el cambio grande)
        with tab4:
            st.subheader(f"Informe Fisicoquímico: {pdb_id}")
            
            # 1. Tarjetas de Métricas (Fila Superior)
            col_m1, col_m2, col_m3, col_m4 = st.columns(4)
            
            mw = analysed_seq.molecular_weight()
            pi = analysed_seq.isoelectric_point()
            instability = analysed_seq.instability_index()
            aromatic = analysed_seq.aromaticity()
            
            col_m1.metric("Peso Molecular", f"{mw/1000:.1f} kDa")
            col_m2.metric("Punto Isoeléctrico", f"{pi:.2f} pH")
            col_m3.metric("Aromaticidad", f"{aromatic*100:.1f}%")
            
            # Lógica de color para inestabilidad
            if instability > 40:
                col_m4.metric("Estabilidad", "Inestable", f"{instability:.2f}", delta_color="inverse")
            else:
                col_m4.metric("Estabilidad", "Estable", f"{instability:.2f}")

            st.divider()

            # 2. Gráficos Organizados (Dos columnas)
            g_col1, g_col2 = st.columns(2)

            # Preparar datos
            aa_counts = {k: seq.count(k) for k in set(seq)}
            df_aa = pd.DataFrame(list(aa_counts.items()), columns=['AA', 'Count'])
            
            # Agregar grupo químico al DataFrame
            df_aa['Grupo'] = df_aa['AA'].map(aa_properties)
            
            with g_col1:
                st.markdown("#### Distribución de Aminoácidos")
                # Gráfico de Barras ordenado
                bar_chart = alt.Chart(df_aa).mark_bar().encode(
                    x=alt.X('AA', sort='-y', title='Aminoácido'),
                    y=alt.Y('Count', title='Frecuencia'),
                    color='Grupo',
                    tooltip=['AA', 'Count', 'Grupo']
                ).interactive()
                st.altair_chart(bar_chart, use_container_width=True)

            with g_col2:
                st.markdown("#### Composición Química")
                # Agrupar por propiedades
                df_props = df_aa.groupby('Grupo').sum(numeric_only=True).reset_index()
                
                # Gráfico de Donas (Donut Chart)
                base = alt.Chart(df_props).encode(theta=alt.Theta("Count", stack=True))
                pie = base.mark_arc(outerRadius=120, innerRadius=80).encode(
                    color=alt.Color("Grupo"),
                    order=alt.Order("Count", sort="descending"),
                    tooltip=["Grupo", "Count"]
                )
                text = base.mark_text(radius=140).encode(
                    text="Count",
                    order=alt.Order("Count", sort="descending"),
                    color=alt.value("black") 
                )
                st.altair_chart(pie + text, use_container_width=True)
            
            # 3. Expander para datos crudos
            with st.expander("Ver secuencia completa FASTA"):
                st.code(f">{pdb_id}\n{seq}")

    except Exception as e:
        st.error(f"Error: {e}")
