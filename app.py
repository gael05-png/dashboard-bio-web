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

# --- Configuración Visual ---
st.set_page_config(page_title="BioSuite Ultimate", layout="wide", page_icon="🧬")

st.markdown("""
<style>
    .metric-card { background-color: #f0f2f6; padding: 15px; border-radius: 10px; }
    .stTabs [data-baseweb="tab-list"] { gap: 10px; }
    .stTabs [data-baseweb="tab"] { height: 50px; white-space: pre-wrap; background-color: #f0f2f6; border-radius: 5px; }
    .stTabs [aria-selected="true"] { background-color: #e6f2ff; border-bottom: 2px solid #4B4B4B; }
</style>
""", unsafe_allow_html=True)

st.title("🧬 BioSuite: Análisis Avanzado de Proteínas")
st.markdown("**Dashboard Científico** | Estructura, Mutagénesis y Mapas de Contacto.")
st.divider()

# --- Barra Lateral ---
with st.sidebar:
    st.header('🎛️ Control')
    pdb_id = st.text_input("ID PDB:", "1A2C").upper()
    
    st.divider()
    st.subheader("Visualización 3D")
    style_3d = st.selectbox("Estilo", ["Cartoon", "Stick", "Sphere", "Line"])
    color_3d = st.selectbox("Color", ["spectrum", "secondary structure", "chain", "residue"])
    
    st.divider()
    st.download_button("📥 Descargar Manual PDF", "Texto de ejemplo", file_name="manual.pdf") # Simulado
    st.info("v3.0. Final Build")

# --- Funciones Científicas ---
def get_structure_data(id_pdb):
    if not os.path.exists('pdb'): os.makedirs('pdb')
    pdbl = PDBList(verbose=False)
    file = pdbl.retrieve_pdb_file(id_pdb, pdir='pdb', file_format='pdb')
    parser = PDBParser(QUIET=True)
    return parser.get_structure(id_pdb, file)

def calculate_contact_map(structure, chain_id=None):
    # Obtener átomos CA (Carbono Alfa)
    model = structure[0]
    atoms = []
    for chain in model:
        if chain_id and chain.id != chain_id: continue
        for res in chain:
            if 'CA' in res:
                atoms.append(res['CA'])
            elif 'C' in res: # Fallback si no hay CA
                atoms.append(res['C'])
                
    # Calcular matriz de distancias
    size = len(atoms)
    contact_map = np.zeros((size, size))
    
    for i in range(size):
        for j in range(size):
            dist = atoms[i] - atoms[j]
            contact_map[i, j] = dist
            
    return contact_map

# Diccionario propiedades
aa_properties = {
    'A': 'Hidrofóbico', 'V': 'Hidrofóbico', 'L': 'Hidrofóbico', 'I': 'Hidrofóbico', 
    'M': 'Hidrofóbico', 'F': 'Hidrofóbico', 'W': 'Hidrofóbico', 'P': 'Hidrofóbico',
    'G': 'Polar', 'S': 'Polar', 'T': 'Polar', 'C': 'Polar', 'Y': 'Polar', 'N': 'Polar', 'Q': 'Polar',
    'D': 'Ácido (-)', 'E': 'Ácido (-)', 'K': 'Básico (+)', 'R': 'Básico (+)', 'H': 'Básico (+)'
}

# --- Lógica Principal ---
if pdb_id:
    try:
        structure = get_structure_data(pdb_id)
        residues = [r for m in structure for c in m for r in c if r.get_id()[0]==' ']
        seq = "".join([seq1(r.get_resname()) for r in residues if seq1(r.get_resname()) != 'X'])
        analysed_seq = ProteinAnalysis(seq)
        
        # Tabs
        t1, t2, t3, t4, t5 = st.tabs(["🧊 3D", "📊 Reporte", "📐 Ramachandran", "🔥 Mapa de Contacto", "⚗️ Mutaciones"])

        # T1: 3D
        with t1:
            c1, c2 = st.columns([3, 1])
            with c1:
                view = py3Dmol.view(query='pdb:'+pdb_id)
                if style_3d == "Cartoon": view.setStyle({'cartoon':{'color':color_3d}})
                elif style_3d == "Stick": view.setStyle({'stick':{'color':color_3d}})
                elif style_3d == "Sphere": view.setStyle({'sphere':{'color':color_3d}})
                elif style_3d == "Line": view.setStyle({'line':{'color':color_3d}})
                view.zoomTo()
                showmol(view, height=500, width=800)
            with c2:
                st.markdown("### Detalles")
                st.write(f"**ID:** {pdb_id}")
                st.write(f"**Longitud:** {len(seq)} aminoácidos")
                
                # Buscador de motivos
                st.markdown("---")
                search_motif = st.text_input("Buscar Motivo (ej: AAA)", "").upper()
                if search_motif:
                    count = seq.count(search_motif)
                    st.write(f"Encontrado {count} veces")
                    if count > 0: st.success("¡Motivo presente!")

        # T2: REPORTE EJECUTIVO
        with t2:
            st.subheader("Informe Fisicoquímico")
            # Métricas
            m1, m2, m3, m4 = st.columns(4)
            mw = analysed_seq.molecular_weight()
            instability = analysed_seq.instability_index()
            
            m1.metric("Peso Molecular", f"{mw/1000:.1f} kDa")
            m2.metric("Punto Isoeléctrico", f"{analysed_seq.isoelectric_point():.2f} pH")
            m3.metric("Aromaticidad", f"{analysed_seq.aromaticity()*100:.1f}%")
            m4.metric("Inestabilidad", f"{instability:.2f}", delta="Estable" if instability < 40 else "Inestable", delta_color="inverse")

            st.divider()
            
            # Gráficos
            g1, g2 = st.columns(2)
            
            # Dataframe para graficos
            aa_counts = {k: seq.count(k) for k in set(seq)}
            df_aa = pd.DataFrame(list(aa_counts.items()), columns=['AA', 'Count'])
            df_aa['Grupo'] = df_aa['AA'].map(aa_properties)
            
            with g1:
                st.markdown("**Distribución por Tipo**")
                base = alt.Chart(df_aa).encode(theta=alt.Theta("Count", stack=True))
                pie = base.mark_arc(outerRadius=100).encode(
                    color=alt.Color("Grupo"),
                    order=alt.Order("Count", sort="descending"),
                    tooltip=["Grupo", "Count"]
                )
                st.altair_chart(pie, use_container_width=True)
                
            with g2:
                st.markdown("**Descargar Datos**")
                st.write("Descarga los datos crudos del análisis para tu reporte.")
                csv = df_aa.to_csv(index=False).encode('utf-8')
                st.download_button(
                    "📄 Descargar CSV de Aminoácidos",
                    csv,
                    "analisis_aa.csv",
                    "text/csv",
                    key='download-csv'
                )

        # T3: RAMACHANDRAN (Geometría)
        with t3:
            st.subheader("Gráfico de Ramachandran")
            st.write("Validación de ángulos de torsión del esqueleto proteico.")
            # Codigo simplificado para velocidad
            ppb = PPBuilder()
            phi_psi = []
            for pp in ppb.build_peptides(structure):
                for phi, psi in pp.get_phi_psi_list():
                    if phi and psi: phi_psi.append([phi*180/np.pi, psi*180/np.pi])
            
            if len(phi_psi) > 0:
                df_rama = pd.DataFrame(phi_psi, columns=['Phi', 'Psi'])
                chart = alt.Chart(df_rama).mark_circle(size=60).encode(
                    x=alt.X('Phi', scale=alt.Scale(domain=[-180, 180])),
                    y=alt.Y('Psi', scale=alt.Scale(domain=[-180, 180])),
                    color=alt.value('teal'),
                    tooltip=['Phi', 'Psi']
                ).interactive()
                st.altair_chart(chart, use_container_width=True)

        # T4: MAPA DE CONTACTOS (NUEVO)
        with t4:
            st.subheader("Mapa de Contactos (Heatmap)")
            st.markdown("Representación matricial de las distancias espaciales entre residuos. Las regiones oscuras indican cercanía (folding).")
            
            with st.spinner("Calculando distancias 3D..."):
                dist_matrix = calculate_contact_map(structure)
                
                # Plot con Matplotlib
                fig, ax = plt.subplots()
                im = ax.imshow(dist_matrix, cmap='viridis_r', origin='lower')
                plt.colorbar(im, label="Distancia (Å)")
                ax.set_title(f"Matriz de Distancias: {pdb_id}")
                ax.set_xlabel("Número de Residuo")
                ax.set_ylabel("Número de Residuo")
                st.pyplot(fig)
                
            st.caption("Eje X e Y representan la secuencia de aminoácidos. La diagonal representa la distancia 0 (residuo consigo mismo).")

        # T5: MUTACIONES
        with t5:
            st.subheader("Simulador")
            c1, c2 = st.columns(2)
            with c1:
                pos = st.number_input("Posición", 1, len(seq), 1)
                new = st.selectbox("Cambiar a:", list(aa_properties.keys()))
                st.write(f"Mutar: {seq[pos-1]} ➝ {new}")
            with c2:
                p_wt = ProteinAnalysis(seq)
                p_mut = ProteinAnalysis(seq[:pos-1] + new + seq[pos:])
                d_mw = p_mut.molecular_weight() - p_wt.molecular_weight()
                st.metric("Δ Peso Molecular", f"{d_mw:.2f}", delta_color="off")

    except Exception as e:
        st.error(f"Error cargando estructura: {e}")
