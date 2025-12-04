import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import matplotlib.pyplot as plt
import requests
from deep_translator import GoogleTranslator
from Bio.PDB import PDBList, PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Align import PairwiseAligner
from stmol import showmol
import py3Dmol
import os
from datetime import date

# --- CONFIGURACIÓN VISUAL ---
st.set_page_config(page_title="BioDashboard Ultimate | C.G.L.S.", layout="wide", page_icon="🧬")

st.markdown("""
<style>
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: rgba(255, 255, 255, 0.05);
        border-radius: 8px;
        border: 1px solid rgba(150, 150, 150, 0.2);
        color: inherit;
    }
    .stTabs [aria-selected="true"] {
        background-color: rgba(100, 100, 100, 0.1);
        border-bottom: 2px solid #FF4B4B;
        font-weight: bold;
    }
    div[data-testid="stMetric"] {
        background-color: rgba(128, 128, 128, 0.05);
        border: 1px solid rgba(128, 128, 128, 0.1);
        padding: 15px;
        border-radius: 10px;
    }
    .definition-box {
        padding: 15px;
        border-left: 5px solid #FF4B4B;
        background-color: rgba(128, 128, 128, 0.1);
        margin-bottom: 20px;
        border-radius: 5px;
    }
</style>
""", unsafe_allow_html=True)

# --- ENCABEZADO ---
st.title("🧬 BioSuite X: Análisis Bioinformático Integral")
st.markdown("**Plataforma de visualización molecular, análisis fisicoquímico y alineamiento de secuencias.**")
st.divider()

# --- BARRA LATERAL ---
with st.sidebar:
    st.header('🎛️ Panel de Control')
    pdb_id = st.text_input("Buscar ID PDB:", "6LU7").upper()
    
    st.divider()
    st.subheader("Visualización 3D")
    style_3d = st.selectbox("Estilo", ["Cartoon", "Stick", "Sphere", "Line"])
    color_3d = st.selectbox("Color", ["spectrum", "chain", "secondary structure", "residue"])
    surface = st.checkbox("Mostrar Superficie Volumétrica", value=False)
    
    st.divider()
    st.markdown("### Desarrollado por:")
    st.markdown("**Cristo Gael Lopezportillo Sánchez**")
    st.caption("Proyecto Final de Bioinformática | 2025")
    
    # Enlaces externos
    if pdb_id:
        st.divider()
        st.markdown("**Enlaces Externos:**")
        st.link_button("Ver en RCSB PDB", f"https://www.rcsb.org/structure/{pdb_id}")
        st.link_button("Ver en UniProt", f"https://www.uniprot.org/")

# --- FUNCIONES ---
def get_pdb_info(pdb_id):
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            title_en = data['struct']['title']
            class_en = data['struct_keywords']['pdbx_keywords']
            translator = GoogleTranslator(source='auto', target='es')
            return translator.translate(title_en), translator.translate(class_en)
    except:
        return None, None
    return "Descripción no disponible", "Desconocido"

def get_structure(id_pdb):
    if not os.path.exists('pdb'): os.makedirs('pdb')
    pdbl = PDBList(verbose=False)
    file = pdbl.retrieve_pdb_file(id_pdb, pdir='pdb', file_format='pdb')
    parser = PDBParser(QUIET=True)
    return parser.get_structure(id_pdb, file)

def get_sequence(structure):
    residues = [r for m in structure for c in m for r in c if r.get_id()[0]==' ']
    return "".join([seq1(r.get_resname()) for r in residues if seq1(r.get_resname()) != 'X'])

def calculate_contact_map(structure):
    model = structure[0]
    atoms = [res['CA'] for chain in model for res in chain if 'CA' in res]
    size = len(atoms)
    contact_map = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            contact_map[i, j] = np.linalg.norm(atoms[i].coord - atoms[j].coord)
    return contact_map

# Datos fijos para evitar errores
aa_props = {'A':'Hidrofóbico','V':'Hidrofóbico','L':'Hidrofóbico','I':'Hidrofóbico','M':'Hidrofóbico','F':'Hidrofóbico','W':'Hidrofóbico','P':'Hidrofóbico','G':'Polar','S':'Polar','T':'Polar','C':'Polar','Y':'Polar','N':'Polar','Q':'Polar','D':'Ácido','E':'Ácido','K':'Básico','R':'Básico','H':'Básico'}
kd_scale = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}

# --- LÓGICA PRINCIPAL ---
if pdb_id:
    try:
        desc_title, desc_class = get_pdb_info(pdb_id)
        
        if desc_title:
            st.markdown(f"""
            <div class="definition-box">
                <h3 style="margin:0; color:#FF4B4B;">{pdb_id} | {desc_class.upper() if desc_class else 'SIN CLASE'}</h3>
                <p style="font-size:18px; margin-top:5px;"><strong>Definición:</strong> {desc_title.capitalize()}</p>
            </div>
            """, unsafe_allow_html=True)
        
        struct1 = get_structure(pdb_id)
        seq1_str = get_sequence(struct1)
        analysed_seq = ProteinAnalysis(seq1_str)
        
        # 6 PESTAÑAS AHORA
        t1, t2, t3, t4, t5, t6 = st.tabs(["🧊 Visor 3D", "📊 Reporte", "📈 Radar", "📉 Hidrofobicidad", "⚔️ Comparador", "🔥 Heatmaps"])

        # TAB 1: 3D
        with t1:
            c1, c2 = st.columns([3, 1])
            with c1:
                view = py3Dmol.view(query='pdb:'+pdb_id)
                if style_3d == "Cartoon": view.setStyle({'cartoon':{'color':color_3d}})
                elif style_3d == "Stick": view.setStyle({'stick':{'color':color_3d}})
                elif style_3d == "Sphere": view.setStyle({'sphere':{'color':color_3d}})
                elif style_3d == "Line": view.setStyle({'line':{'color':color_3d}})
                if surface: view.addSurface(py3Dmol.VDW, {'opacity':0.5, 'color':'white'})
                view.zoomTo()
                showmol(view, height=500, width=800)
            with c2:
                st.markdown("### Detalles")
                st.write(f"**Longitud:** `{len(seq1_str)}` residuos")
                header = struct1.header
                if 'helix' in header: st.write(f"**Hélices:** {len(header['helix'])}")
                if 'sheet' in header: st.write(f"**Láminas:** {len(header['sheet'])}")
                st.info("Renderizado WebGL Activo")

        # TAB 2: REPORTE
        with t2:
            st.subheader("Informe Fisicoquímico")
            c1, c2, c3, c4 = st.columns(4)
            mw = analysed_seq.molecular_weight()
            inst = analysed_seq.instability_index()
            c1.metric("Peso Molecular", f"{mw/1000:.1f} kDa")
            c2.metric("Punto Isoeléctrico", f"{analysed_seq.isoelectric_point():.2f} pH")
            c3.metric("Hidropatía", f"{analysed_seq.gravy():.2f}")
            c4.metric("Estabilidad", f"{inst:.2f}", delta="Inestable" if inst>40 else "Estable", delta_color="inverse")
            
            st.divider()
            
            df_aa = pd.DataFrame(list({k: seq1_str.count(k) for k in set(seq1_str)}.items()), columns=['AA', 'Count'])
            df_aa['Tipo'] = df_aa['AA'].map(aa_props)
            
            gc1, gc2 = st.columns(2)
            with gc1:
                 st.markdown("**Composición de Aminoácidos**")
                 st.altair_chart(alt.Chart(df_aa).mark_bar().encode(x='AA', y='Count', color='Tipo').interactive(), use_container_width=True)
            with gc2:
                 st.markdown("**Generador de Cita (APA)**")
                 today = date.today()
                 citation = f"RCSB Protein Data Bank. (2025). Structure {pdb_id}. Recuperado el {today} de rcsb.org. Análisis realizado mediante BioSuite X por C.G. Lopezportillo."
                 st.code(citation, language='text')

        # TAB 3: RADAR (NUEVO)
        with t3:
            st.subheader("Comparativa Normalizada (Radar)")
            st.write("Comparación de la proteína actual (Azul) vs. Promedio Típico (Rojo).")
            
            # Datos simulados del promedio (benchmark)
            categories = ['Hidrofobicidad', 'Estabilidad', 'Aromaticidad', 'Acidez (pH)']
            
            # Normalización simple para el gráfico (0 a 1)
            # Hidro: rango aprox -2 a 2 -> norm 0.5 base
            val_hydro = (analysed_seq.gravy() + 2) / 4
            # Estabilidad: rango 0 a 100 -> norm
            val_stab = inst / 100
            # Aromaticidad: 0 a 0.2
            val_arom = analysed_seq.aromaticity() * 5
            # pH: 0 a 14 -> norm
            val_ph = analysed_seq.isoelectric_point() / 14
            
            values = [val_hydro, val_stab, val_arom, val_ph]
            values_avg = [0.5, 0.4, 0.4, 0.5] # Valores teóricos promedio
            
            # Matplotlib Radar Chart
            label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(values))
            
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(polar=True)
            fig.patch.set_alpha(0)
            ax.patch.set_alpha(0)
            
            # Plot Data
            ax.plot(label_loc, values, label=f'{pdb_id}', color='cyan')
            ax.fill(label_loc, values, color='cyan', alpha=0.25)
            
            # Plot Avg
            ax.plot(label_loc, values_avg, label='Promedio Global', color='red', linestyle='--')
            
            ax.set_xticks(label_loc)
            ax.set_xticklabels(categories)
            ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
            
            # Estilos oscuros para que se vea bien
            ax.tick_params(colors='gray')
            ax.spines['polar'].set_color('gray')
            
            st.pyplot(fig)
            st.caption("Nota: Los valores están normalizados (0-1) para facilitar la comparación visual.")

        # TAB 4: HIDROFOBICIDAD
        with t4:
            st.subheader("Perfil de Kyte & Doolittle")
            values = [kd_scale.get(aa, 0) for aa in seq1_str]
            window = 9
            weights = np.repeat(1.0, window)/window
            sma = np.convolve(values, weights, 'valid')
            df_kd = pd.DataFrame({'Posición': range(1, len(sma)+1), 'Hidrofobicidad': sma})
            st.altair_chart(alt.Chart(df_kd).mark_line().encode(x='Posición', y='Hidrofobicidad', color=alt.value('#FF4B4B')).properties(height=400).interactive(), use_container_width=True)

        # TAB 5: COMPARADOR
        with t5:
            st.subheader("Alineamiento por Pares")
            id2 = st.text_input("Comparar contra ID (ej: 1CRN):", "").upper()
            if id2:
                try:
                    t2_es, _ = get_pdb_info(id2)
                    st.write(f"VS: **{t2_es}**")
                    struct2 = get_structure(id2)
                    seq2_str = get_sequence(struct2)
                    aligner = PairwiseAligner()
                    aligner.mode = 'global'
                    score = aligner.score(seq1_str, seq2_str)
                    identity = (score / max(len(seq1_str), len(seq2_str))) * 100
                    st.metric("Identidad Genética", f"{identity:.2f}%")
                except: st.error("Error al cargar ID 2.")

        # TAB 6: HEATMAPS
        with t6:
            st.subheader("Mapa de Contactos")
            with st.spinner("Procesando..."):
                matrix = calculate_contact_map(struct1)
                fig, ax = plt.subplots()
                fig.patch.set_alpha(0)
                ax.patch.set_alpha(0)
                im = ax.imshow(matrix, cmap='viridis', origin='lower')
                plt.colorbar(im, label="Å")
                ax.tick_params(colors='gray')
                st.pyplot(fig)

    except Exception as e:
        st.error(f"Error: {e}")
else:
    st.info("Ingresa un ID para iniciar.")
