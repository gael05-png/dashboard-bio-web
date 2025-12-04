import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import matplotlib.pyplot as plt
import requests
from Bio.PDB import PDBList, PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Align import PairwiseAligner
from stmol import showmol
import py3Dmol
import os

# --- CONFIGURACIÓN VISUAL ---
st.set_page_config(page_title="BioDashboard | C.G.L.S.", layout="wide", page_icon="🧬")

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
st.title("🧬 BioSuite X: Plataforma Integral")
st.markdown("**Análisis Estructural, Fisicoquímico y Comparativo de Proteínas**")
st.divider()

# --- BARRA LATERAL ---
with st.sidebar:
    st.header('🎛️ Control Principal')
    pdb_id = st.text_input("Buscar ID PDB:", "6LU7").upper()
    
    st.divider()
    st.subheader("Configuración 3D")
    style_3d = st.selectbox("Estilo", ["Cartoon", "Stick", "Sphere", "Line"])
    color_3d = st.selectbox("Color", ["spectrum", "chain", "secondary structure", "residue"])
    surface = st.checkbox("Mostrar Superficie", value=False)
    
    st.divider()
    st.markdown("### 👨‍🔬 Lead Developer")
    st.markdown("**Cristo Gael Lopezportillo Sánchez**")
    st.caption("Bioinformática Avanzada | 2025")

# --- FUNCIONES ---
def get_pdb_info(pdb_id):
    # Conexión a la API del PDB para sacar la definición
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            title = data['struct']['title']
            classification = data['struct_keywords']['pdbx_keywords']
            return title, classification
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

aa_props = {'A':'Hidrofóbico','V':'Hidrofóbico','L':'Hidrofóbico','I':'Hidrofóbico','M':'Hidrofóbico','F':'Hidrofóbico','W':'Hidrofóbico','P':'Hidrofóbico','G':'Polar','S':'Polar','T':'Polar','C':'Polar','Y':'Polar','N':'Polar','Q':'Polar','D':'Ácido','E':'Ácido','K':'Básico','R':'Básico','H':'Básico'}

# --- LÓGICA PRINCIPAL ---
if pdb_id:
    try:
        # 1. OBTENER DEFINICIÓN (NUEVO)
        desc_title, desc_class = get_pdb_info(pdb_id)
        
        if desc_title:
            st.markdown(f"""
            <div class="definition-box">
                <h3 style="margin:0; color:#FF4B4B;">{pdb_id} | {desc_class}</h3>
                <p style="font-size:18px; margin-top:5px;"><strong>Definición:</strong> {desc_title.capitalize()}</p>
            </div>
            """, unsafe_allow_html=True)
        
        # 2. CARGAR DATOS
        struct1 = get_structure(pdb_id)
        seq1_str = get_sequence(struct1)
        analysed_seq = ProteinAnalysis(seq1_str)
        
        # PESTAÑAS
        t1, t2, t3, t4, t5 = st.tabs(["🧊 Visor 3D", "📊 Reporte", "⚔️ Comparador", "🔥 Heatmaps", "⚗️ Laboratorio"])

        # VISOR 3D
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
                st.markdown("### Detalles Técnicos")
                st.write(f"**Longitud:** `{len(seq1_str)}` residuos")
                st.write(f"**Clasificación:** `{desc_class}`")
                st.info("Renderizado WebGL Activo")

        # REPORTE
        with t2:
            st.subheader(f"Informe Fisicoquímico: {pdb_id}")
            col1, col2, col3, col4 = st.columns(4)
            mw = analysed_seq.molecular_weight()
            inst = analysed_seq.instability_index()
            col1.metric("Peso Molecular", f"{mw/1000:.1f} kDa")
            col2.metric("Punto Isoeléctrico", f"{analysed_seq.isoelectric_point():.2f} pH")
            col3.metric("Hidropatía", f"{analysed_seq.gravy():.2f}")
            col4.metric("Estabilidad", f"{inst:.2f}", delta="Inestable" if inst>40 else "Estable", delta_color="inverse")
            
            st.divider()
            st.markdown(f"**Conclusión Automática:** La estructura identificada como *{desc_title}* presenta un perfil {'ácido' if analysed_seq.isoelectric_point() < 7 else 'básico'} y se clasifica como {desc_class}.")
            
            df_aa = pd.DataFrame(list({k: seq1_str.count(k) for k in set(seq1_str)}.items()), columns=['AA', 'Count'])
            df_aa['Tipo'] = df_aa['AA'].map(aa_props)
            base = alt.Chart(df_aa).encode(theta=alt.Theta("Count", stack=True))
            pie = base.mark_arc(innerRadius=70).encode(color=alt.Color("Tipo"), order=alt.Order("Count", sort="descending"), tooltip=["Tipo", "Count"])
            st.altair_chart(pie, use_container_width=True)
            st.caption(f"Autor: **Cristo Gael Lopezportillo Sánchez**")

        # COMPARADOR
        with t3:
            st.subheader("Alineamiento por Pares")
            id2 = st.text_input("Comparar contra ID:", "").upper()
            if id2:
                try:
                    title2, class2 = get_pdb_info(id2) # Obtener nombre de la segunda también
                    st.write(f"Comparando con: **{title2}**")
                    struct2 = get_structure(id2)
                    seq2_str = get_sequence(struct2)
                    aligner = PairwiseAligner()
                    aligner.mode = 'global'
                    score = aligner.score(seq1_str, seq2_str)
                    identity = (score / max(len(seq1_str), len(seq2_str))) * 100
                    st.metric("Identidad Genética", f"{identity:.2f}%")
                    st.code(f"A ({pdb_id}): {seq1_str[:50]}...\nB ({id2}): {seq2_str[:50]}...")
                except: st.error("Error al cargar la segunda proteína.")

        # HEATMAPS
        with t4:
            st.subheader("Mapa de Contactos")
            with st.spinner("Procesando..."):
                matrix = calculate_contact_map(struct1)
                fig, ax = plt.subplots()
                fig.patch.set_alpha(0)
                ax.patch.set_alpha(0)
                im = ax.imshow(matrix, cmap='plasma', origin='lower')
                plt.colorbar(im, label="Å")
                ax.tick_params(colors='gray')
                st.pyplot(fig)

        # SIMULACIÓN
        with t5:
            st.subheader("Mutagénesis")
            c1, c2 = st.columns(2)
            with c1:
                pos = st.number_input("Residuo", 1, len(seq1_str), 1)
                new = st.selectbox("Cambio", list(aa_props.keys()))
                st.write(f"**{seq1_str[pos-1]}{pos} ➝ {new}**")
            with c2:
                p_wt = ProteinAnalysis(seq1_str)
                p_mut = ProteinAnalysis(seq1_str[:pos-1] + new + seq1_str[pos:])
                d_mw = p_mut.molecular_weight() - p_wt.molecular_weight()
                st.metric("Δ Peso", f"{d_mw:.2f} Da", delta=d_mw)

    except Exception as e:
        st.error(f"Error: {e}")

else:
    st.info("Ingresa un ID para iniciar el sistema.")
