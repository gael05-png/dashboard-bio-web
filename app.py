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

# --- Configuración Avanzada ---
st.set_page_config(page_title="BioSuite Ultimate", layout="wide", page_icon="🧬")

# Estilos CSS para que se vea profesional
st.markdown("""
<style>
    .big-font { font-size:20px !important; font-weight: bold; }
    .metric-container { background-color: #f0f2f6; padding: 10px; border-radius: 10px; }
</style>
""", unsafe_allow_html=True)

st.title("🧬 BioSuite: Análisis Estructural y Mutagénesis")
st.markdown("Plataforma avanzada para visualización molecular, validación geométrica y simulación de mutaciones.")
st.divider()

# --- Barra Lateral ---
with st.sidebar:
    st.header('🎛️ Centro de Control')
    pdb_id = st.text_input("ID PDB:", "1A2C").upper()
    
    st.divider()
    st.subheader("⚙️ Configuración 3D")
    style_3d = st.selectbox("Estilo de Visualización", ["Cartoon", "Stick", "Sphere", "Line"])
    color_3d = st.selectbox("Esquema de Color", ["spectrum", "chain", "residue", "secondary structure"])
    
    st.divider()
    st.write("Desarrollado por: *Cristo Gael Lopezportillo Sánchez*")
    st.write("Equipo de Bioinformática 2024")

# --- Funciones Científicas ---
def get_structure_data(id_pdb):
    if not os.path.exists('pdb'): os.makedirs('pdb')
    pdbl = PDBList(verbose=False)
    file = pdbl.retrieve_pdb_file(id_pdb, pdir='pdb', file_format='pdb')
    parser = PDBParser(QUIET=True)
    return parser.get_structure(id_pdb, file)

def calculate_ramachandran(structure):
    # Calcula ángulos Phi y Psi para validación geométrica
    phi_psi = []
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        for phi, psi in pp.get_phi_psi_list():
            if phi and psi:
                phi_psi.append([phi, psi])
    return np.array(phi_psi)

# --- Lógica Principal ---
if pdb_id:
    try:
        # Carga de datos
        structure = get_structure_data(pdb_id)
        
        # Obtener secuencia limpia
        residues = [r for m in structure for c in m for r in c if r.get_id()[0]==' ']
        seq = "".join([seq1(r.get_resname()) for r in residues if seq1(r.get_resname()) != 'X'])
        
        # Objeto de Análisis
        analysed_seq = ProteinAnalysis(seq)
        
        # --- PESTAÑAS PRINCIPALES ---
        tab1, tab2, tab3, tab4 = st.tabs(["🧊 Visor 3D", "📐 Validación (Ramachandran)", "⚗️ Laboratorio de Mutaciones", "📊 Reporte General"])

        # 1. VISOR 3D MEJORADO
        with tab1:
            c1, c2 = st.columns([3, 1])
            with c1:
                view = py3Dmol.view(query='pdb:'+pdb_id)
                # Aplicar estilos dinámicos
                if style_3d == "Cartoon": view.setStyle({'cartoon':{'color':color_3d}})
                elif style_3d == "Stick": view.setStyle({'stick':{'color':color_3d}})
                elif style_3d == "Sphere": view.setStyle({'sphere':{'color':color_3d}})
                elif style_3d == "Line": view.setStyle({'line':{'color':color_3d}})
                
                view.zoomTo()
                showmol(view, height=500, width=800)
            with c2:
                st.info(f"Visualizando: **{pdb_id}**")
                st.write(f"Residuos totales: {len(seq)}")
                st.write("Usa el menú lateral para cambiar el estilo de renderizado.")

        # 2. GRÁFICO DE RAMACHANDRAN (Nivel Experto)
        with tab2:
            st.subheader("Validación Geométrica: Gráfico de Ramachandran")
            st.markdown("Este gráfico muestra los ángulos de torsión (Phi vs Psi) del esqueleto proteico. Las regiones densas indican estructuras secundarias estables (Hélices Alfa / Láminas Beta).")
            
            angles = calculate_ramachandran(structure)
            if len(angles) > 0:
                # Convertir radianes a grados
                angles = angles * 180 / np.pi
                df_rama = pd.DataFrame(angles, columns=['Phi', 'Psi'])
                
                # Gráfico con Altair
                chart = alt.Chart(df_rama).mark_circle(size=60, opacity=0.5).encode(
                    x=alt.X('Phi', title='Phi (φ) Grados', scale=alt.Scale(domain=[-180, 180])),
                    y=alt.Y('Psi', title='Psi (ψ) Grados', scale=alt.Scale(domain=[-180, 180])),
                    tooltip=['Phi', 'Psi'],
                    color=alt.value('purple')
                ).properties(
                    width=600, height=400, title=f"Distribución Conformacional de {pdb_id}"
                ).interactive()
                
                st.altair_chart(chart, use_container_width=True)
            else:
                st.warning("No se pudieron calcular los ángulos para esta estructura.")

        # 3. LABORATORIO DE MUTACIONES (Interactivo)
        with tab3:
            st.subheader("🧬 Simulador de Mutagénesis Dirigida")
            st.write("Modifica un residuo de la secuencia y observa el impacto en las propiedades fisicoquímicas.")
            
            col_mut1, col_mut2 = st.columns(2)
            
            # Selector de mutación
            with col_mut1:
                st.code(seq[:50] + "...", language="text")
                pos = st.number_input("Posición a mutar (1 - " + str(len(seq)) + ")", min_value=1, max_value=len(seq), value=1)
                new_aa = st.selectbox("Nuevo Aminoácido", ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])
                
                # Crear secuencia mutada
                original_aa = seq[pos-1]
                mutated_seq = seq[:pos-1] + new_aa + seq[pos:]
                
                st.write(f"Cambio: **{original_aa}{pos}{new_aa}**")

            # Comparación de resultados
            with col_mut2:
                # Análisis
                p_wt = ProteinAnalysis(seq) # Wild Type
                p_mut = ProteinAnalysis(mutated_seq) # Mutante
                
                # Métricas
                delta_mw = p_mut.molecular_weight() - p_wt.molecular_weight()
                delta_pi = p_mut.isoelectric_point() - p_wt.isoelectric_point()
                
                st.metric("Δ Peso Molecular", f"{delta_mw:.2f} Da", delta_color="inverse")
                st.metric("Δ Punto Isoeléctrico", f"{delta_pi:.2f} pH", delta_color="normal")
                
                if abs(delta_mw) > 50:
                    st.warning("⚠️ Cambio significativo en masa detectado.")
                else:
                    st.success("✅ Cambio estructural leve.")

        # 4. REPORTE GENERAL
        with tab4:
            st.subheader("Resumen Estadístico")
            # Conteo de aminoácidos
            aa_count = {k: seq.count(k) for k in set(seq)}
            df_aa = pd.DataFrame(list(aa_count.items()), columns=['AA', 'Count'])
            
            # Gráfico de pastel
            fig, ax = plt.subplots()
            ax.pie(df_aa['Count'], labels=df_aa['AA'], autopct='%1.1f%%', startangle=90)
            ax.axis('equal') 
            st.pyplot(fig)
            
            st.write("**Clasificación de Estabilidad:**")
            instability = analysed_seq.instability_index()
            if instability > 40:
                st.error(f"Inestable (Índice: {instability:.2f})")
            else:
                st.success(f"Estable (Índice: {instability:.2f})")

    except Exception as e:
        st.error(f"Error procesando la proteína: {e}")
        st.info("Intenta con otro ID (ej: 6LU7, 4HHB)")

else:
    st.info("Por favor ingresa un ID de PDB en el menú lateral.")
