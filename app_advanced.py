import streamlit as st
import pandas as pd
import plotly.express as px
from edna_utils import parse_fasta, shannon_index, simpson_index, chao1, mock_blast, run_blast, get_taxonomy_ete
from ete3 import NCBITaxa

# ---------------- Streamlit Layout ----------------
st.set_page_config(page_title="eDNA Biodiversity Explorer", layout="wide")
st.title("🧬 eDNA Taxonomy & Biodiversity Assessment")

# ---------------- File Upload ----------------
uploaded_file = st.file_uploader("📂 Upload FASTA file", type=["fa", "fasta", "txt"])

# Mode selector: mock vs real BLAST
blast_mode = st.radio("Choose BLAST mode", ["Mock BLAST (fast)", "NCBI BLAST (real, slow)"], index=0)

# Initialize NCBITaxa (do NOT call update_taxonomy_database() every run)
ncbi = NCBITaxa()

if uploaded_file:
    fasta_str = uploaded_file.read().decode("utf-8")
    sequences = parse_fasta(fasta_str)

    results = []
    for seq_id, seq in sequences.items():
        if blast_mode == "Mock BLAST (fast)":
            hit = mock_blast(seq)
        else:
            hit = run_blast(seq)  # This uses NCBIWWW.qblast

        results.append({
            "seq_id": seq_id,
            "sequence_length": len(seq),
            "species": hit["species"],  # Use species from BLAST result
            "hit_def": hit["hit_def"],
            "accession": hit["accession"]
        })

    # ---------------- Results Table ----------------
    results_df = pd.DataFrame(results)
    st.subheader("📊 Results Table")
    st.dataframe(results_df, use_container_width=True)

    # ---------------- Species Distribution ----------------
    st.subheader("🦠 Species Distribution")
    species_counts = results_df['species'].value_counts()
    st.bar_chart(species_counts)

    # ---------------- Biodiversity Indices ----------------
    counts_dict = species_counts.to_dict()
    shannon = shannon_index(counts_dict)
    simpson = simpson_index(counts_dict)
    chao = chao1(counts_dict)

    st.subheader("🌍 Biodiversity Indices")
    st.markdown(f"""
    - **Shannon Index:** {shannon:.3f}  
    - **Simpson Index:** {simpson:.3f}  
    - **Chao1 Richness Estimator:** {chao:.1f}
    """)

    # ---------------- Taxonomic Hierarchy (Sunburst) ----------------
    st.subheader("🌳 Taxonomic Hierarchy (Sunburst)")
    taxonomy_records = []
    for _, row in results_df.iterrows():
        taxonomy = get_taxonomy_ete(row['species'])
        taxonomy_records.append({
            "kingdom": taxonomy["kingdom"],
            "phylum": taxonomy["phylum"],
            "class": taxonomy["class"],
            "species": row['species']
        })
    tax_df = pd.DataFrame(taxonomy_records)
    fig = px.sunburst(tax_df, path=['kingdom', 'phylum', 'class', 'species'], title="Taxonomic Breakdown")
    st.plotly_chart(fig, use_container_width=True)

    # ---------------- CSV Download ----------------
    csv = results_df.to_csv(index=False).encode("utf-8")
    st.download_button("⬇️ Download Results CSV", csv, "results.csv", "text/csv")

# To run the app: streamlit run app_advanced.py
# If using a virtual environment, activate it first: venv312\Scripts\activate

# Required installations
# !pip install streamlit pandas plotly biopython ete3
