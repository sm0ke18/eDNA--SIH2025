import streamlit as st
import pandas as pd
import plotly.express as px
from edna_utils import parse_fasta, shannon_index, simpson_index, chao1, mock_blast, run_blast, get_taxonomy_ete
# --- Patch for Python 3.13 (cgi removed) ---
import sys, types
if "cgi" not in sys.modules:
    cgi = types.ModuleType("cgi")
    cgi.escape = lambda s, quote=True: s  # old escape replacement
    sys.modules["cgi"] = cgi

from ete3 import NCBITaxa



# ---------------- Streamlit Layout ----------------
st.set_page_config(page_title="eDNA Biodiversity Explorer", layout="wide")
st.title("🧬 eDNA Taxonomy & Biodiversity Assessment")

# ---------------- File Upload ----------------
uploaded_file = st.file_uploader("📂 Upload FASTA file", type=["fa", "fasta", "txt"])

# Mode selector: mock vs real BLAST
blast_mode = st.radio("Choose BLAST mode", ["Mock BLAST (fast)", "NCBI BLAST (real, slow)"], index=0)

# Initialize NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

if uploaded_file:
    fasta_str = uploaded_file.read().decode("utf-8")
    sequences = parse_fasta(fasta_str)

    results = []
    for seq_id, seq in sequences.items():
        if blast_mode == "Mock BLAST (fast)":
            hit = mock_blast(seq)
        else:
            hit = run_blast(seq)

        results.append({
            "seq_id": seq_id,
            "sequence_length": len(seq),
            "common_name": hit["species"],
            "hit_def": hit["hit_def"],
            "accession": hit["accession"]
        })

    # ---------------- Results Table ----------------
    results_df = pd.DataFrame(results)
    st.subheader("📊 Results Table")
    st.dataframe(results_df, use_container_width=True)

    # ---------------- Species Distribution ----------------
    st.subheader("🦠 Species Distribution")
    species_counts = results_df['common_name'].value_counts()
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
        taxonomy = get_taxonomy_ete(row['common_name'])
        taxonomy_records.append({
            "kingdom": taxonomy["kingdom"],
            "phylum": taxonomy["phylum"],
            "class": taxonomy["class"],
            "common_name": row['common_name']
        })
    tax_df = pd.DataFrame(taxonomy_records)
    fig = px.sunburst(tax_df, path=['kingdom', 'phylum', 'class', 'common_name'], title="Taxonomic Breakdown")
    st.plotly_chart(fig, use_container_width=True)

    # ---------------- CSV Download ----------------
    csv = results_df.to_csv(index=False).encode("utf-8")
    st.download_button("⬇️ Download Results CSV", csv, "results.csv", "text/csv")
