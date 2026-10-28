
import streamlit as st
import pandas as pd
import plotly.express as px

st.set_page_config(
    page_title="BacTrace | 16S rRNA Explorer",
    page_icon="🧬",
    layout="wide"
)

# ---------------- PAGE STYLE ----------------
st.markdown("""
<style>
.main {
    background-color: #f7f9fc;
}
.block-container {
    padding-top: 2rem;
    padding-bottom: 3rem;
}
.hero {
    padding: 28px;
    border-radius: 18px;
    background: linear-gradient(120deg, #102a43, #176b87);
    color: white;
    margin-bottom: 20px;
}
.hero h1 {
    color: white;
    font-size: 42px;
    margin-bottom: 8px;
}
.hero p {
    color: #e1f4fa;
    font-size: 17px;
}
.small-note {
    color: #526477;
    font-size: 13px;
}
</style>
""", unsafe_allow_html=True)

# ---------------- PROJECT DATA ----------------
species_data = pd.DataFrame({
    "Species": [
        "K. pneumoniae",
        "K. variicola",
        "K. quasipneumoniae"
    ],
    "Sequences": [4, 4, 3]
})

alignment_data = pd.DataFrame({
    "Parameter": [
        "Fully conserved sites",
        "Variable sites",
        "Parsimony-informative sites"
    ],
    "Count": [714, 43, 15]
})

identity_data = pd.DataFrame({
    "Comparison type": ["Between-species", "Within-species"],
    "Comparisons": [40, 15],
    "Mean identity (%)": [99.0718, 99.5422],
    "Minimum identity (%)": [97.2936, 98.6450],
    "Maximum identity (%)": [99.7874, 100.0000],
    "Mean comparable positions": [1279.85, 1273.9333]
})

pair_data = pd.DataFrame([
    {
        "Sequence 1": "NR_117686.1",
        "Sequence 2": "AB004753.1",
        "Identity (%)": 100.0000,
        "Comparable positions": 1448,
        "Category": "Most similar pair"
    },
    {
        "Sequence 1": "KR092086.1",
        "Sequence 2": "OM780282.1",
        "Identity (%)": 97.2936,
        "Comparable positions": 739,
        "Category": "Least similar pair"
    }
])

# ---------------- SIDEBAR ----------------
st.sidebar.title("🧬 BacTrace")
st.sidebar.caption("Comparative bacterial sequence analysis")

page = st.sidebar.radio(
    "Navigate",
    [
        "Overview",
        "Dataset",
        "Conservation Analysis",
        "Similarity Explorer",
        "Phylogenetic Analysis",
        "Methods & Limitations"
    ]
)

st.sidebar.divider()
st.sidebar.caption("Academic bioinformatics project")
st.sidebar.caption("16S rRNA | Klebsiella species")

# ---------------- HEADER ----------------
st.markdown("""
<div class="hero">
    <h1>🧬 BacTrace</h1>
    <p>A Comparative 16S rRNA Sequence Explorer</p>
    Explore sequence conservation, nucleotide variation,
    pairwise identity and relationships among closely related
    <i>Klebsiella</i> species.
</div>
""", unsafe_allow_html=True)

# ---------------- OVERVIEW ----------------
if page == "Overview":
    st.subheader("Project Overview")
    st.write(
        "BacTrace is a bioinformatics dashboard for examining "
        "16S rRNA sequences from three closely related bacterial "
        "species: K. pneumoniae, K. variicola and "
        "K. quasipneumoniae."
    )

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total sequences", "11")
    c2.metric("Species studied", "3")
    c3.metric("Alignment length", "1,530 bp")
    c4.metric("Unique comparisons", "55")

    st.divider()
    st.subheader("Sequence Dataset")

    st.dataframe(
        species_data,
        use_container_width=True,
        hide_index=True
    )

    fig = px.pie(
        species_data,
        names="Species",
        values="Sequences",
        hole=0.48,
        title="Distribution of Sequences by Species"
    )
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Key Findings")
    a, b = st.columns(2)

    with a:
        st.info(
            "**Conserved sites:** 714\n\n"
            "Positions with the same unambiguous nucleotide "
            "across all 11 sequences."
        )
        st.info(
            "**Variable sites:** 43\n\n"
            "Positions containing at least two distinct "
            "unambiguous nucleotide states."
        )

    with b:
        st.info(
            "**Parsimony-informative sites:** 15\n\n"
            "Positions with at least two nucleotide states, "
            "each occurring at least twice."
        )
        st.info(
            "**Mean identity:** 99.54% within species and "
            "99.07% between species."
        )

    st.warning(
        "These are descriptive sequence-analysis results. "
        "High 16S rRNA identity does not establish reliable "
        "species-level identification."
    )

# ---------------- DATASET ----------------
elif page == "Dataset":
    st.subheader("Sequence Dataset")
    st.write(
        "The dataset contains 11 accession records selected "
        "from the NCBI nucleotide database."
    )

    accession_data = pd.DataFrame({
        "Species": (
            ["K. pneumoniae"] * 4
            + ["K. variicola"] * 4
            + ["K. quasipneumoniae"] * 3
        ),
        "Accession": [
            "X93214.1", "NR_117686.1",
            "AB004753.1", "AF130981.1",
            "NR_025635.1", "KR092086.1",
            "KC853302.1", "MH542319.1",
            "NR_134062.1", "NR_134063.1",
            "OM780282.1"
        ]
    })

    selected_species = st.multiselect(
        "Filter by species",
        options=accession_data["Species"].unique().tolist(),
        default=accession_data["Species"].unique().tolist()
    )

    filtered = accession_data[
        accession_data["Species"].isin(selected_species)
    ]

    st.dataframe(filtered, use_container_width=True, hide_index=True)

    st.download_button(
        "Download accession table (CSV)",
        data=filtered.to_csv(index=False).encode("utf-8"),
        file_name="BacTrace_accessions.csv",
        mime="text/csv"
    )

    st.markdown(
        "[Open NCBI Nucleotide Database](https://www.ncbi.nlm.nih.gov/nuccore/)"
    )

# ---------------- CONSERVATION ----------------
elif page == "Conservation Analysis":
    st.subheader("Conservation and Variation")

    st.write(
        "The aligned dataset contained 1,530 columns. "
        "The following counts were obtained using the "
        "definitions applied in the sequence analysis."
    )

    st.dataframe(
        alignment_data,
        use_container_width=True,
        hide_index=True
    )

    fig = px.bar(
        alignment_data,
        x="Parameter",
        y="Count",
        text="Count",
        title="Observed Alignment-Site Categories"
    )
    fig.update_layout(xaxis_title="", yaxis_title="Number of sites")
    st.plotly_chart(fig, use_container_width=True)

    st.caption(
        "Percentages using 1,530 alignment columns: "
        "46.67% fully conserved, 2.81% variable and "
        "0.98% parsimony-informative. These categories are "
        "defined separately and do not cover every alignment column."
    )

    st.warning(
        "Variable positions are not automatically species-specific "
        "diagnostic markers."
    )

# ---------------- SIMILARITY ----------------
elif page == "Similarity Explorer":
    st.subheader("Pairwise Sequence Identity")

    st.write(
        "The analysis compared all 55 unique sequence pairs, "
        "excluding gaps and ambiguous bases at non-comparable "
        "positions."
    )

    st.dataframe(
        identity_data,
        use_container_width=True,
        hide_index=True
    )

    fig = px.bar(
        identity_data,
        x="Comparison type",
        y="Mean identity (%)",
        text="Mean identity (%)",
        title="Mean Pairwise Identity by Comparison Type",
        range_y=[95, 101]
    )
    fig.update_traces(texttemplate="%{text:.4f}%", textposition="outside")
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Most and Least Similar Pairs")
    st.dataframe(
        pair_data,
        use_container_width=True,
        hide_index=True
    )

    st.download_button(
        "Download identity summary (CSV)",
        data=identity_data.to_csv(index=False).encode("utf-8"),
        file_name="BacTrace_identity_summary.csv",
        mime="text/csv"
    )

    st.caption(
        "The least similar pair had 739 comparable positions. "
        "Identity values should be interpreted alongside sequence coverage."
    )

# ---------------- PHYLOGENETICS ----------------
elif page == "Phylogenetic Analysis":
    st.subheader("Phylogenetic Relationships")

    st.warning(
        "Tree interpretation pending. The IQ-TREE analysis has "
        "been submitted, but the final tree has not yet been "
        "incorporated into this dashboard."
    )

    st.write(
        "Once the tree is available, this page will present the "
        "observed sequence clusters, branch lengths and branch-support "
        "values, where available."
    )

    st.markdown(
        "[Open IQ-TREE Web Server](https://iqtree.cibiv.univie.ac.at/)"
    )

# ---------------- METHODS ----------------
elif page == "Methods & Limitations":
    st.subheader("Methods")

    st.markdown("""
    1. **Sequence retrieval:** NCBI Nucleotide database.
    2. **Multiple sequence alignment:** MAFFT through SeqBench.
    3. **Sequence statistics:** Biopython supplementary calculations.
    4. **Conservation analysis:** Conserved, variable and
       parsimony-informative alignment positions.
    5. **Pairwise identity:** Comparisons using shared,
       unambiguous nucleotide positions.
    6. **Phylogenetic analysis:** IQ-TREE maximum-likelihood workflow.
    """)

    st.subheader("Limitations")

    st.markdown("""
    - Only 11 sequences were analyzed.
    - Sequence coverage varies among accessions.
    - 16S rRNA may not reliably distinguish all closely related
      Klebsiella species.
    - Accession metadata and unusual alignment statistics should
      be verified before final conclusions.
    - The dashboard presents research findings; it is not a
      clinically validated diagnostic tool.
    """)

    st.subheader("Interpretation")
    st.write(
        "The results indicate high sequence similarity in the "
        "selected dataset. Species-level conclusions require "
        "careful assessment of sequence quality, phylogenetic "
        "evidence and, where appropriate, additional genetic markers."
    )

# ---------------- FOOTER ----------------
st.divider()
st.caption(
    "BacTrace | Educational bioinformatics project | "
    "NCBI • MAFFT • IQ-TREE • Biopython"
)
