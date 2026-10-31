
import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO

# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------

st.set_page_config(
    page_title="BacTrace | Bioinformatics Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --------------------------------------------------
# CUSTOM PROFESSIONAL THEME
# --------------------------------------------------

st.markdown("""
<style>
.stApp {
    background-color: #F4F7FB;
    color: #172B4D;
}

[data-testid="stSidebar"] {
    background-color: #102A43;
}

[data-testid="stSidebar"] * {
    color: white;
}

.main-title {
    font-size: 34px;
    font-weight: 800;
    color: #12304A;
}

.subtitle {
    font-size: 14px;
    color: #627D98;
    margin-bottom: 25px;
}

.metric-card {
    background-color: white;
    padding: 20px;
    border-radius: 12px;
    border: 1px solid #D9E2EC;
    min-height: 115px;
}

.metric-label {
    color: #627D98;
    font-size: 13px;
}

.metric-value {
    font-size: 27px;
    font-weight: bold;
    color: #087E8B;
    margin-top: 8px;
}

.section-title {
    color: #12304A;
    font-size: 21px;
    font-weight: 700;
    margin-top: 20px;
}

.stButton button {
    background-color: #087E8B;
    color: white;
    border-radius: 8px;
    border: none;
}

</style>
""", unsafe_allow_html=True)


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown(
    '<div class="main-title">🧬 BacTrace</div>',
    unsafe_allow_html=True
)

st.markdown(
    '<div class="subtitle">'
    'Comparative 16S rRNA Sequence Analysis and '
    'Bacterial Identification Platform'
    '</div>',
    unsafe_allow_html=True
)


# --------------------------------------------------
# SIDEBAR
# --------------------------------------------------

st.sidebar.title("🧬 BacTrace")
st.sidebar.caption("Bioinformatics Analysis Platform")

page = st.sidebar.radio(
    "NAVIGATION",
    [
        "Dashboard",
        "Sequence Input",
        "Sequence Statistics",
        "Alignment Analysis",
        "Conservation Analysis",
        "Similarity Explorer",
        "Phylogenetic Analysis",
        "Report & Export"
    ]
)

st.sidebar.divider()

st.sidebar.info(
    "Upload bacterial nucleotide sequences "
    "in FASTA format to begin analysis."
)


# --------------------------------------------------
# FASTA INPUT
# --------------------------------------------------

if "records" not in st.session_state:
    st.session_state.records = []

if page == "Sequence Input":

    st.markdown("## Sequence Input")

    uploaded_file = st.file_uploader(
        "Upload FASTA File",
        type=["fasta", "fa", "fna", "txt"]
    )

    pasted_sequence = st.text_area(
        "Or paste FASTA sequences here",
        height=200,
        placeholder=">Sequence_1\nATGCGT..."
    )

    if st.button("Process Sequences"):

        records = []

        try:
            if uploaded_file is not None:
                content = uploaded_file.getvalue().decode("utf-8")
            else:
                content = pasted_sequence

            if not content.strip():
                st.error("Please upload or paste a FASTA sequence.")
            else:
                records = list(
                    SeqIO.parse(StringIO(content), "fasta")
                )

                if not records:
                    st.error("No valid FASTA records found.")

                else:
                    valid_records = []

                    for record in records:
                        seq = str(record.seq).upper()

                        if not seq:
                            continue

                        if not set(seq).issubset(set("ACGTN")):
                            st.warning(
                                f"Skipped {record.id}: "
                                "contains unsupported nucleotide symbols."
                            )
                            continue

                        valid_records.append(record)

                    if valid_records:
                        st.session_state.records = valid_records
                        st.success(
                            f"{len(valid_records)} sequences processed successfully."
                        )
                    else:
                        st.error("No valid nucleotide sequences available.")

        except Exception as e:
            st.error(f"Error processing FASTA file: {e}")


# --------------------------------------------------
# CALCULATE REAL SEQUENCE STATISTICS
# --------------------------------------------------

records = st.session_state.records

stats_data = []

for record in records:

    seq = str(record.seq).upper()

    length = len(seq)

    gc_count = seq.count("G") + seq.count("C")

    atgc_count = sum(seq.count(base) for base in "ATGC")

    gc_content = (
        gc_count / atgc_count * 100
        if atgc_count > 0 else 0
    )

    stats_data.append({
        "Sequence": record.id,
        "Length": length,
        "GC_Content": round(gc_content, 2),
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C"),
        "N": seq.count("N")
    })

stats_df = pd.DataFrame(stats_data)


# --------------------------------------------------
# DASHBOARD
# --------------------------------------------------

if page == "Dashboard":

    st.markdown("## Analysis Overview")

    if stats_df.empty:

        st.info(
            "No sequence data available. "
            "Go to Sequence Input and upload a FASTA file."
        )

    else:

        total_sequences = len(stats_df)
        total_bases = int(stats_df["Length"].sum())
        avg_length = stats_df["Length"].mean()
        avg_gc = stats_df["GC_Content"].mean()

        c1, c2, c3, c4 = st.columns(4)

        metrics = [
            (c1, "Total Sequences", total_sequences),
            (c2, "Total Nucleotides", f"{total_bases:,}"),
            (c3, "Average Length", f"{avg_length:,.0f} bp"),
            (c4, "Average GC Content", f"{avg_gc:.2f}%")
        ]

        for col, label, value in metrics:
            with col:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-label">{label}</div>
                    <div class="metric-value">{value}</div>
                </div>
                """, unsafe_allow_html=True)

        st.markdown(
            '<div class="section-title">Sequence Analysis</div>',
            unsafe_allow_html=True
        )

        left, right = st.columns(2)

        with left:

            fig = px.bar(
                stats_df,
                x="Sequence",
                y="Length",
                title="Sequence Length Comparison",
                labels={
                    "Length": "Length (bp)",
                    "Sequence": "Sequence ID"
                }
            )

            fig.update_layout(
                plot_bgcolor="white",
                paper_bgcolor="white"
            )

            st.plotly_chart(fig, use_container_width=True)

        with right:

            fig = px.bar(
                stats_df,
                x="Sequence",
                y="GC_Content",
                title="GC Content Distribution",
                labels={
                    "GC_Content": "GC Content (%)",
                    "Sequence": "Sequence ID"
                }
            )

            fig.update_layout(
                plot_bgcolor="white",
                paper_bgcolor="white"
            )

            st.plotly_chart(fig, use_container_width=True)

        st.markdown("### Uploaded Sequence Summary")

        st.dataframe(
            stats_df,
            use_container_width=True,
            hide_index=True
        )


# --------------------------------------------------
# SEQUENCE STATISTICS
# --------------------------------------------------

elif page == "Sequence Statistics":

    st.header("Sequence Statistics")

    if stats_df.empty:
        st.warning("Upload sequences first.")

    else:
        st.dataframe(stats_df, use_container_width=True)

        st.markdown("### Nucleotide Composition")

        composition = {
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [
                int(stats_df["A"].sum()),
                int(stats_df["T"].sum()),
                int(stats_df["G"].sum()),
                int(stats_df["C"].sum())
            ]
        }

        composition_df = pd.DataFrame(composition)

        fig = px.pie(
            composition_df,
            names="Nucleotide",
            values="Count",
            title="Overall Nucleotide Composition"
        )

        st.plotly_chart(fig, use_container_width=True)


# --------------------------------------------------
# ALIGNMENT ANALYSIS
# --------------------------------------------------

elif page == "Alignment Analysis":

    st.header("Sequence Alignment Analysis")

    st.info(
        "This section is reserved for the multiple sequence "
        "alignment output generated by the BacTrace alignment module."
    )

    st.write(
        "Alignment algorithms: MAFFT or a documented "
        "educational fallback."
    )


# --------------------------------------------------
# CONSERVATION ANALYSIS
# --------------------------------------------------

elif page == "Conservation Analysis":

    st.header("Conservation Analysis")

    st.info(
        "This section displays nucleotide conservation "
        "across aligned sequences."
    )

    st.write(
        "Conservation analysis requires multiple sequence "
        "alignment before calculating conserved positions."
    )


# --------------------------------------------------
# SIMILARITY EXPLORER
# --------------------------------------------------

elif page == "Similarity Explorer":

    st.header("Similarity Explorer")

    st.info(
        "Pairwise sequence similarity analysis will be "
        "displayed here using the BacTrace comparison module."
    )


# --------------------------------------------------
# PHYLOGENETIC ANALYSIS
# --------------------------------------------------

elif page == "Phylogenetic Analysis":

    st.header("Phylogenetic Analysis")

    st.info(
        "Phylogenetic tree generation requires aligned "
        "sequences and an appropriate distance matrix."
    )

    st.write(
        "Neighbor-Joining can be used for preliminary "
        "educational phylogenetic analysis."
    )


# --------------------------------------------------
# REPORT AND EXPORT
# --------------------------------------------------

elif page == "Report & Export":

    st.header("Report & Export")

    if stats_df.empty:

        st.warning("No sequence data available for export.")

    else:

        csv = stats_df.to_csv(index=False).encode("utf-8")

        st.download_button(
            "Download Sequence Statistics (CSV)",
            data=csv,
            file_name="BacTrace_Sequence_Statistics.csv",
            mime="text/csv"
        )

        fasta_output = ""

        for record in records:
            fasta_output += (
                f">{record.id}\n{str(record.seq)}\n"
            )

        st.download_button(
            "Download FASTA Sequences",
            data=fasta_output,
            file_name="BacTrace_Sequences.fasta",
            mime="text/plain"
        )


# --------------------------------------------------
# FOOTER
# --------------------------------------------------

st.divider()

st.caption(
    "BacTrace | Comparative 16S rRNA Sequence Analysis "
    "| Academic Bioinformatics Project"
)

st.caption(
    "Disclaimer: 16S rRNA sequence similarity provides "
    "preliminary taxonomic evidence and does not independently "
    "confirm bacterial species identity."
)
