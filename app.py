
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO
from io import StringIO

# ==================================================
# BACTrace - PAGE CONFIGURATION
# ==================================================

st.set_page_config(
    page_title="BacTrace | 16S rRNA Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ==================================================
# CUSTOM CSS
# ==================================================

st.markdown("""
<style>

.stApp {
    background-color: #F5F8FC;
    color: #18324B;
}

header[data-testid="stHeader"] {
    background-color: #F5F8FC;
}

.block-container {
    padding-top: 1.5rem;
    padding-bottom: 2rem;
    max-width: 1500px;
}

.main-header {
    background: linear-gradient(110deg, #102A43, #087E8B);
    padding: 30px 35px;
    border-radius: 15px;
    color: white;
    margin-bottom: 22px;
}

.main-header h1 {
    color: white;
    font-size: 34px;
    margin: 0;
    font-weight: 800;
}

.main-header p {
    color: #D9F3F1;
    font-size: 14px;
    margin-top: 8px;
}

.metric-card {
    background: white;
    border: 1px solid #E1E8EF;
    border-radius: 13px;
    padding: 20px;
    min-height: 120px;
    box-shadow: 0px 3px 12px rgba(16,42,67,0.04);
}

.metric-label {
    font-size: 13px;
    color: #627D98;
    font-weight: 500;
}

.metric-value {
    font-size: 27px;
    font-weight: 750;
    color: #087E8B;
    margin-top: 10px;
}

.section-heading {
    font-size: 21px;
    font-weight: 700;
    color: #18324B;
    margin-top: 20px;
    margin-bottom: 12px;
}

.info-card {
    background: white;
    border: 1px solid #E1E8EF;
    border-radius: 12px;
    padding: 18px;
    min-height: 115px;
}

.info-card h4 {
    color: #087E8B;
    margin-bottom: 8px;
}

.info-card p {
    color: #627D98;
    font-size: 13px;
}

div[data-testid="stTabs"] button {
    font-weight: 600;
    font-size: 14px;
}

div[data-testid="stFileUploader"] {
    background: white;
    border-radius: 12px;
}

div[data-testid="stDataFrame"] {
    border-radius: 10px;
}

.footer {
    text-align: center;
    color: #829AB1;
    font-size: 12px;
    padding: 20px;
}

</style>
""", unsafe_allow_html=True)


# ==================================================
# SESSION STATE
# ==================================================

if "records" not in st.session_state:
    st.session_state.records = []

if "active_page" not in st.session_state:
    st.session_state.active_page = "Dashboard"


# ==================================================
# HEADER
# ==================================================

st.markdown("""
<div class="main-header">
    <h1>🧬 BacTrace</h1>
    <p>
    Comparative 16S rRNA Sequence Analysis and
    Bacterial Identification Platform
    </p>
</div>
""", unsafe_allow_html=True)


# ==================================================
# TOP NAVIGATION
# ==================================================

pages = [
    "Dashboard",
    "Sequence Input",
    "Statistics",
    "Alignment",
    "Conservation",
    "Similarity",
    "Phylogeny",
    "Reports"
]

selected_page = st.radio(
    "Navigation",
    pages,
    horizontal=True,
    label_visibility="collapsed",
    key="navigation"
)

st.session_state.active_page = selected_page

st.divider()


# ==================================================
# FASTA INPUT FUNCTION
# ==================================================

def load_fasta(content):

    records = list(
        SeqIO.parse(StringIO(content), "fasta")
    )

    valid_records = []

    for record in records:

        sequence = str(record.seq).upper()

        if not sequence:
            continue

        if not set(sequence).issubset(set("ACGTN")):
            continue

        valid_records.append(record)

    return valid_records


# ==================================================
# SEQUENCE STATISTICS FUNCTION
# ==================================================

def calculate_statistics(records):

    data = []

    for record in records:

        seq = str(record.seq).upper()

        a = seq.count("A")
        t = seq.count("T")
        g = seq.count("G")
        c = seq.count("C")
        n = seq.count("N")

        length = len(seq)

        valid_bases = a + t + g + c

        gc = (
            (g + c) / valid_bases * 100
            if valid_bases else 0
        )

        data.append({
            "Sequence ID": record.id,
            "Length (bp)": length,
            "GC Content (%)": round(gc, 2),
            "A": a,
            "T": t,
            "G": g,
            "C": c,
            "N": n
        })

    return pd.DataFrame(data)


# ==================================================
# SEQUENCE INPUT PAGE
# ==================================================

if selected_page == "Sequence Input":

    st.markdown(
        '<div class="section-heading">Sequence Input</div>',
        unsafe_allow_html=True
    )

    st.write(
        "Upload or paste bacterial 16S rRNA nucleotide "
        "sequences in FASTA format."
    )

    col1, col2 = st.columns(2)

    with col1:

        uploaded_file = st.file_uploader(
            "Upload FASTA file",
            type=["fasta", "fa", "fna", "txt"]
        )

    with col2:

        st.info(
            "Supported formats: FASTA (.fasta, .fa, .fna, .txt)"
        )

    pasted = st.text_area(
        "Or paste FASTA sequence",
        height=180,
        placeholder=">Sequence_1\nATGCGT..."
    )

    if st.button("Process Sequences", type="primary"):

        try:

            if uploaded_file is not None:

                content = uploaded_file.getvalue().decode("utf-8")

            else:

                content = pasted

            if not content.strip():

                st.error("Please provide sequence data.")

            else:

                new_records = load_fasta(content)

                if not new_records:

                    st.error(
                        "No valid nucleotide FASTA records found. "
                        "Check the headers and nucleotide symbols."
                    )

                else:

                    st.session_state.records = new_records

                    st.success(
                        f"{len(new_records)} sequences loaded successfully."
                    )

        except Exception as error:

            st.error(f"Sequence processing error: {error}")


# ==================================================
# COMMON DATA
# ==================================================

records = st.session_state.records

if records:

    stats_df = calculate_statistics(records)

else:

    stats_df = pd.DataFrame()


# ==================================================
# DASHBOARD PAGE
# ==================================================

if selected_page == "Dashboard":

    st.markdown(
        '<div class="section-heading">Analysis Overview</div>',
        unsafe_allow_html=True
    )

    if stats_df.empty:

        st.info(
            "No sequences uploaded yet. Open the Sequence Input "
            "tab to upload your bacterial FASTA data."
        )

        st.markdown("### Platform Modules")

        modules = [
            ("Sequence Statistics",
             "Length, nucleotide composition and GC content."),
            ("Multiple Alignment",
             "Compare nucleotide sequences across samples."),
            ("Conservation Analysis",
             "Identify conserved and variable nucleotide positions."),
            ("Similarity Explorer",
             "Explore pairwise sequence similarity."),
            ("Phylogenetic Analysis",
             "Investigate preliminary evolutionary relationships."),
            ("Report Generation",
             "Export sequence analysis results.")
        ]

        cols = st.columns(3)

        for i, (title, description) in enumerate(modules):

            with cols[i % 3]:

                st.markdown(f"""
                <div class="info-card">
                    <h4>{title}</h4>
                    <p>{description}</p>
                </div>
                """, unsafe_allow_html=True)

    else:

        total_sequences = len(stats_df)

        total_bases = int(
            stats_df["Length (bp)"].sum()
        )

        average_length = stats_df["Length (bp)"].mean()

        average_gc = stats_df["GC Content (%)"].mean()

        st.caption(
            f"Loaded dataset: {total_sequences} valid sequences"
        )

        # METRIC CARDS

        c1, c2, c3, c4 = st.columns(4)

        metrics = [
            (c1, "Total Sequences", total_sequences),
            (c2, "Total Nucleotides", f"{total_bases:,}"),
            (c3, "Average Length", f"{average_length:,.0f} bp"),
            (c4, "Average GC Content", f"{average_gc:.2f}%")
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
            '<div class="section-heading">Sequence Data Visualization</div>',
            unsafe_allow_html=True
        )

        chart1, chart2 = st.columns(2)

        with chart1:

            fig = px.bar(
                stats_df,
                x="Sequence ID",
                y="Length (bp)",
                title="Sequence Length Comparison",
                color="Length (bp)",
                color_continuous_scale="Teal"
            )

            fig.update_layout(
                plot_bgcolor="white",
                paper_bgcolor="white",
                coloraxis_showscale=False
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )

        with chart2:

            fig = px.bar(
                stats_df,
                x="Sequence ID",
                y="GC Content (%)",
                title="GC Content Distribution",
                color="GC Content (%)",
                color_continuous_scale="Teal"
            )

            fig.update_layout(
                plot_bgcolor="white",
                paper_bgcolor="white",
                coloraxis_showscale=False
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )

        st.markdown(
            '<div class="section-heading">Nucleotide Composition</div>',
            unsafe_allow_html=True
        )

        total_a = int(stats_df["A"].sum())
        total_t = int(stats_df["T"].sum())
        total_g = int(stats_df["G"].sum())
        total_c = int(stats_df["C"].sum())

        composition_df = pd.DataFrame({
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [
                total_a,
                total_t,
                total_g,
                total_c
            ]
        })

        fig = px.pie(
            composition_df,
            names="Nucleotide",
            values="Count",
            hole=0.55,
            title="Overall Nucleotide Composition"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

        st.markdown("### Uploaded Sequence Summary")

        st.dataframe(
            stats_df,
            use_container_width=True,
            hide_index=True
        )


# ==================================================
# STATISTICS PAGE
# ==================================================

elif selected_page == "Statistics":

    st.title("Sequence Statistics")

    if stats_df.empty:

        st.warning("Upload FASTA sequences first.")

    else:

        st.dataframe(
            stats_df,
            use_container_width=True,
            hide_index=True
        )

        st.markdown("### Individual Sequence Composition")

        selected_sequence = st.selectbox(
            "Select sequence",
            stats_df["Sequence ID"].tolist()
        )

        row = stats_df[
            stats_df["Sequence ID"] == selected_sequence
        ].iloc[0]

        nucleotide_df = pd.DataFrame({
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [
                row["A"],
                row["T"],
                row["G"],
                row["C"]
            ]
        })

        fig = px.bar(
            nucleotide_df,
            x="Nucleotide",
            y="Count",
            title=f"Nucleotide Composition: {selected_sequence}"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )


# ==================================================
# ALIGNMENT PAGE
# ==================================================

elif selected_page == "Alignment":

    st.title("Multiple Sequence Alignment")

    st.write(
        "Alignment analysis compares nucleotide positions "
        "across the uploaded 16S rRNA sequences."
    )

    st.info(
        "Connect this section to the MAFFT or Biopython "
        "alignment function in your full BacTrace application."
    )


# ==================================================
# CONSERVATION PAGE
# ==================================================

elif selected_page == "Conservation":

    st.title("Conservation Analysis")

    st.write(
        "Conservation analysis identifies nucleotide positions "
        "that remain similar or vary among aligned sequences."
    )

    st.info(
        "Conservation scores require aligned sequence data."
    )


# ==================================================
# SIMILARITY PAGE
# ==================================================

elif selected_page == "Similarity":

    st.title("Pairwise Similarity Explorer")

    st.write(
        "Compare sequence pairs to investigate nucleotide "
        "identity and sequence similarity."
    )

    st.info(
        "Connect this section to the pairwise identity "
        "calculation from your BacTrace analysis module."
    )


# ==================================================
# PHYLOGENY PAGE
# ==================================================

elif selected_page == "Phylogeny":

    st.title("Phylogenetic Analysis")

    st.write(
        "Phylogenetic analysis provides a preliminary view "
        "of evolutionary relationships among bacterial sequences."
    )

    st.info(
        "Connect this section to the Neighbor-Joining "
        "or IQ-TREE analysis module."
    )


# ==================================================
# REPORT PAGE
# ==================================================

elif selected_page == "Reports":

    st.title("Reports & Export")

    if stats_df.empty:

        st.warning("Upload sequences before generating reports.")

    else:

        st.success("Sequence statistics are ready for export.")

        csv_data = stats_df.to_csv(
            index=False
        ).encode("utf-8")

        st.download_button(
            "Download Statistics CSV",
            data=csv_data,
            file_name="BacTrace_Statistics.csv",
            mime="text/csv",
            type="primary"
        )

        fasta_text = ""

        for record in records:

            fasta_text += (
                f">{record.id}\n{str(record.seq)}\n"
            )

        st.download_button(
            "Download FASTA Dataset",
            data=fasta_text,
            file_name="BacTrace_Sequences.fasta",
            mime="text/plain"
        )


# ==================================================
# FOOTER
# ==================================================

st.divider()

st.markdown("""
<div class="footer">
BacTrace | Comparative 16S rRNA Sequence Analysis Platform<br>
Academic Bioinformatics Project<br>
16S rRNA similarity alone does not establish definitive species identity.
</div>
""", unsafe_allow_html=True)
