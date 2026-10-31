
import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO

# ==================================================
# PAGE CONFIGURATION
# ==================================================

st.set_page_config(
    page_title="BacTrace - 16S rRNA Analysis",
    page_icon="🧬",
    layout="wide"
)

# ==================================================
# SESSION STATE
# ==================================================

if "records" not in st.session_state:
    st.session_state.records = []

# ==================================================
# HEADER
# ==================================================

st.title("🧬 BacTrace")

st.subheader(
    "Comparative 16S rRNA Sequence Analysis "
    "and Bacterial Identification Platform"
)

st.caption(
    "A bioinformatics platform for nucleotide sequence "
    "statistics, comparison and preliminary bacterial identification."
)

st.divider()

# ==================================================
# NAVIGATION
# ==================================================

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "Dashboard",
    "Sequence Input",
    "Statistics",
    "Alignment",
    "Phylogeny",
    "Reports"
])

# ==================================================
# FASTA PROCESSING
# ==================================================

def process_fasta(content):

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
# STATISTICS CALCULATION
# ==================================================

def calculate_statistics(records):

    data = []

    for record in records:

        sequence = str(record.seq).upper()

        a = sequence.count("A")
        t = sequence.count("T")
        g = sequence.count("G")
        c = sequence.count("C")
        n = sequence.count("N")

        length = len(sequence)

        valid_bases = a + t + g + c

        gc_content = (
            ((g + c) / valid_bases) * 100
            if valid_bases > 0 else 0
        )

        data.append({
            "Sequence ID": record.id,
            "Length (bp)": length,
            "GC Content (%)": round(gc_content, 2),
            "A": a,
            "T": t,
            "G": g,
            "C": c,
            "N": n
        })

    return pd.DataFrame(data)


# ==================================================
# SEQUENCE INPUT TAB
# ==================================================

with tab2:

    st.header("Sequence Input")

    st.write(
        "Upload your bacterial 16S rRNA sequences "
        "or paste FASTA data below."
    )

    uploaded_file = st.file_uploader(
        "Choose FASTA file",
        type=["fasta", "fa", "fna", "txt"]
    )

    pasted_sequence = st.text_area(
        "Paste FASTA sequence",
        height=180,
        placeholder=">Sequence_1\nATGCGT..."
    )

    if st.button("Process Sequences", type="primary"):

        try:

            if uploaded_file is not None:

                content = uploaded_file.getvalue().decode("utf-8")

            else:

                content = pasted_sequence

            if not content.strip():

                st.error("Please upload or paste a FASTA sequence.")

            else:

                new_records = process_fasta(content)

                if not new_records:

                    st.error(
                        "No valid FASTA sequences found. "
                        "Check the nucleotide sequence format."
                    )

                else:

                    st.session_state.records = new_records

                    st.success(
                        f"{len(new_records)} sequences processed successfully."
                    )

        except Exception as error:

            st.error(f"Processing error: {error}")


# ==================================================
# COMMON DATA
# ==================================================

records = st.session_state.records

if records:

    stats_df = calculate_statistics(records)

else:

    stats_df = pd.DataFrame()


# ==================================================
# DASHBOARD TAB
# ==================================================

with tab1:

    st.header("Analysis Dashboard")

    if stats_df.empty:

        st.info(
            "Welcome to BacTrace. Upload FASTA sequences "
            "using the Sequence Input tab to begin your analysis."
        )

        st.subheader("Available Analysis Modules")

        c1, c2, c3 = st.columns(3)

        with c1:
            st.markdown("**Sequence Statistics**")
            st.write(
                "Calculate sequence length, GC content "
                "and nucleotide composition."
            )

        with c2:
            st.markdown("**Sequence Alignment**")
            st.write(
                "Compare nucleotide sequences "
                "and identify matching regions."
            )

        with c3:
            st.markdown("**Phylogenetic Analysis**")
            st.write(
                "Explore preliminary evolutionary "
                "relationships between sequences."
            )

    else:

        total_sequences = len(stats_df)

        total_bases = int(
            stats_df["Length (bp)"].sum()
        )

        average_length = stats_df["Length (bp)"].mean()

        average_gc = stats_df["GC Content (%)"].mean()

        st.caption(
            f"Dataset contains {total_sequences} valid sequences."
        )

        # METRIC CARDS

        c1, c2, c3, c4 = st.columns(4)

        c1.metric(
            "Total Sequences",
            total_sequences
        )

        c2.metric(
            "Total Nucleotides",
            f"{total_bases:,}"
        )

        c3.metric(
            "Average Length",
            f"{average_length:,.0f} bp"
        )

        c4.metric(
            "Average GC Content",
            f"{average_gc:.2f}%"
        )

        st.divider()

        # SEQUENCE LENGTH CHART

        st.subheader("Sequence Length Comparison")

        fig1 = px.bar(
            stats_df,
            x="Sequence ID",
            y="Length (bp)",
            title="Length of Uploaded Sequences"
        )

        st.plotly_chart(
            fig1,
            use_container_width=True
        )

        # GC CONTENT CHART

        st.subheader("GC Content Distribution")

        fig2 = px.bar(
            stats_df,
            x="Sequence ID",
            y="GC Content (%)",
            title="GC Content Across Sequences"
        )

        st.plotly_chart(
            fig2,
            use_container_width=True
        )

        # NUCLEOTIDE COMPOSITION

        st.subheader("Nucleotide Composition")

        composition_df = pd.DataFrame({
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [
                int(stats_df["A"].sum()),
                int(stats_df["T"].sum()),
                int(stats_df["G"].sum()),
                int(stats_df["C"].sum())
            ]
        })

        fig3 = px.pie(
            composition_df,
            names="Nucleotide",
            values="Count",
            title="Overall Nucleotide Composition",
            hole=0.4
        )

        st.plotly_chart(
            fig3,
            use_container_width=True
        )

        st.subheader("Uploaded Sequence Summary")

        st.dataframe(
            stats_df,
            use_container_width=True,
            hide_index=True
        )


# ==================================================
# STATISTICS TAB
# ==================================================

with tab3:

    st.header("Sequence Statistics")

    if stats_df.empty:

        st.info("Upload FASTA sequences to view statistics.")

    else:

        st.dataframe(
            stats_df,
            use_container_width=True,
            hide_index=True
        )

        selected_sequence = st.selectbox(
            "Select a sequence",
            stats_df["Sequence ID"].tolist()
        )

        selected_row = stats_df[
            stats_df["Sequence ID"] == selected_sequence
        ].iloc[0]

        nucleotide_df = pd.DataFrame({
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [
                selected_row["A"],
                selected_row["T"],
                selected_row["G"],
                selected_row["C"]
            ]
        })

        st.subheader("Individual Nucleotide Composition")

        fig = px.bar(
            nucleotide_df,
            x="Nucleotide",
            y="Count"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )


# ==================================================
# ALIGNMENT TAB
# ==================================================

with tab4:

    st.header("Multiple Sequence Alignment")

    st.write(
        "This module compares uploaded sequences "
        "to identify conserved and variable nucleotide positions."
    )

    st.info(
        "Connect your existing MAFFT or Biopython alignment "
        "function here to display actual alignment results."
    )


# ==================================================
# PHYLOGENY TAB
# ==================================================

with tab5:

    st.header("Phylogenetic Analysis")

    st.write(
        "Phylogenetic analysis helps explore evolutionary "
        "relationships among bacterial sequences."
    )

    st.info(
        "Connect your Neighbor-Joining or IQ-TREE "
        "analysis function here to generate a phylogenetic tree."
    )


# ==================================================
# REPORTS TAB
# ==================================================

with tab6:

    st.header("Reports and Data Export")

    if stats_df.empty:

        st.info("Upload sequences to generate downloadable reports.")

    else:

        st.subheader("Download Sequence Statistics")

        csv_data = stats_df.to_csv(
            index=False
        ).encode("utf-8")

        st.download_button(
            "Download CSV Report",
            data=csv_data,
            file_name="BacTrace_Sequence_Statistics.csv",
            mime="text/csv"
        )

        fasta_output = ""

        for record in records:

            fasta_output += (
                f">{record.id}\n"
                f"{str(record.seq)}\n"
            )

        st.download_button(
            "Download FASTA Dataset",
            data=fasta_output,
            file_name="BacTrace_Sequences.fasta",
            mime="text/plain"
        )


# ==================================================
# FOOTER
# ==================================================

st.divider()

st.caption(
    "BacTrace | Academic Bioinformatics Platform | "
    "Comparative 16S rRNA Sequence Analysis"
)

st.caption(
    "Note: 16S rRNA similarity alone does not establish "
    "definitive bacterial species identity."
)
