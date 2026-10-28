
import streamlit as st
import pandas as pd
import plotly.express as px

# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------
st.set_page_config(
    page_title="BacTrace | 16S rRNA Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --------------------------------------------------
# PROJECT DATA
# --------------------------------------------------
PAGES = [
    "Home",
    "Research Overview",
    "Dataset",
    "Conservation",
    "Similarity",
    "Phylogenetic Analysis",
    "Methods & References"
]

ACCESSIONS = [
    {"Species": "Klebsiella pneumoniae", "Accession": "X93214.1"},
    {"Species": "Klebsiella pneumoniae", "Accession": "NR_117686.1"},
    {"Species": "Klebsiella pneumoniae", "Accession": "AB004753.1"},
    {"Species": "Klebsiella pneumoniae", "Accession": "AF130981.1"},
    {"Species": "Klebsiella variicola", "Accession": "NR_025635.1"},
    {"Species": "Klebsiella variicola", "Accession": "KR092086.1"},
    {"Species": "Klebsiella variicola", "Accession": "KC853302.1"},
    {"Species": "Klebsiella variicola", "Accession": "MH542319.1"},
    {"Species": "Klebsiella quasipneumoniae", "Accession": "NR_134062.1"},
    {"Species": "Klebsiella quasipneumoniae", "Accession": "NR_134063.1"},
    {"Species": "Klebsiella quasipneumoniae", "Accession": "OM780282.1"},
]

accession_df = pd.DataFrame(ACCESSIONS)

# Results from the completed alignment analysis
alignment_length = 1530
sequence_count = 11
fully_conserved = 714
variable_sites = 43
parsimony_informative = 15

similarity_summary = pd.DataFrame([
    {
        "Comparison": "Within species",
        "Number of pairs": 15,
        "Mean identity (%)": 99.5422,
        "Minimum identity (%)": 98.6450,
        "Maximum identity (%)": 100.0000,
        "Mean comparable positions": 1273.93,
    },
    {
        "Comparison": "Between species",
        "Number of pairs": 40,
        "Mean identity (%)": 99.0718,
        "Minimum identity (%)": 97.2936,
        "Maximum identity (%)": 99.7874,
        "Mean comparable positions": 1279.85,
    },
])

# --------------------------------------------------
# SESSION STATE / PAGE NAVIGATION
# --------------------------------------------------
if "page_index" not in st.session_state:
    st.session_state.page_index = 0


def go_to_page(index):
    st.session_state.page_index = index


# --------------------------------------------------
# CUSTOM DESIGN
# --------------------------------------------------
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Manrope:wght@400;500;600;700;800&display=swap');

:root {
    --navy: #10243A;
    --teal: #087E8B;
    --teal-light: #E7F5F5;
    --ink: #203449;
    --muted: #64748B;
    --line: #DFE7ED;
    --paper: #F6F9FB;
}

html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
}

.stApp {
    background: var(--paper);
    color: var(--ink);
}

.block-container {
    max-width: 1220px;
    padding-top: 1.3rem;
    padding-bottom: 3rem;
}

header[data-testid="stHeader"] {
    background: transparent;
}

#MainMenu, footer {
    visibility: hidden;
}

.brand {
    font-family: 'Manrope', sans-serif;
    font-size: 1.55rem;
    font-weight: 800;
    color: var(--navy);
    letter-spacing: -0.8px;
}

.brand span {
    color: var(--teal);
}

.brand-subtitle {
    color: var(--muted);
    font-size: 0.76rem;
    margin-top: -3px;
    letter-spacing: 0.3px;
}

.topline {
    border-bottom: 1px solid var(--line);
    margin: 12px 0 22px 0;
}

.eyebrow {
    color: var(--teal);
    text-transform: uppercase;
    letter-spacing: 2px;
    font-size: 0.74rem;
    font-weight: 700;
    margin-bottom: 10px;
}

.hero {
    background: linear-gradient(125deg, #10243A 0%, #173D54 65%, #087E8B 120%);
    padding: 52px 48px;
    border-radius: 24px;
    color: white;
    margin: 12px 0 25px 0;
    border: 1px solid rgba(255,255,255,0.1);
}

.hero h1 {
    color: white;
    font-family: 'Manrope', sans-serif;
    font-size: clamp(2.2rem, 5vw, 3.8rem);
    font-weight: 800;
    line-height: 1.12;
    letter-spacing: -1.8px;
    margin: 10px 0 18px 0;
}

.hero p {
    color: #E1EDF3;
    font-size: 1.05rem;
    line-height: 1.8;
    max-width: 760px;
}

.hero-tag {
    display: inline-block;
    border: 1px solid rgba(255,255,255,0.35);
    border-radius: 30px;
    padding: 7px 13px;
    color: #DDF7F6;
    font-size: 0.78rem;
    letter-spacing: 0.6px;
}

.section-title {
    font-family: 'Manrope', sans-serif;
    font-size: 1.65rem;
    font-weight: 800;
    color: var(--navy);
    letter-spacing: -0.6px;
    margin-bottom: 7px;
}

.section-desc {
    color: var(--muted);
    line-height: 1.75;
    margin-bottom: 22px;
}

.info-card {
    background: white;
    border: 1px solid var(--line);
    border-radius: 16px;
    padding: 22px;
    height: 100%;
}

.info-card h3 {
    color: var(--navy);
    font-family: 'Manrope', sans-serif;
    font-size: 1.05rem;
    margin: 8px 0;
}

.info-card p {
    color: var(--muted);
    font-size: 0.92rem;
    line-height: 1.7;
}

.card-icon {
    font-size: 1.6rem;
}

.metric-card {
    background: white;
    border: 1px solid var(--line);
    border-radius: 15px;
    padding: 20px;
}

.metric-label {
    color: var(--muted);
    font-size: 0.81rem;
    font-weight: 600;
}

.metric-value {
    font-family: 'Manrope', sans-serif;
    color: var(--navy);
    font-size: 2rem;
    font-weight: 800;
    margin: 6px 0 3px 0;
}

.metric-note {
    color: var(--teal);
    font-size: 0.77rem;
}

.note-box {
    background: var(--teal-light);
    border-left: 4px solid var(--teal);
    padding: 15px 18px;
    border-radius: 7px;
    color: #244C58;
    line-height: 1.7;
    margin: 14px 0;
}

.warning-box {
    background: #FFF7E8;
    border-left: 4px solid #C58A27;
    padding: 15px 18px;
    border-radius: 7px;
    color: #664C1C;
    line-height: 1.7;
    margin: 14px 0;
}

.page-footer {
    border-top: 1px solid var(--line);
    margin-top: 42px;
    padding-top: 18px;
    color: var(--muted);
    font-size: 0.82rem;
}

.stButton > button {
    border-radius: 10px;
    font-weight: 600;
    border: 1px solid var(--line);
    min-height: 42px;
    transition: all 0.2s ease;
}

.stButton > button[kind="primary"] {
    background: var(--teal);
    border-color: var(--teal);
    color: white;
}

.stButton > button:hover {
    border-color: var(--teal);
    color: var(--teal);
}

div[data-testid="stDataFrame"] {
    border: 1px solid var(--line);
    border-radius: 12px;
    overflow: hidden;
}

.stTabs [data-baseweb="tab-list"] {
    gap: 12px;
}

.stTabs [data-baseweb="tab"] {
    font-weight: 600;
}

@media (max-width: 700px) {
    .block-container {
        padding-left: 1rem;
        padding-right: 1rem;
    }
    .hero {
        padding: 30px 23px;
        border-radius: 17px;
    }
    .hero h1 {
        font-size: 2.2rem;
    }
    .info-card {
        padding: 17px;
    }
}
</style>
""", unsafe_allow_html=True)


# --------------------------------------------------
# HEADER
# --------------------------------------------------
header_left, header_right = st.columns([3, 1])

with header_left:
    st.markdown("""
    <div class="brand">Bac<span>Trace</span> 🧬</div>
    <div class="brand-subtitle">COMPARATIVE MICROBIAL SEQUENCE ANALYSIS</div>
    """, unsafe_allow_html=True)

with header_right:
    st.markdown(
        f"<div style='text-align:right;color:#64748B;font-size:0.82rem;padding-top:8px;'>"
        f"RESEARCH PROJECT<br><b style='color:#10243A;'>Page {st.session_state.page_index + 1} of {len(PAGES)}</b></div>",
        unsafe_allow_html=True
    )

st.markdown("<div class='topline'></div>", unsafe_allow_html=True)

# Horizontal page navigation
nav_cols = st.columns(7)
for i, label in enumerate(PAGES):
    with nav_cols[i]:
        st.button(
            label,
            key=f"top_nav_{i}",
            on_click=go_to_page,
            args=(i,),
            type="primary" if st.session_state.page_index == i else "secondary",
            use_container_width=True
        )

st.progress((st.session_state.page_index + 1) / len(PAGES))


# --------------------------------------------------
# REUSABLE DISPLAY HELPERS
# --------------------------------------------------
def section_heading(title, description):
    st.markdown(f"""
    <div class="eyebrow">BacTrace Research Workspace</div>
    <div class="section-title">{title}</div>
    <div class="section-desc">{description}</div>
    """, unsafe_allow_html=True)


def metric_card(label, value, note):
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-label">{label}</div>
        <div class="metric-value">{value}</div>
        <div class="metric-note">{note}</div>
    </div>
    """, unsafe_allow_html=True)


def info_card(icon, title, description):
    st.markdown(f"""
    <div class="info-card">
        <div class="card-icon">{icon}</div>
        <h3>{title}</h3>
        <p>{description}</p>
    </div>
    """, unsafe_allow_html=True)


# --------------------------------------------------
# PAGE 1 — HOME
# --------------------------------------------------
page = st.session_state.page_index

if page == 0:
    st.markdown("""
    <div class="hero">
        <div class="hero-tag">BIOINFORMATICS • 16S rRNA • COMPARATIVE ANALYSIS</div>
        <h1>Understanding bacterial similarity,<br>one sequence at a time.</h1>
        <p>
        BacTrace is a comparative sequence analysis platform designed to explore
        16S ribosomal RNA sequences from three closely related bacterial species:
        <i>Klebsiella pneumoniae</i>, <i>Klebsiella variicola</i>, and
        <i>Klebsiella quasipneumoniae</i>.
        </p>
        <p>
        Explore the dataset, conserved and variable positions, pairwise sequence
        identity, and the role of phylogenetic analysis in bacterial identification.
        </p>
    </div>
    """, unsafe_allow_html=True)

    start_col, empty_col = st.columns([1, 2])
    with start_col:
        st.button(
            "Start Exploring  →",
            key="start_exploring",
            type="primary",
            use_container_width=True,
            on_click=go_to_page,
            args=(1,)
        )

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("### Inside the platform")

    cards = st.columns(3)
    with cards[0]:
        info_card("🧬", "Sequence dataset",
                  "Explore 11 selected 16S rRNA accession records across three Klebsiella species.")
    with cards[1]:
        info_card("📊", "Conservation & similarity",
                  "Review alignment statistics, conserved positions, variable sites, and pairwise identity.")
    with cards[2]:
        info_card("🌿", "Phylogenetic interpretation",
                  "Use evolutionary relationships alongside sequence evidence, while considering limitations.")

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("""
    <div class="note-box">
    <b>Research focus:</b> Closely related Klebsiella species can show high 16S rRNA
    sequence similarity. BacTrace examines how useful this marker is for comparison
    and why species-level identification may require additional evidence.
    </div>
    """, unsafe_allow_html=True)


# --------------------------------------------------
# PAGE 2 — RESEARCH OVERVIEW
# --------------------------------------------------
elif page == 1:
    section_heading(
        "Research overview",
        "The scientific question, study design, and objectives behind BacTrace."
    )

    st.markdown("### Background")
    st.write(
        "The genus Klebsiella includes closely related bacterial species whose "
        "identification can be challenging when relying on a single molecular marker. "
        "The 16S ribosomal RNA gene is widely used in bacterial sequence analysis, "
        "but closely related species may share highly similar sequences."
    )

    st.markdown("### Research question")
    st.markdown("""
    <div class="note-box">
    How effectively can comparative 16S rRNA sequence analysis describe the
    similarities and differences among <i>K. pneumoniae</i>,
    <i>K. variicola</i>, and <i>K. quasipneumoniae</i>?
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### Aim")
    st.write(
        "To compare selected 16S rRNA sequences from three closely related "
        "Klebsiella species using multiple sequence alignment, conservation "
        "analysis, pairwise identity, and phylogenetic analysis."
    )

    st.markdown("### Objectives")
    objectives = [
        "Collect selected 16S rRNA reference or database sequences from NCBI.",
        "Align the sequences and examine alignment length and sequence coverage.",
        "Identify fully conserved positions, variable sites, and parsimony-informative sites.",
        "Compare within-species and between-species pairwise sequence identity.",
        "Interpret phylogenetic relationships and explain the limitations of 16S rRNA identification."
    ]
    for objective in objectives:
        st.markdown(f"- {objective}")

    st.markdown("### Species included")
    species_cols = st.columns(3)
    with species_cols[0]:
        info_card("01", "Klebsiella pneumoniae",
                  "A clinically important species included for comparative sequence analysis.")
    with species_cols[1]:
        info_card("02", "Klebsiella variicola",
                  "A closely related species that may be difficult to distinguish using limited marker information.")
    with species_cols[2]:
        info_card("03", "Klebsiella quasipneumoniae",
                  "A related species included to examine sequence similarity across the group.")

    st.markdown("""
    <div class="warning-box">
    <b>Interpretation boundary:</b> This is an in-silico sequence comparison study.
    Similarity values alone do not establish clinical identity, pathogenicity,
    antibiotic resistance, or a definitive species diagnosis.
    </div>
    """, unsafe_allow_html=True)


# --------------------------------------------------
# PAGE 3 — DATASET
# --------------------------------------------------
elif page == 2:
    section_heading(
        "Sequence dataset",
        "Browse the selected accession records used in the comparative analysis."
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        metric_card("Selected sequences", "11", "Across three species")
    with c2:
        metric_card("Species represented", "3", "Closely related Klebsiella")
    with c3:
        metric_card("Alignment length", "1,530", "Multiple sequence alignment columns")

    st.markdown("### Accession browser")
    species_filter = st.selectbox(
        "Filter by species",
        ["All species"] + sorted(accession_df["Species"].unique().tolist())
    )

    filtered_df = accession_df.copy()
    if species_filter != "All species":
        filtered_df = filtered_df[filtered_df["Species"] == species_filter]

    st.dataframe(
        filtered_df,
        use_container_width=True,
        hide_index=True
    )

    csv_data = filtered_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download accession table (CSV)",
        data=csv_data,
        file_name="BacTrace_accessions.csv",
        mime="text/csv"
    )

    st.markdown("### Database source")
    st.write(
        "The accession identifiers refer to sequence records selected from NCBI. "
        "Before final submission, confirm each record's organism name, sequence "
        "type, completeness, and record status against the original NCBI entry."
    )

    st.markdown("""
    <div class="warning-box">
    <b>Sequence coverage matters:</b> The alignment includes sequences of different
    lengths. KR092086.1 is a shorter partial sequence, and some records may contain
    gaps or ambiguous bases. Pairwise comparisons should therefore report the number
    of positions actually compared, not identity percentage alone.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### Open NCBI records")
    st.markdown(
        "[Search NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) · "
        "[NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)"
    )


# --------------------------------------------------
# PAGE 4 — CONSERVATION ANALYSIS
# --------------------------------------------------
elif page == 3:
    section_heading(
        "Conservation analysis",
        "A summary of conserved and variable positions in the 1,530-column alignment."
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        metric_card("Alignment length", "1,530", "Alignment columns")
    with c2:
        metric_card("Fully conserved", "714", "46.67% of alignment")
    with c3:
        metric_card("Variable sites", "43", "2.81% of alignment")

    st.markdown("### Reported site statistics")

    conservation_df = pd.DataFrame([
        {"Measure": "Fully conserved sites", "Count": 714, "Percentage": 46.67},
        {"Measure": "Variable sites", "Count": 43, "Percentage": 2.81},
        {"Measure": "Parsimony-informative sites", "Count": 15, "Percentage": 0.98},
    ])

    fig = px.bar(
        conservation_df,
        x="Measure",
        y="Count",
        text="Count",
        color="Measure",
        color_discrete_sequence=["#087E8B", "#E3A44D", "#52749A"],
        title="Reported alignment-site statistics"
    )
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="DM Sans", color="#203449"),
        showlegend=False,
        xaxis_title="Site category",
        yaxis_title="Number of sites",
        margin=dict(l=20, r=20, t=60, b=20)
    )
    fig.update_traces(textposition="outside")
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("### How to interpret these measures")
    tab1, tab2, tab3 = st.tabs([
        "Conserved sites", "Variable sites", "Parsimony-informative"
    ])

    with tab1:
        st.write(
            "A fully conserved site was defined as a column containing the same "
            "unambiguous A, C, G, or T nucleotide in every sequence, with no gaps "
            "or ambiguous characters. The analysis reported 714 such sites."
        )
    with tab2:
        st.write(
            "A variable site was defined as a column containing at least two "
            "different unambiguous nucleotide states among the sequences. "
            "Gaps and ambiguous characters were ignored when counting states."
        )
    with tab3:
        st.write(
            "A parsimony-informative site has at least two nucleotide states, "
            "each represented by at least two sequences. The 15 reported "
            "parsimony-informative sites are a subset of variable sites."
        )

    st.markdown("""
    <div class="note-box">
    <b>Important:</b> These categories overlap: parsimony-informative sites are
    included among variable sites. Do not add all three counts together. The 714
    fully conserved sites represent 46.67% of the 1,530 alignment columns under
    the strict definition used in this analysis.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### Alignment tool")
    st.write(
        "The sequences were aligned using MAFFT through SeqBench. The conservation "
        "statistics displayed here come from the separately calculated alignment "
        "analysis. They may differ from the alignment tool's overall conservation "
        "summary because the tools can use different definitions."
    )


# --------------------------------------------------
# PAGE 5 — SIMILARITY ANALYSIS
# --------------------------------------------------
elif page == 4:
    section_heading(
        "Similarity explorer",
        "Compare reported pairwise identity within species and between species."
    )

    st.markdown("### Main findings")

    c1, c2, c3 = st.columns(3)
    with c1:
        metric_card("Mean within-species identity", "99.54%", "15 pairwise comparisons")
    with c2:
        metric_card("Mean between-species identity", "99.07%", "40 pairwise comparisons")
    with c3:
        metric_card("Mean difference", "0.47 pp", "Within minus between")

    st.markdown("### Within-species vs between-species identity")

    chart_df = similarity_summary.melt(
        id_vars=["Comparison"],
        value_vars=[
            "Mean identity (%)",
            "Minimum identity (%)",
            "Maximum identity (%)"
        ],
        var_name="Statistic",
        value_name="Identity (%)"
    )

    fig = px.bar(
        chart_df,
        x="Comparison",
        y="Identity (%)",
        color="Statistic",
        barmode="group",
        title="Pairwise identity summary",
        color_discrete_sequence=["#087E8B", "#E3A44D", "#52749A"],
        hover_data={"Identity (%)": ":.4f"}
    )
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="DM Sans", color="#203449"),
        yaxis=dict(title="Identity (%)", range=[96.5, 100.3]),
        xaxis_title="Comparison group",
        margin=dict(l=20, r=20, t=60, b=20),
        legend_title_text=""
    )
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("### Summary table")
    display_summary = similarity_summary.copy()
    display_summary["Mean identity (%)"] = display_summary["Mean identity (%)"].map(lambda x: f"{x:.4f}")
    display_summary["Minimum identity (%)"] = display_summary["Minimum identity (%)"].map(lambda x: f"{x:.4f}")
    display_summary["Maximum identity (%)"] = display_summary["Maximum identity (%)"].map(lambda x: f"{x:.4f}")
    display_summary["Mean comparable positions"] = display_summary["Mean comparable positions"].map(lambda x: f"{x:.2f}")

    st.dataframe(display_summary, use_container_width=True, hide_index=True)

    st.markdown("### Selected pairwise comparisons")

    st.markdown("""
    <div class="info-card">
        <h3>Highest reported identity</h3>
        <p><b>NR_117686.1 vs AB004753.1</b><br>
        100% identity across 1,448 comparable unambiguous nucleotide positions.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("""
    <div class="info-card">
        <h3>Lowest reported identity</h3>
        <p><b>KR092086.1 vs OM780282.1</b><br>
        Approximately 97.29% identity across 739 comparable unambiguous nucleotide positions.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="warning-box">
    <b>Interpretation:</b> Identity was calculated over shared unambiguous A/C/G/T
    positions, excluding gaps and ambiguous bases. The values are descriptive and
    do not prove that 16S rRNA can reliably distinguish these species. The pairwise
    coverage and sequence metadata should be checked before making a diagnostic claim.
    </div>
    """, unsafe_allow_html=True)

    st.download_button(
        "Download similarity summary (CSV)",
        data=similarity_summary.to_csv(index=False).encode("utf-8"),
        file_name="BacTrace_similarity_summary.csv",
        mime="text/csv"
    )


# --------------------------------------------------
# PAGE 6 — PHYLOGENETIC ANALYSIS
# --------------------------------------------------
elif page == 5:
    section_heading(
        "Phylogenetic analysis",
        "Explore evolutionary relationships after the IQ-TREE analysis is finalized."
    )

    st.markdown("""
    <div class="hero">
        <div class="hero-tag">EVOLUTIONARY COMPARISON</div>
        <h1 style="font-size:2.5rem;">From alignment to tree</h1>
        <p>
        A phylogenetic tree summarizes the relationships inferred from the aligned
        sequences. Its interpretation should consider branch lengths, support values,
        sequence coverage, and the limitations of the selected marker.
        </p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### Tree status")
    st.markdown("""
    <div class="warning-box">
    <b>Tree image and topology pending.</b> The final IQ-TREE tree has not yet been
    added to this website. No species-clustering result or bootstrap-support value
    is assumed here.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### What to examine in the final tree")
    tree_cols = st.columns(3)
    with tree_cols[0]:
        info_card("01", "Clustering",
                  "Check whether sequences assigned to the same species cluster together.")
    with tree_cols[1]:
        info_card("02", "Branch lengths",
                  "Review the estimated amount of sequence change represented by branches.")
    with tree_cols[2]:
        info_card("03", "Branch support",
                  "Use support values to assess confidence in individual inferred groupings.")

    st.markdown("### Add your final tree")
    st.write(
        "When your IQ-TREE result is ready, upload a tree image or a supported "
        "tree file to your project repository. The tree can then be embedded here "
        "and discussed using the actual topology and support values."
    )

    st.markdown("Useful tools:")
    st.markdown(
        "- [IQ-TREE](http://www.iqtree.org/)\n"
        "- [iTOL — Interactive Tree Of Life](https://itol.embl.de/)"
    )

    st.markdown("""
    <div class="note-box">
    <b>Scientific caution:</b> A tree inferred from 16S rRNA sequences may not
    resolve closely related Klebsiella species reliably. Avoid interpreting a
    branch as definitive species identification without additional supporting evidence.
    </div>
    """, unsafe_allow_html=True)


# --------------------------------------------------
# PAGE 7 — METHODS AND REFERENCES
# --------------------------------------------------
elif page == 6:
    section_heading(
        "Methods & references",
        "The computational workflow, analysis definitions, and resources supporting the study."
    )

    st.markdown("### Workflow")

    workflow = [
        ("01", "Sequence retrieval", "Selected 16S rRNA records from NCBI Nucleotide."),
        ("02", "Multiple sequence alignment", "Aligned sequences using MAFFT through SeqBench."),
        ("03", "Alignment statistics", "Summarized alignment length, gaps, and ambiguous bases."),
        ("04", "Conservation analysis", "Counted fully conserved, variable, and parsimony-informative positions."),
        ("05", "Pairwise identity", "Compared shared unambiguous nucleotide positions."),
        ("06", "Phylogenetic inference", "Prepared an IQ-TREE analysis; final tree interpretation is pending.")
    ]

    for number, title, description in workflow:
        st.markdown(f"""
        <div class="info-card" style="margin-bottom:10px;">
            <div style="color:#087E8B;font-weight:800;font-size:0.8rem;">STEP {number}</div>
            <h3>{title}</h3>
            <p>{description}</p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("### Analysis definitions")
    st.write(
        "**Pairwise identity:** Number of matching unambiguous nucleotide positions "
        "divided by the number of comparable unambiguous nucleotide positions, "
        "multiplied by 100."
    )
    st.write(
        "**Fully conserved site:** A column with the same unambiguous A/C/G/T base "
        "in every sequence and no gaps or ambiguous characters."
    )
    st.write(
        "**Variable site:** A column with at least two distinct unambiguous nucleotide "
        "states, ignoring gaps and ambiguous characters."
    )
    st.write(
        "**Parsimony-informative site:** A column with at least two nucleotide states, "
        "each occurring in at least two sequences."
    )

    st.markdown("### Tools and resources")
    references = [
        (
            "NCBI Nucleotide",
            "Sequence retrieval and accession verification.",
            "https://www.ncbi.nlm.nih.gov/nuccore/"
        ),
        (
            "MAFFT",
            "Katoh, K. & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7.",
            "https://doi.org/10.1093/molbev/mst010"
        ),
        (
            "Biopython",
            "Cock, P. J. A. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics.",
            "https://doi.org/10.1093/bioinformatics/btp163"
        ),
        (
            "IQ-TREE 2",
            "Minh, B. Q. et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era.",
            "https://doi.org/10.1093/molbev/msaa015"
        ),
        (
            "Klebsiella taxonomy",
            "Boye, K. & Hansen, D. S. (2003). Sequencing of 16S rDNA of Klebsiella: taxonomic relations within the genus and to other Enterobacteriaceae.",
            "https://pubmed.ncbi.nlm.nih.gov/12635932/"
        ),
        (
            "Klebsiella species identification",
            "Study discussing misidentification of K. pneumoniae, K. variicola, and K. quasipneumoniae.",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5541162/"
        )
    ]

    for title, description, url in references:
        st.markdown(f"**{title}**  \n{description}  \n[Open resource]({url})")
        st.markdown("---")

    st.markdown("### Limitations")
    limitations = [
        "The dataset contains only 11 selected sequences and may not represent the full diversity of each species.",
        "Sequence lengths and coverage differ, affecting which positions can be compared.",
        "Accession metadata and record completeness should be verified against NCBI before final reporting.",
        "16S rRNA sequence identity alone may not reliably separate closely related Klebsiella species.",
        "The current website reports descriptive comparisons and does not establish diagnostic performance."
    ]
    for item in limitations:
        st.markdown(f"- {item}")

    st.markdown("""
    <div class="note-box">
    <b>Next research improvement:</b> Validate accession metadata, add the finalized
    phylogenetic tree and support values, and consider complementary molecular markers
    or whole-genome evidence for improved species-level discrimination.
    </div>
    """, unsafe_allow_html=True)


# --------------------------------------------------
# FOOTER + PREVIOUS / NEXT
# --------------------------------------------------
st.markdown("""
<div class="page-footer">
    <b style="color:#10243A;">BacTrace</b> · Comparative 16S rRNA Sequence Analysis
    <br>
    Educational bioinformatics project · Findings should be interpreted within the
    limitations of the dataset and selected molecular marker.
</div>
""", unsafe_allow_html=True)

st.markdown("<br>", unsafe_allow_html=True)

prev_col, page_col, next_col = st.columns([1, 2, 1])

with prev_col:
    st.button(
        "← Previous",
        key="previous_page",
        use_container_width=True,
        disabled=(page == 0),
        on_click=go_to_page,
        args=(max(0, page - 1),)
    )

with page_col:
    st.markdown(
        f"<div style='text-align:center;color:#64748B;padding-top:10px;'>"
        f"<b style='color:#10243A;'>{PAGES[page]}</b>"
        f" &nbsp; · &nbsp; {page + 1} / {len(PAGES)}"
        f"</div>",
        unsafe_allow_html=True
    )

with next_col:
    st.button(
        "Next →",
        key="next_page",
        use_container_width=True,
        disabled=(page == len(PAGES) - 1),
        type="primary",
        on_click=go_to_page,
        args=(min(len(PAGES) - 1, page + 1),)
    )
