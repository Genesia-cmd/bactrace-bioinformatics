import io
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from itertools import combinations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

st.set_page_config(
    page_title="BacTrace | Comparative Sequence Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
.stApp { background: linear-gradient(180deg,#f5f9fb 0%,#ffffff 44%); }
[data-testid="stSidebar"] { background:#102b3f; }
[data-testid="stSidebar"] * { color:#f4f8fb !important; }
.hero { padding:2rem 2.2rem; border-radius:20px; background:linear-gradient(120deg,#102b3f,#176b75); color:white; margin-bottom:1.3rem; }
.hero h1 { color:white; margin:0; font-size:2.3rem; }
.hero p { color:#e0f2f1; font-size:1.02rem; margin:.55rem 0 0; }
.panel { padding:1rem 1.2rem; border:1px solid #dbe5ec; border-radius:14px; background:white; margin:.4rem 0 1rem; }
.small-note { color:#526575; font-size:.9rem; }
div[data-testid="stMetric"] { background:white; border:1px solid #dbe5ec; padding:.8rem; border-radius:12px; }
.stButton button,.stDownloadButton button { border-radius:9px; font-weight:600; }
</style>
""", unsafe_allow_html=True)

IUPAC = set("ACGTRYSWKMBDHVN")
UNAMBIGUOUS = set("ACGT")
PAGES = [
    "Home & Input",
    "Analysis Overview",
    "Sequence Alignment",
    "Conservation Analysis",
    "Similarity Explorer",
    "Phylogenetic Analysis",
    "Report & Export",
]

if "analysis" not in st.session_state:
    st.session_state.analysis = None


# ------------------------------- Input ----------------------------------------
def parse_fasta_text(text):
    """Parse multi-FASTA, validate unique IDs and IUPAC DNA characters."""
    text = (text or "").strip()
    if not text:
        raise ValueError("Upload a FASTA file or paste at least two FASTA records.")
    if not text.startswith(">"):
        raise ValueError("FASTA input must begin with a header line starting with '>'.")
    records, seen = [], set()
    current_id, current_desc, chunks = None, "", []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, current_desc, "".join(chunks).upper()))
            header = line[1:].strip()
            if not header:
                raise ValueError("A FASTA header is empty.")
            parts = header.split(maxsplit=1)
            current_id = parts[0]
            current_desc = parts[1] if len(parts) > 1 else ""
            if current_id in seen:
                raise ValueError(f"Duplicate FASTA ID '{current_id}'. Use unique IDs.")
            seen.add(current_id)
            chunks = []
        else:
            if current_id is None:
                raise ValueError("Sequence data was found before the first FASTA header.")
            chunks.append(re.sub(r"\s+", "", line))
    if current_id is not None:
        records.append((current_id, current_desc, "".join(chunks).upper()))
    if len(records) < 2:
        raise ValueError("At least two sequences are required for comparative analysis.")
    output, issues = [], []
    for rid, desc, seq in records:
        if not seq:
            issues.append(f"{rid}: empty sequence")
            continue
        invalid = sorted(set(seq) - IUPAC)
        if invalid:
            issues.append(f"{rid}: unsupported character(s) {', '.join(invalid)}. Use DNA IUPAC symbols only.")
            continue
        output.append(SeqRecord(Seq(seq), id=rid, description=desc))
    if issues:
        raise ValueError("Fix the following input issue(s):\n- " + "\n- ".join(issues))
    return output


def get_species(rec):
    text = f"{rec.id} {rec.description}".lower()
    patterns = [
        ("K. pneumoniae", ["klebsiella pneumoniae", "k. pneumoniae", "k_pneumoniae", "k-pneumoniae"]),
        ("K. variicola", ["klebsiella variicola", "k. variicola", "k_variicola", "k-variicola"]),
        ("K. quasipneumoniae", ["klebsiella quasipneumoniae", "k. quasipneumoniae", "k_quasipneumoniae", "k-quasipneumoniae"]),
    ]
    for label, terms in patterns:
        if any(term in text for term in terms):
            return label
    return "Unassigned"


def fasta_string(records):
    handle = io.StringIO()
    SeqIO.write(records, handle, "fasta")
    return handle.getvalue()


# ------------------------------- Alignment ------------------------------------
def star_alignment_fallback(records):
    """
    Reference-based star alignment with consistent insertion slots.
    This is a documented fallback, not a replacement for a full MSA program.
    """
    ref = str(records[0].seq).upper()
    aligner = PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    per_record = []
    max_insertions = [0] * (len(ref) + 1)

    for rec in records:
        seq = str(rec.seq).upper()
        insertions = [""] * (len(ref) + 1)
        mapped = ["-"] * len(ref)
        if rec.id == records[0].id:
            mapped = list(ref)
        else:
            alignment = aligner.align(ref, seq)[0]
            coords = alignment.coordinates
            for j in range(coords.shape[1] - 1):
                r0, r1 = int(coords[0, j]), int(coords[0, j + 1])
                s0, s1 = int(coords[1, j]), int(coords[1, j + 1])
                dr, ds = r1 - r0, s1 - s0
                if dr > 0 and ds > 0:
                    if dr != ds:
                        raise ValueError("Unexpected pairwise alignment segment; use MAFFT for this input.")
                    for offset in range(dr):
                        mapped[r0 + offset] = seq[s0 + offset]
                elif dr > 0 and ds == 0:
                    for pos in range(r0, r1):
                        mapped[pos] = "-"
                elif dr == 0 and ds > 0:
                    insertions[r0] += seq[s0:s1]
        for boundary, ins in enumerate(insertions):
            max_insertions[boundary] = max(max_insertions[boundary], len(ins))
        per_record.append((rec.id, insertions, mapped))

    aligned_records = []
    for rec, (rid, insertions, mapped) in zip(records, per_record):
        pieces = []
        for boundary in range(len(ref) + 1):
            ins = insertions[boundary]
            pieces.append(ins + "-" * (max_insertions[boundary] - len(ins)))
            if boundary < len(ref):
                pieces.append(mapped[boundary])
        aligned_records.append(SeqRecord(Seq("".join(pieces)), id=rid, description=""))
    return MultipleSeqAlignment(aligned_records)


def align_sequences(records):
    mafft = shutil.which("mafft")
    if mafft:
        with tempfile.TemporaryDirectory() as td:
            inp = os.path.join(td, "input.fasta")
            SeqIO.write(records, inp, "fasta")
            try:
                proc = subprocess.run(
                    [mafft, "--auto", "--quiet", inp],
                    capture_output=True, text=True, timeout=240
                )
                if proc.returncode == 0 and proc.stdout.strip():
                    alignment = AlignIO.read(io.StringIO(proc.stdout), "fasta")
                    if len(alignment) == len(records):
                        return alignment, "MAFFT (--auto)"
            except Exception:
                pass
    return star_alignment_fallback(records), "Reference-based Biopython star-alignment fallback (not equivalent to MAFFT)"


# ------------------------------- Analysis -------------------------------------
def sequence_statistics(records):
    rows = []
    for rec in records:
        seq = str(rec.seq).upper()
        n = len(seq)
        rows.append({
            "Sequence ID": rec.id,
            "Description": rec.description,
            "Assigned group": get_species(rec),
            "Length (nt)": n,
            "A": seq.count("A"),
            "C": seq.count("C"),
            "G": seq.count("G"),
            "T": seq.count("T"),
            "Ambiguous bases": sum(c not in UNAMBIGUOUS for c in seq),
            "GC content (%)": round(100 * (seq.count("G") + seq.count("C")) / n, 3) if n else np.nan,
        })
    return pd.DataFrame(rows)


def conservation_analysis(alignment):
    rows = []
    for pos in range(alignment.get_alignment_length()):
        chars = [str(rec.seq[pos]).upper() for rec in alignment]
        bases = [c for c in chars if c in UNAMBIGUOUS]
        counts = pd.Series(bases).value_counts().to_dict() if bases else {}
        states = sorted(counts)
        fully_conserved = all(c in UNAMBIGUOUS for c in chars) and len(set(chars)) == 1
        variable = len(states) >= 2
        informative = sum(v >= 2 for v in counts.values()) >= 2
        rows.append({
            "Alignment position": pos + 1,
            "Nucleotide states": ", ".join(states) if states else "No unambiguous bases",
            "Distinct states": len(states),
            "Gap count": chars.count("-"),
            "Ambiguous count": sum(c not in UNAMBIGUOUS and c != "-" for c in chars),
            "Fully conserved": fully_conserved,
            "Variable": variable,
            "Parsimony-informative": informative,
        })
    df = pd.DataFrame(rows)
    summary = {
        "Alignment columns": alignment.get_alignment_length(),
        "Fully conserved sites": int(df["Fully conserved"].sum()),
        "Variable sites": int(df["Variable"].sum()),
        "Parsimony-informative sites": int(df["Parsimony-informative"].sum()),
    }
    return df, summary


def pairwise_table(alignment, records):
    seqs = {rec.id: str(rec.seq).upper() for rec in alignment}
    groups = {rec.id: get_species(rec) for rec in records}
    rows = []
    for id1, id2 in combinations(seqs, 2):
        a, b = seqs[id1], seqs[id2]
        comparable = matches = 0
        for x, y in zip(a, b):
            if x in UNAMBIGUOUS and y in UNAMBIGUOUS:
                comparable += 1
                matches += (x == y)
        identity = (100 * matches / comparable) if comparable else np.nan
        g1, g2 = groups[id1], groups[id2]
        if "Unassigned" in (g1, g2):
            kind = "Group unassigned"
        elif g1 == g2:
            kind = "Within-species"
        else:
            kind = "Between-species"
        rows.append({
            "Sequence 1": id1,
            "Sequence 2": id2,
            "Identity (%)": round(identity, 5) if pd.notna(identity) else np.nan,
            "Matching positions": int(matches),
            "Comparable positions": int(comparable),
            "Group 1": g1,
            "Group 2": g2,
            "Comparison type": kind,
        })
    return pd.DataFrame(rows)


def make_tree(alignment):
    # Biopython's identity distance with neighbour joining; no bootstrap support.
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    tree = DistanceTreeConstructor().nj(distance_matrix)
    return tree, "Neighbour-joining using Biopython identity distances (no bootstrap support)"


def newick_string(tree):
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue()


def tree_figure(tree):
    terminals = tree.get_terminals()
    y = {term.name: i for i, term in enumerate(terminals)}
    x = {}
    def set_x(clade, value):
        x[clade] = value
        for child in clade.clades:
            set_x(child, value + max(child.branch_length or 0, 0))
    set_x(tree.root, 0)
    fig = go.Figure()
    for parent in tree.find_clades(order="preorder"):
        if parent.clades:
            py = (min(y[t.name] for t in parent.get_terminals()) + max(y[t.name] for t in parent.get_terminals())) / 2
            for child in parent.clades:
                cy = (min(y[t.name] for t in child.get_terminals()) + max(y[t.name] for t in child.get_terminals())) / 2
                fig.add_trace(go.Scatter(x=[x[parent], x[parent], x[child]], y=[py, cy, cy],
                                         mode="lines", line=dict(width=2), hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(
        x=[x[t] for t in terminals], y=[y[t.name] for t in terminals],
        text=[t.name for t in terminals], mode="markers+text",
        textposition="middle right", marker=dict(size=7), showlegend=False
    ))
    fig.update_layout(
        title="Neighbour-joining tree",
        xaxis_title="Cumulative distance",
        yaxis=dict(showticklabels=False, title=""),
        height=max(430, len(terminals) * 38),
        margin=dict(l=15, r=180, t=60, b=40),
    )
    return fig


def build_analysis(records):
    alignment, alignment_method = align_sequences(records)
    conservation, cons_summary = conservation_analysis(alignment)
    pairs = pairwise_table(alignment, records)
    tree, tree_method = make_tree(alignment)
    return {
        "records": records,
        "stats": sequence_statistics(records),
        "alignment": alignment,
        "alignment_method": alignment_method,
        "aligned_fasta": fasta_string(list(alignment)),
        "conservation": conservation,
        "cons_summary": cons_summary,
        "pairs": pairs,
        "tree": tree,
        "tree_method": tree_method,
        "newick": newick_string(tree),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }


def create_pdf(analysis):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle

    output = io.BytesIO()
    doc = SimpleDocTemplate(output, pagesize=A4, rightMargin=1.6*cm, leftMargin=1.6*cm,
                            topMargin=1.5*cm, bottomMargin=1.5*cm)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="SectionBlue", parent=styles["Heading2"],
                              textColor=colors.HexColor("#176b75"), spaceBefore=10))
    s = analysis["cons_summary"]
    n = s["Alignment columns"]
    story = [
        Paragraph("BacTrace", styles["Title"]),
        Paragraph("Comparative 16S rRNA Sequence Analysis Report", styles["Heading2"]),
        Paragraph(f"Generated {analysis['timestamp']}", styles["Normal"]),
        Spacer(1, 12),
        Paragraph("1. Dataset and methods", styles["SectionBlue"]),
        Paragraph(
            f"Sequences analysed: {len(analysis['records'])}<br/>"
            f"Alignment columns: {n}<br/>"
            f"Alignment method: {analysis['alignment_method']}<br/>"
            f"Tree method: {analysis['tree_method']}",
            styles["BodyText"]),
        Paragraph("2. Conservation analysis", styles["SectionBlue"]),
    ]
    table_data = [["Metric", "Count", "Percentage"]]
    for label, key in [
        ("Fully conserved sites", "Fully conserved sites"),
        ("Variable sites", "Variable sites"),
        ("Parsimony-informative sites", "Parsimony-informative sites"),
    ]:
        table_data.append([label, str(s[key]), f"{100*s[key]/n:.2f}%" if n else "N/A"])
    table = Table(table_data, colWidths=[8.5*cm, 2.5*cm, 3*cm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#102b3f")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("GRID", (0, 0), (-1, -1), .4, colors.HexColor("#ccd6dd")),
        ("PADDING", (0, 0), (-1, -1), 6),
    ]))
    story += [table, Paragraph("3. Pairwise identity", styles["SectionBlue"])]
    pairs = analysis["pairs"].dropna(subset=["Identity (%)"])
    if len(pairs):
        story.append(Paragraph(
            f"Unique pairs: {len(pairs)}<br/>Mean identity: {pairs['Identity (%)'].mean():.3f}%<br/>"
            f"Minimum identity: {pairs['Identity (%)'].min():.3f}%<br/>"
            f"Maximum identity: {pairs['Identity (%)'].max():.3f}%",
            styles["BodyText"]))
    else:
        story.append(Paragraph("No pair had comparable unambiguous nucleotide positions.", styles["BodyText"]))
    story += [
        Paragraph("4. Interpretation and limitations", styles["SectionBlue"]),
        Paragraph(
            "Identity excludes gaps and ambiguous characters; comparable-position counts should be considered "
            "alongside identity percentages. Species labels are inferred from FASTA headers and must be checked. "
            "This report is for research and education, not clinical diagnosis. The fallback alignment is not "
            "equivalent to a full MAFFT multiple sequence alignment. The neighbour-joining tree has no bootstrap "
            "support. Confirm results with appropriate tools before publication or diagnostic use.",
            styles["BodyText"]),
    ]
    doc.build(story)
    return output.getvalue()


# -------------------------------- Navigation ----------------------------------
with st.sidebar:
    st.markdown("## 🧬 BacTrace")
    st.caption("Comparative 16S rRNA analysis")
    st.divider()
    page = st.radio("NAVIGATION", PAGES, key="page")
    st.divider()
    if st.session_state.analysis:
        st.success("Analysis ready")
        st.caption("Last run: " + st.session_state.analysis["timestamp"])
    else:
        st.info("Upload sequences to begin")
    st.caption("Research and educational use only.")

a = st.session_state.analysis

# -------------------------------- Page 1 --------------------------------------
if page == "Home & Input":
    st.markdown("""
    <div class="hero">
      <h1>🧬 BacTrace</h1>
      <p>Comparative 16S rRNA sequence analysis, visualization and reporting.</p>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("### Start a new analysis")
    st.write("Upload a multi-FASTA file or paste nucleotide sequences. Analysis runs on the records supplied in this session.")
    col1, col2 = st.columns(2)
    with col1:
        upload = st.file_uploader("Upload FASTA", type=["fasta", "fa", "fna", "txt"])
    with col2:
        pasted = st.text_area("Or paste multi-FASTA", height=190,
                              placeholder=">seq1 Klebsiella pneumoniae\nACGT...\n>seq2 Klebsiella variicola\nACGA...")
    st.caption("Supported sequence symbols: A, C, G, T and standard IUPAC ambiguity codes. Use unique FASTA IDs.")
    if st.button("Validate and run full analysis", type="primary", use_container_width=True):
        try:
            if upload is not None:
                raw = upload.getvalue().decode("utf-8", errors="replace")
                records = parse_fasta_text(raw)
            else:
                records = parse_fasta_text(pasted)
            with st.spinner("Validating, aligning, calculating metrics and building the tree..."):
                st.session_state.analysis = build_analysis(records)
            st.success(f"Analysis completed for {len(records)} sequences.")
            st.info("Use the navigation menu to explore results and download your report.")
            st.rerun()
        except Exception as exc:
            st.error(str(exc))
    st.markdown("### Workflow")
    cols = st.columns(4)
    for c, num, title, desc in zip(
        cols, ["01", "02", "03", "04"],
        ["Input", "Analyse", "Explore", "Export"],
        ["Upload or paste FASTA", "Run calculations", "Review charts and tree", "Download data and PDF"]
    ):
        with c:
            st.markdown(f"<div class='panel'><h3>{num} · {title}</h3><p>{desc}</p></div>", unsafe_allow_html=True)
    st.warning("16S rRNA similarity alone may not reliably distinguish closely related Klebsiella species.")

# --------------------------- Pages need analysis -------------------------------
elif a is None:
    st.markdown("<div class='hero'><h1>BacTrace</h1><p>Run an analysis to open this section.</p></div>", unsafe_allow_html=True)
    st.info("Go to Home & Input, upload a FASTA file or paste at least two sequences, and run the analysis.")

# -------------------------------- Page 2 --------------------------------------
elif page == "Analysis Overview":
    st.title("Analysis Overview")
    st.caption("Metrics are calculated from the current input.")
    cols = st.columns(4)
    cols[0].metric("Sequences", len(a["records"]))
    cols[1].metric("Alignment columns", a["alignment"].get_alignment_length())
    cols[2].metric("Variable sites", a["cons_summary"]["Variable sites"])
    cols[3].metric("Pairwise comparisons", len(a["pairs"]))
    st.markdown("### Sequence statistics")
    st.dataframe(a["stats"], use_container_width=True, hide_index=True)
    st.markdown("### Methods used")
    st.write("**Alignment:** " + a["alignment_method"])
    st.write("**Tree:** " + a["tree_method"])
    st.write("**Run timestamp:** " + a["timestamp"])
    st.download_button("Download sequence statistics CSV", a["stats"].to_csv(index=False).encode(), "bactrace_sequence_statistics.csv", "text/csv")
    st.warning("Species names are inferred from FASTA headers. Verify labels and accession metadata before interpretation.")

# -------------------------------- Page 3 --------------------------------------
elif page == "Sequence Alignment":
    st.title("Sequence Alignment")
    st.write("**Method:** " + a["alignment_method"])
    st.write(f"**Alignment size:** {len(a['alignment'])} sequences × {a['alignment'].get_alignment_length()} columns")
    st.text_area("Aligned FASTA preview", a["aligned_fasta"], height=360)
    st.download_button("Download aligned FASTA", a["aligned_fasta"].encode(), "bactrace_aligned.fasta", "text/plain")
    if "not equivalent to MAFFT" in a["alignment_method"]:
        st.warning("MAFFT was not available. A reference-based star-alignment fallback was used; for report-quality biological interpretation, install MAFFT and rerun.")
    else:
        st.success("MAFFT alignment completed.")
    st.caption("Inspect sequence completeness and alignment quality before interpreting positional differences.")

# -------------------------------- Page 4 --------------------------------------
elif page == "Conservation Analysis":
    st.title("Conservation Analysis")
    s = a["cons_summary"]
    n = s["Alignment columns"]
    cols = st.columns(4)
    cols[0].metric("Alignment columns", n)
    cols[1].metric("Fully conserved", s["Fully conserved sites"], f"{100*s['Fully conserved sites']/n:.2f}%")
    cols[2].metric("Variable sites", s["Variable sites"], f"{100*s['Variable sites']/n:.2f}%")
    cols[3].metric("Parsimony-informative", s["Parsimony-informative sites"], f"{100*s['Parsimony-informative sites']/n:.2f}%")
    chart = pd.DataFrame({
        "Category": ["Fully conserved", "Variable", "Parsimony-informative"],
        "Count": [s["Fully conserved sites"], s["Variable sites"], s["Parsimony-informative sites"]],
    })
    fig = px.bar(chart, x="Category", y="Count", text="Count", title="Conserved and variable positions")
    fig.update_layout(xaxis_title="", yaxis_title="Alignment columns", showlegend=False)
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("### Position-level results")
    st.dataframe(a["conservation"], use_container_width=True, hide_index=True)
    st.download_button("Download conservation CSV", a["conservation"].to_csv(index=False).encode(), "bactrace_conservation.csv", "text/csv")
    st.info("Definitions: fully conserved means the same unambiguous A/C/G/T in every sequence; variable means at least two unambiguous states; parsimony-informative means at least two states each occur at least twice. Categories are not exhaustive.")

# -------------------------------- Page 5 --------------------------------------
elif page == "Similarity Explorer":
    st.title("Similarity Explorer")
    pairs = a["pairs"]
    valid = pairs.dropna(subset=["Identity (%)"])
    cols = st.columns(3)
    cols[0].metric("Comparable pairs", len(valid))
    cols[1].metric("Mean identity", f"{valid['Identity (%)'].mean():.3f}%" if len(valid) else "N/A")
    cols[2].metric("Mean comparable positions", f"{valid['Comparable positions'].mean():.1f}" if len(valid) else "N/A")
    ids = [r.id for r in a["records"]]
    matrix = pd.DataFrame(np.nan, index=ids, columns=ids)
    for _, row in valid.iterrows():
        matrix.loc[row["Sequence 1"], row["Sequence 2"]] = row["Identity (%)"]
        matrix.loc[row["Sequence 2"], row["Sequence 1"]] = row["Identity (%)"]
    for rid in ids:
        matrix.loc[rid, rid] = 100.0
    fig = px.imshow(matrix, text_auto=".2f", aspect="auto", color_continuous_scale="Teal",
                    title="Pairwise identity (%)", labels={"color": "Identity (%)"})
    fig.update_layout(height=max(450, len(ids)*38))
    st.plotly_chart(fig, use_container_width=True)
    min_positions = st.slider("Minimum comparable positions", 0, max(1, int(valid["Comparable positions"].max()) if len(valid) else 1), 0, 10)
    st.dataframe(pairs[pairs["Comparable positions"] >= min_positions].sort_values("Identity (%)", ascending=False),
                 use_container_width=True, hide_index=True)
    st.download_button("Download pairwise identity CSV", pairs.to_csv(index=False).encode(), "bactrace_pairwise_identity.csv", "text/csv")
    labelled = valid[valid["Comparison type"].isin(["Within-species", "Between-species"])]
    if len(labelled):
        grouped = labelled.groupby("Comparison type")["Identity (%)"].agg(["count", "mean", "min", "max"]).reset_index()
        grouped.columns = ["Comparison type", "Pair count", "Mean identity (%)", "Minimum (%)", "Maximum (%)"]
        st.markdown("### Within-species vs between-species")
        st.dataframe(grouped.round(4), use_container_width=True, hide_index=True)
        fig2 = px.box(labelled, x="Comparison type", y="Identity (%)", points="all", title="Identity distribution")
        st.plotly_chart(fig2, use_container_width=True)
    else:
        st.info("Species comparison requires species labels in FASTA headers, such as 'Klebsiella pneumoniae'.")
    st.caption("Identity excludes gaps and ambiguous characters. Always interpret identity with comparable-position coverage.")

# -------------------------------- Page 6 --------------------------------------
elif page == "Phylogenetic Analysis":
    st.title("Phylogenetic Analysis")
    st.write("**Method:** " + a["tree_method"])
    try:
        st.plotly_chart(tree_figure(a["tree"]), use_container_width=True)
    except Exception as exc:
        st.error(f"Tree display failed: {exc}")
    st.download_button("Download Newick tree", a["newick"].encode(), "bactrace_tree.newick", "text/plain")
    st.markdown("### Interpretation")
    st.write("Interpret the actual topology and branch lengths alongside sequence completeness, alignment quality, and accession metadata.")
    st.warning("The built-in tree is a neighbour-joining tree with no bootstrap support. It is not an IQ-TREE analysis. For a publication-quality tree, use a validated alignment, a suitable evolutionary model, and support estimation.")

# -------------------------------- Page 7 --------------------------------------
elif page == "Report & Export":
    st.title("Report & Export")
    st.write("Every file below is generated from the current analysis.")
    st.download_button("Sequence statistics CSV", a["stats"].to_csv(index=False).encode(), "bactrace_sequence_statistics.csv", "text/csv", use_container_width=True)
    st.download_button("Conservation analysis CSV", a["conservation"].to_csv(index=False).encode(), "bactrace_conservation.csv", "text/csv", use_container_width=True)
    st.download_button("Pairwise identity CSV", a["pairs"].to_csv(index=False).encode(), "bactrace_pairwise_identity.csv", "text/csv", use_container_width=True)
    st.download_button("Aligned FASTA", a["aligned_fasta"].encode(), "bactrace_aligned.fasta", "text/plain", use_container_width=True)
    st.download_button("Newick tree", a["newick"].encode(), "bactrace_tree.newick", "text/plain", use_container_width=True)
    st.markdown("### PDF report")
    try:
        pdf_bytes = create_pdf(a)
        st.download_button("Download PDF analysis report", pdf_bytes, "bactrace_analysis_report.pdf",
                           "application/pdf", type="primary", use_container_width=True)
    except Exception as exc:
        st.error(f"Could not generate PDF: {exc}")
        st.caption("Check that reportlab is installed in requirements.txt.")
    st.markdown("### Scientific notes")
    st.write("- Species labels are inferred from FASTA headers and must be verified.")
    st.write("- Alignment method and tree method are stated explicitly.")
    st.write("- This dashboard is intended for research and education, not clinical diagnosis.")

st.markdown("---")
st.markdown("<p class='small-note'>BacTrace · Comparative 16S rRNA sequence analysis · Research and educational use only</p>", unsafe_allow_html=True)
