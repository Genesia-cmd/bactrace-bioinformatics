def make_pdf(results: dict) -> bytes:
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer, pagesize=A4, rightMargin=17 * mm, leftMargin=17 * mm,
        topMargin=17 * mm, bottomMargin=17 * mm, title="BacTrace Analysis Report"
    )
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportTitle", parent=styles["Title"], alignment=TA_CENTER,
                              textColor=colors.HexColor("#102A43"), fontSize=20, leading=24, spaceAfter=8))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading2"],
                              textColor=colors.HexColor("#087E8B"), spaceBefore=10, spaceAfter=6))
    styles.add(ParagraphStyle(name="SmallBody", parent=styles["BodyText"], fontSize=8.5, leading=11))
    story = []
    stats = results["stats"]
    cons = results["conservation"]
    pairs = results["pairs"]
    alignment = results["alignment"]
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    story += [
        Paragraph("BacTrace", styles["ReportTitle"]),
        Paragraph("Comparative 16S rRNA Sequence Analysis Report", styles["Heading2"]),
        Paragraph(f"Generated: {generated} | Application version: {APP_VERSION}", styles["SmallBody"]),
        Spacer(1, 8),
        Paragraph("1. Executive summary", styles["Section"]),
    ]
    mean_identity = pairs["Identity (%)"].mean() if not pairs.empty else np.nan
    high = pairs.loc[pairs["Identity (%)"].idxmax()] if not pairs.empty and pairs["Identity (%)"].notna().any() else None
    low = pairs.loc[pairs["Identity (%)"].idxmin()] if not pairs.empty and pairs["Identity (%)"].notna().any() else None
    summary = (
        f"The active dataset contains {len(stats)} sequences and an alignment of "
        f"{alignment.get_alignment_length()} columns using {results['alignment_method']}. "
        f"Mean pairwise identity across calculable pairs is "
        f"{mean_identity:.2f}%." if pd.notna(mean_identity) else
        f"The active dataset contains {len(stats)} sequences. No pair had comparable unambiguous positions."
    )
    story.append(Paragraph(summary, styles["SmallBody"]))
    story += [Paragraph("2. Input dataset and sequence statistics", styles["Section"])]
    story.append(dataframe_table(stats, styles))
    story += [
        Paragraph("3. Methods", styles["Section"]),
        Paragraph(
            "Sequences were parsed as FASTA or raw nucleotide input. U was normalized to T. "
            "GC content uses unambiguous A/C/G/T positions as its denominator. Pairwise identity "
            "counts only aligned positions where both sequences contain A, C, G, or T; gaps and "
            "ambiguity codes are excluded. The multiple alignment method used for this run is "
            + results["alignment_method"] + ". The phylogeny, if available, uses neighbor joining "
            "with distance 1 − identity calculated over comparable aligned bases. No bootstrap support was calculated.",
            styles["SmallBody"]
        ),
        Paragraph("4. Alignment and conservation", styles["Section"]),
    ]
    ncol = alignment.get_alignment_length()
    fully = int(cons["Fully conserved"].sum())
    variable = int(cons["Variable site"].sum())
    informative = int(cons["Parsimony-informative"].sum())
    gaps = int(cons["Gap-containing"].sum())
    story.append(Paragraph(
        f"Alignment columns: {ncol}; fully conserved sites: {fully} ({fully/ncol*100:.2f}%); "
        f"variable sites: {variable} ({variable/ncol*100:.2f}%); parsimony-informative sites: "
        f"{informative} ({informative/ncol*100:.2f}%); gap-containing columns: {gaps}. "
        "These categories are not mutually exclusive and do not necessarily sum to the alignment length.",
        styles["SmallBody"]
    ))
    story += [Paragraph("5. Pairwise identity", styles["Section"])]
    if high is not None:
        story.append(Paragraph(
            f"Highest identity: {high['Sequence 1']} vs {high['Sequence 2']} "
            f"({high['Identity (%)']:.2f}%, {int(high['Comparable positions'])} comparable positions).",
            styles["SmallBody"]))
    if low is not None:
        story.append(Paragraph(
            f"Lowest identity: {low['Sequence 1']} vs {low['Sequence 2']} "
            f"({low['Identity (%)']:.2f}%, {int(low['Comparable positions'])} comparable positions).",
            styles["SmallBody"]))
    story.append(dataframe_table(pairs, styles, max_rows=24))
    if not results["groups"].empty:
        story += [Paragraph("6. Species-group comparisons (user-supplied labels)", styles["Section"]),
                  dataframe_table(results["groups"], styles)]
    story += [Paragraph("7. Phylogenetic analysis", styles["Section"])]
    if results.get("newick"):
        story.append(Paragraph(
            "A neighbor-joining tree was generated from the active alignment. The tree is exploratory, "
            "has no bootstrap support values, and is not sufficient by itself for definitive species identification.",
            styles["SmallBody"]))
        story.append(Paragraph("Newick: " + results["newick"][:1200], styles["SmallBody"]))
    else:
        story.append(Paragraph("Tree generation was not successful or was not available for this dataset: " +
                               results.get("tree_error", "no tree generated"), styles["SmallBody"]))
    story += [
        Paragraph("8. Interpretation and quality checks", styles["Section"]),
        Paragraph(results["interpretation"].replace("&", "&amp;"), styles["SmallBody"]),
        Paragraph("9. Limitations and recommendations", styles["Section"]),
        Paragraph(
            "Results depend on sequence quality, homology, length, and alignment accuracy. Short or non-overlapping "
            "regions can make comparisons misleading. The progressive fallback is a lightweight educational method "
            "and may be less accurate than MAFFT. 16S rRNA may not resolve closely related taxa, including some "
            "Klebsiella species. Validate results with curated references, additional markers, or genome-based methods "
            "before publication or diagnostic use. This tool is educational and not a clinical diagnostic system.",
            styles["SmallBody"]),
        Paragraph("10. Conclusion", styles["Section"]),
        Paragraph(
            "BacTrace summarizes sequence composition, alignment-derived conservation, strict pairwise identity, "
            "and exploratory phylogenetic relationships for the supplied dataset. Interpret all outputs in the context "
            "of sequence coverage and the stated methodological limitations.", styles["SmallBody"]),
    ]
    doc.build(story, onFirstPage=pdf_page, onLaterPages=pdf_page)
    return buffer.getvalue()


def dataframe_table(df: pd.DataFrame, styles, max_rows: int = 35):
    if df is None or df.empty:
        return Paragraph("No results available for this section.", styles["SmallBody"])
    view = df.head(max_rows).copy()
    headers = [str(c) for c in view.columns]
    data = [[Paragraph(str(h), styles["SmallBody"]) for h in headers]]
    for _, row in view.iterrows():
        data.append([Paragraph(str(v)[:110], styles["SmallBody"]) for v in row.tolist()])
    available_width = A4[0] - 34 * mm
    widths = [available_width / len(headers)] * len(headers)
    table = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#102A43")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D9E2EC")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7FA")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    elements = [table]
    if len(df) > max_rows:
        elements.append(Paragraph(f"Showing first {max_rows} of {len(df)} rows in the PDF. Download CSV for all rows.", styles["SmallBody"]))
    return KeepTogether(elements)


def pdf_page(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D9E2EC"))
    canvas.line(17 * mm, 13 * mm, A4[0] - 17 * mm, 13 * mm)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#52606D"))
    canvas.drawString(17 * mm, 8 * mm, "BacTrace | Educational bioinformatics analysis")
    canvas.drawRightString(A4[0] - 17 * mm, 8 * mm, f"Page {doc.page}")
    canvas.restoreState()


def build_interpretation(stats, cons, pairs, warnings, alignment_method):
    lines = []
    lines.append(f"Analysed {len(stats)} sequences with alignment method: {alignment_method}.")
    if not pairs.empty and pairs["Identity (%)"].notna().any():
        hi = pairs.loc[pairs["Identity (%)"].idxmax()]
        lo = pairs.loc[pairs["Identity (%)"].idxmin()]
        lines.append(f"Highest observed pairwise identity: {hi['Sequence 1']} vs {hi['Sequence 2']} = "
                     f"{hi['Identity (%)']:.2f}% across {int(hi['Comparable positions'])} comparable positions.")
        lines.append(f"Lowest observed pairwise identity: {lo['Sequence 1']} vs {lo['Sequence 2']} = "
                     f"{lo['Identity (%)']:.2f}% across {int(lo['Comparable positions'])} comparable positions.")
        lines.append(f"Mean pairwise identity: {pairs['Identity (%)'].mean():.2f}%. Identity is not proof of species identity.")
    else:
        lines.append("No pair had enough comparable unambiguous positions for an identity percentage.")
    n = len(cons)
    if n:
        lines.append(f"Fully conserved sites: {int(cons['Fully conserved'].sum())}/{n}; "
                     f"variable sites: {int(cons['Variable site'].sum())}/{n}; "
                     f"parsimony-informative sites: {int(cons['Parsimony-informative'].sum())}/{n}.")
        gap_cols = int(cons["Gap-containing"].sum())
        amb_cols = int(cons["Ambiguity-containing"].sum())
        lines.append(f"Columns containing gaps: {gap_cols}; columns containing ambiguity codes: {amb_cols}.")
    min_len, max_len = int(stats["Length (nt)"].min()), int(stats["Length (nt)"].max())
    if min_len != max_len:
        lines.append(f"Input lengths range from {min_len} to {max_len} nt; uneven lengths can indicate partial coverage or non-homologous regions.")
    if not warnings:
        lines.append("No validation warnings were raised, but sequence homology and biological provenance still require user confirmation.")
    else:
        lines.append("Quality warnings: " + " ".join(warnings))
    lines.append("Recommendation: verify accession metadata and homologous-region coverage; use curated reference sequences and additional markers or genome-based methods when species-level resolution is needed.")
    return " ".join(lines)


def perform_analysis(records, species_text):
    stats = sequence_statistics(records)
    alignment, alignment_method = run_mafft(records)
    cons = conservation_analysis(alignment)
    pairs = pairwise_identity(alignment)
    labels = parse_species_labels(species_text, [r.id for r in records])
    groups, group_note = group_comparisons(pairs, labels)
    tree, newick, tree_error = None, "", ""
    try:
        tree, _ = build_tree(alignment)
        newick = tree_to_newick(tree)
    except Exception as exc:
        tree_error = str(exc)
    interpretation = build_interpretation(stats, cons, pairs, [], alignment_method)
    return {
        "records": records, "stats": stats, "alignment": alignment,
        "alignment_method": alignment_method, "conservation": cons, "pairs": pairs,
        "labels": labels, "groups": groups, "group_note": group_note,
        "tree": tree, "newick": newick, "tree_error": tree_error,
        "interpretation": interpretation,
        "warnings": [],
    }


def render_sidebar():
    st.sidebar.markdown("## 🧬 BacTrace")
    st.sidebar.caption("Comparative 16S rRNA sequence analysis")
    current = st.session_state.get("page", PAGES[0])
    selected = st.sidebar.radio("Navigate pages", PAGES, index=PAGES.index(current), key="page_radio")
    if selected != current:
        st.session_state.page = selected
    st.sidebar.markdown("---")
    st.sidebar.caption("**Workflow**")
    if st.session_state.get("results"):
        st.sidebar.success("Analysis available in this session")
    else:
        st.sidebar.info("Submit at least two sequences to begin")
    if st.sidebar.button("＋ Start new analysis", use_container_width=True):
        for key in ["results", "input_records", "validation_warnings", "species_text"]:
            st.session_state.pop(key, None)
        st.session_state.page = PAGES[0]
        st.rerun()
    st.sidebar.markdown("---")
    st.sidebar.caption("Educational tool • Not a diagnostic system")


def render_home():
    st.markdown("# BacTrace")
    st.markdown("### Comparative 16S rRNA Sequence Analysis and Bacterial Identification Platform")
    st.markdown(
        "Transform nucleotide sequences into alignment, conservation, similarity and exploratory phylogenetic results. "
        "All reported metrics are calculated from the active input dataset."
    )
    st.info("For meaningful comparison, submit sequences from the same gene or homologous region. 16S rRNA similarity alone does not confirm species identity.")
    uploaded = st.file_uploader(
        "Option A — Upload FASTA file(s)",
        type=["fasta", "fa", "fas", "fna", "txt"],
        accept_multiple_files=True,
        help="Each file may contain one or more FASTA records. If multiple files are selected, records are combined."
    )
    pasted = st.text_area(
        "Option B — Paste multi-FASTA or a single raw nucleotide sequence",
        height=220,
        placeholder=">Sequence_A\\nAGAGTTTGATCCTGGCTCAG...\\n>Sequence_B\\nAGAGTTTGATCCTGGCTCAG..."
    )
    col1, col2 = st.columns([1, 1])
    with col1:
        if st.button("Load clearly labelled example dataset", use_container_width=True):
            st.session_state.example_loaded = True
            st.session_state.pasted_input = EXAMPLE_FASTA
            st.rerun()
    with col2:
        st.caption("The example sequences are demonstration inputs, not curated reference sequences.")
    if "pasted_input" in st.session_state and not pasted:
        pasted = st.session_state.pasted_input
    species_text = st.text_area(
        "Optional species labels (one per line: sequence_id<TAB>species name)",
        height=90,
        placeholder="Example_16S_A\\tKlebsiella pneumoniae\\nExample_16S_B\\tKlebsiella variicola"
    )
    st.caption("Labels are used only for grouped comparisons; they are not inferred or verified.")
    if st.button("Analyse Sequences", type="primary", use_container_width=True):
        try:
            records = []
            if uploaded:
                for file in uploaded:
                    raw = file.getvalue().decode("utf-8-sig", errors="replace")
                    records.extend(parse_fasta_text(raw))
            elif pasted.strip():
                records = parse_fasta_text(pasted)
            else:
                st.error("Upload a FASTA file or paste nucleotide sequences first.")
                return
            records, warnings = validate_records(records)
            # A second check catches duplicate IDs across multiple uploaded files.
            ids = [r.id for r in records]
            if len(ids) != len(set(ids)):
                raise ValueError("Duplicate sequence identifiers occur across uploaded files. Rename them before analysis.")
            with st.spinner("Calculating sequence statistics, multiple alignment, conservation, identity and tree..."):
                progress = st.progress(0, text="Validating input")
                progress.progress(15, text="Calculating sequence statistics")
                results = perform_analysis(records, species_text)
                progress.progress(100, text="Analysis complete")
            results["warnings"] = warnings
            results["interpretation"] = build_interpretation(
                results["stats"], results["conservation"], results["pairs"], warnings, results["alignment_method"]
            )
            st.session_state.results = results
            st.session_state.validation_warnings = warnings
            st.session_state.species_text = species_text
            st.session_state.page = "Analysis Overview"
            st.success("Analysis completed from the submitted sequences.")
            st.rerun()
        except Exception as exc:
            st.error(str(exc))
    st.markdown("---")
    st.markdown("**Input rules**")
    st.markdown(
        "- At least two nucleotide sequences are required for comparative analysis.\n"
        "- Standard IUPAC ambiguity codes are accepted; ambiguity is excluded from strict identity calculations.\n"
        "- Limit: 40 sequences and 250,000 total bases per run.\n"
        "- RNA U characters are normalized to T; original input is retained in the active records as normalized DNA."
    )


def require_results():
    results = st.session_state.get("results")
    if not results:
        st.warning("No analysis is available yet. Go to **Home & Input** and analyse at least two sequences.")
        if st.button("Go to Home & Input"):
            st.session_state.page = PAGES[0]
            st.rerun()
        return None
    return results


def render_overview(r):
    st.title("Analysis Overview")
    stats, cons, pairs = r["stats"], r["conservation"], r["pairs"]
    nseq = len(stats)
    ncol = r["alignment"].get_alignment_length()
    mean_gc = stats["GC content (%)"].mean()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Sequences analysed", nseq)
    c2.metric("Alignment length", f"{ncol:,} nt")
    c3.metric("Mean GC content", f"{mean_gc:.2f}%" if pd.notna(mean_gc) else "N/A")
    c4.metric("Mean pairwise identity", f"{pairs['Identity (%)'].mean():.2f}%" if not pairs.empty and pairs["Identity (%)"].notna().any() else "N/A")
    st.caption(f"Alignment method actually used: **{r['alignment_method']}**")
    if r["warnings"]:
        with st.expander("Input quality warnings", expanded=True):
            for w in r["warnings"]:
                st.warning(w)
    a, b, c = st.columns(3)
    a.metric("Fully conserved sites", int(cons["Fully conserved"].sum()))
    b.metric("Variable sites", int(cons["Variable site"].sum()))
    c.metric("Parsimony-informative sites", int(cons["Parsimony-informative"].sum()))
    st.subheader("Sequence statistics")
    st.dataframe(stats, use_container_width=True, hide_index=True)
    st.download_button("Download sequence statistics CSV", stats.to_csv(index=False).encode("utf-8"), "BacTrace_Sequence_Statistics.csv", "text/csv")
    chart = px.bar(stats, x="Sequence ID", y="GC content (%)", title="GC content by sequence")
    chart.update_layout(margin=dict(l=10, r=10, t=50, b=10))
    st.plotly_chart(chart, use_container_width=True)
    st.subheader("Calculated interpretation")
    st.write(r["interpretation"])


def render_alignment(r):
    st.title("Sequence Alignment")
    aln = r["alignment"]
    st.caption(f"Method: {r['alignment_method']} • {len(aln)} sequences • {aln.get_alignment_length()} alignment columns")
    st.warning("The built-in fallback is a lightweight progressive profile alignment and may be less accurate than MAFFT. Inspect alignments before drawing biological conclusions.")
    width = st.slider("Alignment preview width (columns)", min_value=20, max_value=200, value=80, step=10)
    start = st.number_input("Start position (1-based)", min_value=1, max_value=max(1, aln.get_alignment_length()), value=1, step=1)
    end = min(aln.get_alignment_length(), int(start) + width - 1)
    st.caption(f"Showing alignment positions {int(start)}–{end}")
    lines = []
    for rec in aln:
        seq = str(rec.seq)[int(start)-1:end]
        # Wrap visually with a monospaced code block; the user can inspect gaps directly.
        lines.append(f"{rec.id[:28]:28}  {seq}")
    st.code("\n".join(lines), language="text")
    st.download_button("Download aligned FASTA", aligned_fasta_string(aln).encode("utf-8"),
                       "BacTrace_Aligned_Sequences.fasta", "text/plain")
    st.markdown("**How to read this view:** `-` represents an alignment gap. The displayed alignment is computed from the submitted sequences.")


def render_conservation(r):
    st.title("Conservation Analysis")
    df = r["conservation"]
    n = len(df)
    a, b, c, d = st.columns(4)
    a.metric("Alignment columns", n)
    b.metric("Fully conserved", f"{int(df['Fully conserved'].sum())} ({df['Fully conserved'].mean()*100:.1f}%)" if n else "0")
    c.metric("Variable", f"{int(df['Variable site'].sum())} ({df['Variable site'].mean()*100:.1f}%)" if n else "0")
    d.metric("Parsimony-informative", f"{int(df['Parsimony-informative'].sum())} ({df['Parsimony-informative'].mean()*100:.1f}%)" if n else "0")
    st.caption("Site categories can overlap or leave some positions unclassified; do not sum their counts as if they were mutually exclusive.")
    chart_df = pd.DataFrame({
        "Site category": ["Fully conserved", "Variable", "Parsimony-informative", "Gap-containing", "Ambiguity-containing"],
        "Count": [int(df[c].sum()) for c in ["Fully conserved", "Variable site", "Parsimony-informative", "Gap-containing", "Ambiguity-containing"]]
    })
    st.plotly_chart(px.bar(chart_df, x="Site category", y="Count", title="Alignment site categories"), use_container_width=True)
    filt = st.multiselect("Filter positions by property", ["Fully conserved", "Variable site", "Parsimony-informative", "Gap-containing", "Ambiguity-containing"])
    view = df.copy()
    if filt:
        view = view[view[filt].any(axis=1)]
    search = st.text_input("Search nucleotide pattern or position", placeholder="e.g. 120 or ACG")
    if search.strip():
        if search.strip().isdigit():
            view = view[view["Alignment position"].astype(str).str.contains(search.strip())]
        else:
            view = view[view["Nucleotide pattern"].str.contains(search.strip().upper(), regex=False)]
    st.dataframe(view, use_container_width=True, hide_index=True)
    st.download_button("Download conservation CSV", df.to_csv(index=False).encode("utf-8"),
                       "BacTrace_Conservation_Results.csv", "text/csv")


def render_similarity(r):
    st.title("Similarity Explorer")
    pairs = r["pairs"]
    if pairs.empty:
        st.info("No unique sequence pairs are available.")
    else:
        valid = pairs.dropna(subset=["Identity (%)"])
        if not valid.empty:
            high = valid.loc[valid["Identity (%)"].idxmax()]
            low = valid.loc[valid["Identity (%)"].idxmin()]
            x, y, z = st.columns(3)
            x.metric("Highest identity", f"{high['Identity (%)']:.2f}%")
            y.metric("Lowest identity", f"{low['Identity (%)']:.2f}%")
            z.metric("Mean identity", f"{valid['Identity (%)'].mean():.2f}%")
            st.caption(f"Highest: {high['Sequence 1']} ↔ {high['Sequence 2']} ({int(high['Comparable positions'])} comparable positions)")
            st.caption(f"Lowest: {low['Sequence 1']} ↔ {low['Sequence 2']} ({int(low['Comparable positions'])} comparable positions)")
            chart = px.bar(valid, x="Sequence 1", y="Identity (%)", color="Sequence 2",
                           hover_data=["Comparable positions", "Comparable coverage (%)"],
                           title="Pairwise identity (comparable positions shown in hover)")
            st.plotly_chart(chart, use_container_width=True)
        q = st.text_input("Filter by sequence identifier")
        view = pairs
        if q.strip():
            view = pairs[pairs["Sequence 1"].str.contains(q, case=False, regex=False) |
                         pairs["Sequence 2"].str.contains(q, case=False, regex=False)]
        st.dataframe(view, use_container_width=True, hide_index=True)
        st.download_button("Download pairwise identity CSV", pairs.to_csv(index=False).encode("utf-8"),
                           "BacTrace_Pairwise_Identity.csv", "text/csv")
    st.markdown("---")
    st.subheader("Within-species and between-species comparisons")
    st.caption("Only user-supplied labels are used. Labels are not inferred from sequence names and are not independently verified.")
    if r["labels"]:
        st.write(f"Labels supplied for {len(r['labels'])} sequence(s):")
        st.dataframe(pd.DataFrame([{"Sequence ID": k, "Species label": v} for k, v in r["labels"].items()]), hide_index=True)
    if not r["groups"].empty:
        st.dataframe(r["groups"], use_container_width=True, hide_index=True)
    else:
        st.info(r["group_note"] or "No group summary is available. Add at least two labelled sequences in one species group for within-species identity, and multiple groups for between-species comparisons.")
    if r["group_note"]:
        st.caption(r["group_note"])
    st.warning("High 16S rRNA identity alone does not establish species identity. This is especially important for closely related Klebsiella taxa.")


def render_phylogeny(r):
    st.title("Phylogenetic Analysis")
    if not r["tree"]:
        st.error("Tree generation was not successful for this analysis.")
        st.code(r["tree_error"] or "No further error detail was provided.")
        return
    st.caption("Method: neighbor joining using distance = 1 − strict pairwise identity across comparable aligned A/C/G/T positions.")
    st.warning("Exploratory tree only. No bootstrap support was calculated. Branch lengths and topology depend on the alignment and chosen distance definition.")
    segments, tips = tree_to_plot_df(r["tree"])
    import plotly.graph_objects as go
    fig = go.Figure()
    for x0, y0, x1, y1 in segments:
        fig.add_trace(go.Scatter(x=[x0, x1], y=[y0, y0], mode="lines",
                                 line=dict(color="#087E8B", width=2), hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter(x=[x1, x1], y=[y0, y1], mode="lines",
                                 line=dict(color="#087E8B", width=2), hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=tips["x"], y=tips["y"], mode="markers+text",
                             text=tips["Sequence"], textposition="middle right",
                             marker=dict(size=7, color="#102A43"), hovertemplate="%{text}<extra></extra>",
                             showlegend=False))
    fig.update_layout(title="Neighbor-joining tree", xaxis_title="Distance from root (relative)",
                      yaxis=dict(visible=False, autorange="reversed"),
                      showlegend=False, height=max(350, 55 * len(tips)),
                      margin=dict(l=10, r=150, t=50, b=40))
    st.plotly_chart(fig, use_container_width=True)
    st.download_button("Download Newick tree", r["newick"].encode("utf-8"),
                       "BacTrace_Phylogenetic_Tree.nwk", "text/plain")
    st.code(r["newick"], language="text")
    st.markdown("**Interpretation note:** A preliminary 16S rRNA tree is not sufficient to definitively identify closely related bacterial species. Consider curated references and additional markers or genome-based methods.")


def render_exports(r):
    st.title("Report & Export")
    st.markdown("Generate a report from the **currently active analysis**. Starting a new analysis replaces the active results.")
    c1, c2 = st.columns(2)
    with c1:
        st.download_button("Sequence statistics CSV", r["stats"].to_csv(index=False).encode("utf-8"),
                           "BacTrace_Sequence_Statistics.csv", "text/csv", use_container_width=True)
        st.download_button("Conservation results CSV", r["conservation"].to_csv(index=False).encode("utf-8"),
                           "BacTrace_Conservation_Results.csv", "text/csv", use_container_width=True)
        st.download_button("Pairwise identity CSV", r["pairs"].to_csv(index=False).encode("utf-8"),
                           "BacTrace_Pairwise_Identity.csv", "text/csv", use_container_width=True)
    with c2:
        st.download_button("Aligned sequences FASTA", aligned_fasta_string(r["alignment"]).encode("utf-8"),
                           "BacTrace_Aligned_Sequences.fasta", "text/plain", use_container_width=True)
        if r["newick"]:
            st.download_button("Phylogenetic tree Newick", r["newick"].encode("utf-8"),
                               "BacTrace_Phylogenetic_Tree.nwk", "text/plain", use_container_width=True)
        else:
            st.caption("Newick download unavailable because tree generation did not succeed.")
    st.markdown("---")
    st.subheader("Complete PDF report")
    if st.button("Generate PDF report", type="primary"):
        try:
            with st.spinner("Preparing report from the active results..."):
                pdf = make_pdf(r)
            st.session_state.report_pdf = pdf
            st.success("PDF report generated.")
        except Exception as exc:
            st.error(f"PDF report generation failed: {exc}")
    if st.session_state.get("report_pdf"):
        st.download_button("Download BacTrace Analysis Report PDF", st.session_state.report_pdf,
                           "BacTrace_Analysis_Report.pdf", "application/pdf", use_container_width=True)
    st.markdown("---")
    st.caption("The PDF contains calculated statistics, method descriptions, pairwise comparisons, a tree summary when available, interpretation and limitations.")


def main():
    st.set_page_config(page_title="BacTrace | 16S rRNA Analysis", page_icon="🧬", layout="wide",
                       initial_sidebar_state="expanded")
    st.markdown("""
    <style>
    .stApp { background: #F7F9FC; color: #102A43; }
    [data-testid="stSidebar"] { background: #102A43; }
    [data-testid="stSidebar"] * { color: #F0F4F8 !important; }
    h1, h2, h3 { color: #102A43; }
    div.stButton > button[kind="primary"] { background: #087E8B; border-color: #087E8B; color: white; }
    div[data-testid="stMetric"] { background: white; padding: 14px; border-radius: 10px; border: 1px solid #D9E2EC; }
    div[data-testid="stDataFrame"] { background: white; border-radius: 8px; }
    </style>
    """, unsafe_allow_html=True)
    if "page" not in st.session_state:
        st.session_state.page = PAGES[0]
    render_sidebar()
    page = st.session_state.page
    if page == "Home & Input":
        render_home()
    else:
        r = require_results()
        if r is None:
            return
        if page == "Analysis Overview":
            render_overview(r)
        elif page == "Sequence Alignment":
            render_alignment(r)
        elif page == "Conservation Analysis":
            render_conservation(r)
        elif page == "Similarity Explorer":
            render_similarity(r)
        elif page == "Phylogenetic Analysis":
            render_phylogeny(r)
        elif page == "Report & Export":
            render_exports(r)
    st.markdown("---")
    current_index = PAGES.index(st.session_state.page)
    prev_col, label_col, next_col = st.columns([1, 2, 1])
    with prev_col:
        if st.button("← Previous", disabled=current_index == 0, use_container_width=True):
            st.session_state.page = PAGES[current_index - 1]
            st.rerun()
    with label_col:
        st.markdown(f"<div style='text-align:center;padding-top:8px;color:#52606D'>Page {current_index + 1} of {len(PAGES)} · {PAGES[current_index]}</div>", unsafe_allow_html=True)
    with next_col:
        if st.button("Next →", disabled=current_index == len(PAGES) - 1, use_container_width=True):
            st.session_state.page = PAGES[current_index + 1]
            st.rerun()


if __name__ == "__main__":
    main()
