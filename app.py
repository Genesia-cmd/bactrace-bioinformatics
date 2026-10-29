        st.rerun()

# ---------------------------- Page 2 ------------------------------------------
elif selected_page == "Analysis Overview":
    st.title("Analysis Overview")
    st.caption("Summary metrics calculated from the current input.")
    c = st.columns(4)
    c[0].metric("Sequences", len(analysis["records"]))
    c[1].metric("Alignment columns", analysis["alignment"].get_alignment_length())
    c[2].metric("Variable sites", analysis["cons_summary"]["Variable sites"])
    c[3].metric("Pairwise comparisons", len(analysis["pairs"]))
    st.markdown("### Sequence summary")
    st.dataframe(analysis["summary"], use_container_width=True, hide_index=True)
    st.markdown("### Analysis methods")
    st.write(f"**Alignment:** {analysis['alignment_method']}")
    st.write(f"**Phylogenetic tree:** {analysis['tree_method']}")
    st.write(f"**Analysis timestamp:** {analysis['timestamp']}")
    st.download_button("Download sequence summary (CSV)", analysis["summary"].to_csv(index=False).encode("utf-8"), "bactrace_sequence_summary.csv", "text/csv")
    st.warning("Check sequence descriptions and accession metadata. Species labels are inferred from FASTA headers only and may be unassigned or incorrect.")

# ---------------------------- Page 3 ------------------------------------------
elif selected_page == "Sequence Alignment":
    st.title("Sequence Alignment")
    st.caption(f"Method: {analysis['alignment_method']}")
    st.write(f"**Alignment dimensions:** {len(analysis['alignment'])} sequences × {analysis['alignment'].get_alignment_length()} columns")
    st.markdown("### Aligned FASTA")
    st.text_area("Aligned sequences (read-only preview)", analysis["aligned_fasta"], height=350)
    st.download_button("Download aligned FASTA", analysis["aligned_fasta"].encode("utf-8"), "bactrace_aligned.fasta", "text/plain")
    st.markdown("### Alignment note")
    if "fallback" in analysis["alignment_method"].lower():
        st.warning("MAFFT was not found or did not run. A simplified pairwise fallback was used; this is not equivalent to a full multiple sequence alignment. For report-quality results, install MAFFT and rerun.")
    else:
        st.success("MAFFT alignment completed.")
    st.caption("Inspect alignment quality and sequence completeness before interpreting conserved or variable positions.")

# ---------------------------- Page 4 ------------------------------------------
elif selected_page == "Conservation Analysis":
    st.title("Conservation Analysis")
    cons = analysis["cons_summary"]
    ncols = cons["Alignment columns"]
    cols = st.columns(4)
    cols[0].metric("Alignment columns", ncols)
    cols[1].metric("Fully conserved", cons["Fully conserved sites"], f"{100*cons['Fully conserved sites']/ncols:.2f}%")
    cols[2].metric("Variable sites", cons["Variable sites"], f"{100*cons['Variable sites']/ncols:.2f}%")
    cols[3].metric("Parsimony-informative", cons["Parsimony-informative sites"], f"{100*cons['Parsimony-informative sites']/ncols:.2f}%")
    chart_df = pd.DataFrame({
        "Category": ["Fully conserved", "Variable", "Parsimony-informative"],
        "Count": [cons["Fully conserved sites"], cons["Variable sites"], cons["Parsimony-informative sites"]],
    })
    fig = px.bar(chart_df, x="Category", y="Count", text="Count", title="Conservation and variation summary")
    fig.update_layout(yaxis_title="Number of alignment columns", xaxis_title="", showlegend=False)
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("### Position-level results")
    st.dataframe(analysis["cons_df"], use_container_width=True, hide_index=True)
    st.download_button("Download conservation table (CSV)", analysis["cons_df"].to_csv(index=False).encode("utf-8"), "bactrace_conservation.csv", "text/csv")
    st.info("Definitions: fully conserved = same unambiguous A/C/G/T in every sequence; variable = at least two unambiguous nucleotide states; parsimony-informative = at least two states each appearing at least twice. These categories do not necessarily sum to the alignment length.")

# ---------------------------- Page 5 ------------------------------------------
elif selected_page == "Similarity Explorer":
    st.title("Similarity Explorer")
    pairs = analysis["pairs"].copy()
    valid = pairs.dropna(subset=["Identity (%)"])
    c = st.columns(3)
    c[0].metric("Unique pairs", len(valid))
    c[1].metric("Mean identity", f"{valid['Identity (%)'].mean():.3f}%")
    c[2].metric("Mean comparable positions", f"{valid['Comparable positions'].mean():.1f}")
    ids = [r.id for r in analysis["records"]]
    matrix = pd.DataFrame(np.nan, index=ids, columns=ids)
    for _, row in valid.iterrows():
        matrix.loc[row["Sequence 1"], row["Sequence 2"]] = row["Identity (%)"]
        matrix.loc[row["Sequence 2"], row["Sequence 1"]] = row["Identity (%)"]
    np.fill_diagonal(matrix.values, 100.0)
    fig = px.imshow(matrix, text_auto=".2f", aspect="auto", color_continuous_scale="Teal",
                    title="Pairwise sequence identity (%)", labels=dict(color="Identity (%)"))
    fig.update_layout(height=max(480, len(ids)*38))
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("### Pairwise comparison table")
    min_cov = st.slider("Filter by minimum comparable positions", min_value=0, max_value=max(1, int(valid["Comparable positions"].max())), value=0, step=10)
    filtered = pairs[pairs["Comparable positions"] >= min_cov]
    st.dataframe(filtered.sort_values("Identity (%)", ascending=False), use_container_width=True, hide_index=True)
    st.download_button("Download pairwise table (CSV)", pairs.to_csv(index=False).encode("utf-8"), "bactrace_pairwise_identity.csv", "text/csv")
    labelled = valid[valid["Comparison type"].isin(["Within-species", "Between-species"])]
    if not labelled.empty:
        st.markdown("### Within-species vs between-species")
        grouped = labelled.groupby("Comparison type")["Identity (%)"].agg(["count", "mean", "min", "max"]).reset_index()
        grouped.columns = ["Comparison type", "Pair count", "Mean identity (%)", "Minimum (%)", "Maximum (%)"]
        st.dataframe(grouped.round(4), use_container_width=True, hide_index=True)
        fig2 = px.box(labelled, x="Comparison type", y="Identity (%)", points="all", title="Identity distribution by comparison type")
        st.plotly_chart(fig2, use_container_width=True)
    else:
        st.info("Species-group comparison is unavailable because species names were not detected in the FASTA headers. Include reliable species labels in record descriptions to enable it.")
    st.caption("Identity excludes gaps and ambiguous characters. Always interpret the identity percentage alongside comparable-position coverage.")

# ---------------------------- Page 6 ------------------------------------------
elif selected_page == "Phylogenetic Analysis":
    st.title("Phylogenetic Analysis")
    st.write(f"**Method:** {analysis['tree_method']}")
    try:
        st.plotly_chart(draw_tree(analysis["tree"]), use_container_width=True)
    except Exception as exc:
        st.error(f"Tree visualization could not be rendered: {exc}")
    st.download_button("Download Newick tree", analysis["newick"].encode("utf-8"), "bactrace_tree.newick", "text/plain")
    st.markdown("### Interpretation guidance")
    st.write("The tree summarizes relationships inferred from the current alignment and distance method. Examine the actual topology, branch lengths, sequence coverage, and accession metadata before drawing biological conclusions.")
    st.warning("This local neighbour-joining tree has no bootstrap support values. It is not an IQ-TREE result. For a publication-quality phylogeny, use a validated alignment and a suitable evolutionary model with support estimation.")

# ---------------------------- Page 7 ------------------------------------------
elif selected_page == "Report & Export":
    st.title("Report & Export")
    st.write("Download the outputs generated from the current analysis.")
    st.markdown("### Available files")
    st.download_button("Sequence summary (CSV)", analysis["summary"].to_csv(index=False).encode("utf-8"), "bactrace_sequence_summary.csv", "text/csv", use_container_width=True)
    st.download_button("Conservation table (CSV)", analysis["cons_df"].to_csv(index=False).encode("utf-8"), "bactrace_conservation.csv", "text/csv", use_container_width=True)
    st.download_button("Pairwise identity (CSV)", analysis["pairs"].to_csv(index=False).encode("utf-8"), "bactrace_pairwise_identity.csv", "text/csv", use_container_width=True)
    st.download_button("Aligned FASTA", analysis["aligned_fasta"].encode("utf-8"), "bactrace_aligned.fasta", "text/plain", use_container_width=True)
    st.download_button("Phylogenetic tree (Newick)", analysis["newick"].encode("utf-8"), "bactrace_tree.newick", "text/plain", use_container_width=True)
    st.markdown("### PDF report")
    try:
        pdf = make_pdf_report(analysis)
        st.download_button("Generate / download PDF report", pdf, "bactrace_analysis_report.pdf", "application/pdf", type="primary", use_container_width=True)
    except Exception as exc:
        st.error(f"PDF report unavailable: {exc}")
        st.caption("Ensure reportlab is included in requirements.txt.")
    st.markdown("### Important interpretation notes")
    st.write("- Results are calculated from the current input dataset.")
    st.write("- Species names are inferred from FASTA headers; verify labels against the original database records.")
    st.write("- Alignment and tree methods are identified explicitly in the report and interface.")
    st.write("- This is a research and educational tool, not a validated clinical diagnostic application.")
