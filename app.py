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
        fig.add_trace(go.Scatter(x=[x0, x1], y=[y0, y1], mode="lines",
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
