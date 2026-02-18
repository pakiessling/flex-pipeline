"""
07_report.py — Generate a self-contained HTML summary report.

Embeds all plots as base64 so the report is a single portable .html file.
Sections:
  1. Run summary (config, software versions, cell counts)
  2. Per-sample QC (cell counts, doublet rates)
  3. Integration overview (UMAP by sample + cluster)
  4. CyteType annotation (if JSON available)
  5. Top marker genes per cluster (table)
"""

import argparse
import base64
import datetime
import io
import json
import logging
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _img_tag(b64: str, alt: str = "", style: str = "max-width:100%;") -> str:
    return f'<img src="data:image/png;base64,{b64}" alt="{alt}" style="{style}">'


def _df_to_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    return df.head(max_rows).to_html(
        classes="table", border=0, index=True, escape=True
    )


# ---------------------------------------------------------------------------
# Plot generators
# ---------------------------------------------------------------------------

def _plot_umap_by(adata, color_col: str, title: str):
    if color_col not in adata.obs.columns and color_col not in adata.var_names:
        return None
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=color_col, ax=ax, show=False, title=title)
    return fig


def _plot_qc_violin(adata):
    metrics = [c for c in ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_MT"]
               if c in adata.obs.columns]
    if not metrics:
        return None
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics):
        sc.pl.violin(adata, keys=m, groupby="Sample", ax=ax, show=False, rotation=45)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _section_summary(adata) -> str:
    pipeline_log = adata.uns.get("pipeline_log", {})
    rows = [
        ("Total cells", adata.n_obs),
        ("Total genes", adata.n_vars),
        ("Samples", adata.obs["Sample"].nunique() if "Sample" in adata.obs else "?"),
        ("Report generated", datetime.datetime.now().strftime("%Y-%m-%d %H:%M")),
    ]
    rows_html = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    sw_sections = ""
    for step, info in pipeline_log.items():
        sw = info.get("software", {})
        if sw:
            sw_str = ", ".join(f"{k}={v}" for k, v in sw.items())
            sw_sections += f"<tr><td>{step}</td><td>{sw_str}</td></tr>"
    return f"""
    <section id="summary">
      <h2>Run Summary</h2>
      <table class="table"><tbody>{rows_html}</tbody></table>
      {"<h3>Software Versions</h3><table class='table'><thead><tr><th>Step</th><th>Versions</th></tr></thead><tbody>" + sw_sections + "</tbody></table>" if sw_sections else ""}
    </section>
    """


def _section_qc(adata) -> str:
    if "Sample" not in adata.obs.columns:
        return "<section id='qc'><h2>QC</h2><p>No sample metadata found.</p></section>"

    stats = []
    for sample in sorted(adata.obs["Sample"].unique()):
        sub = adata[adata.obs["Sample"] == sample]
        row = {"Sample": sample, "Cells": sub.n_obs}
        if "cell_quality" in sub.obs.columns:
            n_low = (sub.obs["cell_quality"] == "low-quality").sum()
            row["Low-quality (%)"] = f"{100 * n_low / sub.n_obs:.1f}"
        if "predicted_doublet" in sub.obs.columns:
            n_dbl = sub.obs["predicted_doublet"].astype(bool).sum()
            row["Doublets (%)"] = f"{100 * n_dbl / sub.n_obs:.1f}"
        stats.append(row)

    table_html = pd.DataFrame(stats).to_html(classes="table", border=0, index=False)

    # QC violin
    fig = _plot_qc_violin(adata)
    violin_html = _img_tag(_fig_to_b64(fig), "QC metrics by sample") if fig else ""

    return f"""
    <section id="qc">
      <h2>Per-sample QC</h2>
      {table_html}
      {violin_html}
    </section>
    """


def _section_integration(adata) -> str:
    imgs = []
    for col in ["Sample", "leiden_1_5", "leiden_3", "singler_label", "phase"]:
        if col in adata.obs.columns:
            fig = _plot_umap_by(adata, col, col.replace("_", " ").title())
            if fig:
                imgs.append(_img_tag(_fig_to_b64(fig), col, style="width:48%; margin:1%;"))

    imgs_html = "\n".join(imgs)
    return f"""
    <section id="integration">
      <h2>Integration &amp; Clustering</h2>
      <div style="display:flex; flex-wrap:wrap;">
        {imgs_html}
      </div>
    </section>
    """


def _section_cytetype(cytetype_json_path: str) -> str:
    if not cytetype_json_path or not os.path.exists(cytetype_json_path):
        return ""

    with open(cytetype_json_path) as fh:
        data = json.load(fh)

    json_str = json.dumps(data, ensure_ascii=False)

    # Inline the cytetype_report.py rendering logic
    cards_js = """
    function renderClusters(data) {
        const container = document.getElementById('cytetype-container');
        const annotations = data.annotations || [];
        const raw = data.raw_annotations || {};
        let html = '';
        annotations.forEach(basic => {
            const id = basic.clusterId;
            const details = raw[id]?.latest || {};
            const fullOut = details.annotation?.fullOutput || {};
            const cellType = fullOut.cellType || {};
            const review = details.review || {};
            const conf = review.confidence || 'Unknown';
            let confClass = 'badge-neutral';
            if(conf === 'High') confClass = 'badge-high';
            if(conf === 'Moderate') confClass = 'badge-mod';
            if(conf === 'Low') confClass = 'badge-low';
            html += `<div class="cluster-card">
              <div class="card-header" onclick="this.classList.toggle('active');this.nextElementSibling.classList.toggle('show')">
                <div><span class="id-badge">ID: ${id}</span><span style="font-weight:600">${basic.annotation}</span></div>
                <span class="badge ${confClass}" style="margin-left:10px">Confidence: ${conf}</span>
              </div>
              <div class="card-body">
                <p><strong>Ontology:</strong> ${basic.ontologyTerm || ''} (${basic.ontologyTermID || ''})</p>
                <p><strong>Granular annotation:</strong> ${basic.granularAnnotation || ''}</p>
                <p><strong>Conclusion:</strong> ${cellType.conclusion || basic.justification || ''}</p>
                <p><strong>Key supporting genes:</strong> ${(cellType.keySupportingGenes||[]).join(', ')}</p>
                <p><strong>Missing genes:</strong> ${(cellType.missingGenes||[]).join(', ') || 'None'}</p>
                <p><strong>Unexpected genes:</strong> ${(cellType.unexpectedGenes||[]).join(', ') || 'None'}</p>
              </div>
            </div>`;
        });
        container.innerHTML = html;
    }
    const cytedata = DATA_PLACEHOLDER;
    renderClusters(cytedata);
    """

    return f"""
    <section id="cytetype">
      <h2>CyteType Cluster Annotation</h2>
      <div id="cytetype-container"></div>
      <script>{cards_js.replace("DATA_PLACEHOLDER", json_str)}</script>
    </section>
    """


def _section_markers(adata) -> str:
    if "rank_genes_groups" not in adata.uns:
        return ""

    rgg = adata.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)[:20]  # cap at 20 clusters for readability
    group_key = rgg["params"].get("groupby", "cluster")

    cards = []
    for g in groups:
        genes = [str(x) for x in rgg["names"][g] if str(x) not in ("", "nan")][:15]
        if not genes:
            continue
        gene_tags = " ".join(f'<span class="gene-tag">{gene}</span>' for gene in genes)
        cards.append(f"""
        <div class="cluster-card">
          <div class="card-header" onclick="this.classList.toggle('active');this.nextElementSibling.classList.toggle('show')">
            <span class="id-badge">{group_key}: {g}</span>
            <span style="font-weight:600">Top marker genes</span>
          </div>
          <div class="card-body">
            <div style="line-height:2">{gene_tags}</div>
          </div>
        </div>
        """)

    return f"""
    <section id="markers">
      <h2>Marker Genes</h2>
      {"".join(cards)}
    </section>
    """


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_report(adata, cytetype_json_path: str, output_dir: str):
    logger.info("Building HTML report sections …")

    html_parts = {
        "summary":     _section_summary(adata),
        "qc":          _section_qc(adata),
        "integration": _section_integration(adata),
        "cytetype":    _section_cytetype(cytetype_json_path),
        "markers":     _section_markers(adata),
    }

    body = "\n".join(html_parts.values())

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>flex-pipeline Report</title>
  <style>
    :root {{
      --primary: #2c3e50; --accent: #3498db; --bg: #f8f9fa;
      --card-bg: #fff; --text: #333; --border: #e9ecef;
      --success: #28a745; --warning: #ffc107; --danger: #dc3545;
      --info: #17a2b8; --text-light: #666;
    }}
    * {{ margin:0; padding:0; box-sizing:border-box; }}
    body {{ font-family:'Segoe UI',Roboto,Arial,sans-serif; background:var(--bg); color:var(--text); padding:24px; }}
    .container {{ max-width:1400px; margin:0 auto; }}
    nav {{ background:var(--primary); padding:12px 20px; border-radius:6px; margin-bottom:24px; display:flex; gap:20px; flex-wrap:wrap; }}
    nav a {{ color:#fff; text-decoration:none; font-size:0.9rem; opacity:0.85; }}
    nav a:hover {{ opacity:1; }}
    h1 {{ color:var(--primary); border-bottom:3px solid var(--accent); padding-bottom:10px; margin-bottom:24px; }}
    h2 {{ color:var(--primary); margin:32px 0 16px; font-size:1.4rem; border-left:4px solid var(--accent); padding-left:10px; }}
    h3 {{ color:var(--primary); margin:16px 0 8px; font-size:1.1rem; }}
    section {{ margin-bottom:40px; }}
    .table {{ width:100%; border-collapse:collapse; background:var(--card-bg); border-radius:6px; overflow:hidden; box-shadow:0 1px 4px rgba(0,0,0,.06); }}
    .table th, .table td {{ padding:10px 14px; text-align:left; border-bottom:1px solid var(--border); }}
    .table th {{ background:var(--primary); color:#fff; font-weight:600; }}
    .table tr:last-child td {{ border-bottom:none; }}
    .cluster-card {{ background:var(--card-bg); border:1px solid var(--border); border-radius:8px; margin-bottom:12px; overflow:hidden; box-shadow:0 1px 4px rgba(0,0,0,.05); }}
    .card-header {{ background:var(--primary); color:#fff; padding:12px 16px; cursor:pointer; display:flex; justify-content:space-between; align-items:center; transition:background .2s; }}
    .card-header:hover {{ background:#34495e; }}
    .card-body {{ padding:16px; display:none; }}
    .card-body.show {{ display:block; }}
    .id-badge {{ background:rgba(255,255,255,.2); padding:3px 9px; border-radius:4px; margin-right:12px; font-family:monospace; }}
    .badge {{ padding:3px 8px; border-radius:4px; font-size:.8em; font-weight:600; color:#fff; }}
    .badge-high {{ background:var(--success); }}
    .badge-mod  {{ background:var(--warning); color:#333; }}
    .badge-low  {{ background:var(--danger); }}
    .badge-neutral {{ background:var(--text-light); }}
    .gene-tag {{ display:inline-block; padding:2px 8px; border-radius:12px; background:#d4edda; color:#155724; font-size:.82em; font-weight:600; margin:2px; }}
    img {{ border-radius:6px; }}
  </style>
</head>
<body>
  <div class="container">
    <h1>flex-pipeline Analysis Report</h1>
    <nav>
      <a href="#summary">Summary</a>
      <a href="#qc">QC</a>
      <a href="#integration">Integration</a>
      {"<a href='#cytetype'>CyteType</a>" if cytetype_json_path else ""}
      <a href="#markers">Markers</a>
    </nav>
    {body}
  </div>
</body>
</html>"""

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "report.html")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info(f"Report saved → {out_path}")
    return out_path


def main(args):
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input not found: {args.input_file}")

    logger.info(f"Loading {args.input_file} …")
    adata = sc.read_h5ad(args.input_file)
    logger.info(f"  {adata.n_obs} cells × {adata.n_vars} genes")

    build_report(adata, args.cytetype_json, args.output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate HTML pipeline report.")
    parser.add_argument("--input_file",   required=True, help="Final .h5ad path")
    parser.add_argument("--output_dir",   required=True, help="Directory to write report.html")
    parser.add_argument("--cytetype_json", default="", help="Path to cytetype_annotation.json (optional)")
    args = parser.parse_args()
    main(args)
