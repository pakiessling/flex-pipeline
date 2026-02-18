import json
import sys

with open("./cytetype_annotation.json", "r") as f:
    data = json.load(f)
    
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Deep Cell Annotation Report</title>
    <style>
        :root {
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #333333;
            --text-light: #666666;
            --border: #e9ecef;
            --success: #28a745;
            --warning: #ffc107;
            --danger: #dc3545;
            --info: #17a2b8;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }
        
        body {
            font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
            line-height: 1.6;
            color: var(--text);
            background: var(--bg);
            padding: 20px;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
        }

        /* Headers */
        h1 { color: var(--primary); margin-bottom: 30px; border-bottom: 3px solid var(--accent); padding-bottom: 10px; }
        h2 { color: var(--primary); margin-top: 30px; font-size: 1.5rem; }
        h3 { color: var(--primary); font-size: 1.1rem; margin-top: 15px; margin-bottom: 8px; font-weight: 600; }
        h4 { color: var(--text-light); font-size: 0.95rem; text-transform: uppercase; letter-spacing: 0.5px; margin-top: 15px; }

        /* Cluster Card */
        .cluster-card {
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            margin-bottom: 20px;
            border: 1px solid var(--border);
            overflow: hidden;
        }

        /* Collapsible Header */
        .card-header {
            background: var(--primary);
            color: white;
            padding: 15px 20px;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.2s;
        }
        .card-header:hover { background: #34495e; }
        .card-header .title { font-size: 1.1rem; font-weight: 600; }
        .card-header .id-badge {
            background: rgba(255,255,255,0.2);
            padding: 4px 10px;
            border-radius: 4px;
            margin-right: 15px;
            font-family: monospace;
        }

        .card-body {
            padding: 20px;
            display: none; /* Hidden by default */
        }
        .card-body.show { display: block; }

        /* Summary Grid */
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
            background: #f1f8ff;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid var(--accent);
        }

        .meta-item { margin-bottom: 5px; }
        .meta-label { font-weight: 600; color: var(--text-light); font-size: 0.9em; }

        /* Evidence Sections */
        .evidence-section {
            border-top: 1px solid var(--border);
            margin-top: 20px;
            padding-top: 20px;
        }

        .evidence-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }
        @media (max-width: 900px) { .evidence-grid { grid-template-columns: 1fr; } }

        .gene-box {
            background: #fff;
            border: 1px solid var(--border);
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 10px;
        }
        .gene-box.supporting { border-left: 4px solid var(--success); }
        .gene-box.missing { border-left: 4px solid var(--warning); }
        .gene-box.unexpected { border-left: 4px solid var(--danger); }

        .gene-list { margin-bottom: 8px; }
        .gene-tag {
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            background: #e9ecef;
            font-size: 0.85em;
            font-weight: 600;
            margin-right: 5px;
            margin-bottom: 5px;
        }
        .supporting .gene-tag { background: #d4edda; color: #155724; }
        .missing .gene-tag { background: #fff3cd; color: #856404; }
        .unexpected .gene-tag { background: #f8d7da; color: #721c24; }

        .reasoning-text {
            font-size: 0.95em;
            color: #555;
            font-style: italic;
        }

        /* Metadata & Context */
        .context-box {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 10px;
        }

        /* Review Badges */
        .badge {
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: 600;
            color: white;
            display: inline-block;
        }
        .badge-high { background-color: var(--success); }
        .badge-mod { background-color: var(--warning); color: #333; }
        .badge-low { background-color: var(--danger); }
        .badge-neutral { background-color: var(--text-light); }

        .expert-quote {
            background: #fff;
            border-left: 4px solid var(--info);
            padding: 10px 15px;
            font-style: italic;
            color: #444;
            margin-top: 5px;
        }

        /* Utility */
        .toggle-icon { font-weight: bold; transition: transform 0.3s; }
        .card-header.active .toggle-icon { transform: rotate(180deg); }
        
        code { background: #eee; padding: 2px 5px; border-radius: 3px; font-family: monospace; }
        pre { white-space: pre-wrap; background: #eee; padding: 10px; border-radius: 4px; overflow-x: auto; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Detailed Cell Annotation Report</h1>
        <div id="summary-container"></div>
        <div id="clusters-container"></div>
    </div>

    <script>
        // PASTE YOUR JSON HERE
        const data = DATA_PLACEHOLDER;

        function init() {
            renderSummary(data);
            renderClusters(data);
        }

        function renderSummary(data) {
            const container = document.getElementById('summary-container');
            if(!data.summary || Object.keys(data.summary).length === 0) return;

            let html = `<div class="cluster-card"><div class="card-header active" onclick="toggleCard(this)">`;
            html += `<span class="title">Dataset Summary</span><span class="toggle-icon">▼</span></div>`;
            html += `<div class="card-body show">`;
            // Add summary content here if available in JSON
            html += `<p><em>See individual cluster details below.</em></p>`;
            html += `</div></div>`;
            container.innerHTML = html;
        }

        function renderClusters(data) {
            const container = document.getElementById('clusters-container');
            const annotations = data.annotations || [];
            const raw = data.raw_annotations || {};

            let html = '';

            annotations.forEach(basic => {
                const id = basic.clusterId;
                const details = raw[id]?.latest || {};
                const fullOut = details.annotation?.fullOutput || {};
                const cellType = fullOut.cellType || {};
                const review = details.review || {};
                
                // Determine Confidence Badge
                const conf = review.confidence || 'Unknown';
                let confClass = 'badge-neutral';
                if(conf === 'High') confClass = 'badge-high';
                if(conf === 'Moderate') confClass = 'badge-mod';
                if(conf === 'Low') confClass = 'badge-low';

                html += `
                <div class="cluster-card">
                    <div class="card-header" onclick="toggleCard(this)">
                        <div>
                            <span class="id-badge">ID: ${id}</span>
                            <span class="title">${basic.annotation}</span>
                        </div>
                        <div>
                            <span class="badge ${confClass}">Confidence: ${conf}</span>
                            <span class="toggle-icon" style="margin-left:15px">▼</span>
                        </div>
                    </div>
                    
                    <div class="card-body">
                        <!-- Top Summary Grid -->
                        <div class="summary-grid">
                            <div>
                                <div class="meta-item"><span class="meta-label">Ontology:</span> ${basic.ontologyTerm} (${basic.ontologyTermID})</div>
                                <div class="meta-item"><span class="meta-label">Cell State:</span> ${basic.cellState || 'N/A'}</div>
                                <div class="meta-item"><span class="meta-label">Heterogeneity:</span> ${review.isHeterogeneous ? 'Heterogeneous' : 'Homogeneous'}</div>
                            </div>
                            <div>
                                <div class="meta-item"><span class="meta-label">Granular Description:</span></div>
                                <p style="font-size:0.95em">${basic.granularAnnotation}</p>
                            </div>
                        </div>

                        <!-- Gene Evidence -->
                        <h3>Molecular Evidence</h3>
                        <div class="evidence-grid">
                            <div class="gene-box supporting">
                                <h4>Key Supporting Genes</h4>
                                <div class="gene-list">
                                    ${(cellType.keySupportingGenes || []).map(g => `<span class="gene-tag">${g}</span>`).join('')}
                                </div>
                                <p class="reasoning-text">${cellType.keySupportingGenesReasoning || 'No reasoning provided.'}</p>
                            </div>

                            <div>
                                <div class="gene-box missing">
                                    <h4>Missing Expected Genes</h4>
                                    <div class="gene-list">
                                        ${(cellType.missingGenes || []).map(g => `<span class="gene-tag">${g}</span>`).join('')}
                                    </div>
                                    <p class="reasoning-text">${cellType.missingGenesReasoning || 'None.'}</p>
                                </div>
                                <div class="gene-box unexpected">
                                    <h4>Unexpected Genes</h4>
                                    <div class="gene-list">
                                        ${(cellType.unexpectedGenes || []).map(g => `<span class="gene-tag">${g}</span>`).join('')}
                                    </div>
                                    <p class="reasoning-text">${cellType.unexpectedGenesReasoning || 'None.'}</p>
                                </div>
                            </div>
                        </div>

                        <!-- Pathway & Context -->
                        <div class="evidence-section">
                            <h3>Biological Context & Pathways</h3>
                            <div class="context-box">
                                <div class="meta-label">Pathway Evidence:</div>
                                <p class="reasoning-text" style="margin-bottom:10px">${cellType.pathwaysSupportingEvidenceReasoning || 'N/A'}</p>
                                
                                <div class="meta-label">Tissue Context Fit:</div>
                                <p class="reasoning-text" style="margin-bottom:10px">${cellType.tissueContextFit || 'N/A'}</p>
                                
                                <div class="meta-label">Sample Composition:</div>
                                <p class="reasoning-text">${cellType.clusterContextFit || 'N/A'}</p>
                                ${(details.annotation?.clusterContext?.metadataKeywords || []).length > 0 
                                    ? `<div style="margin-top:5px">${details.annotation.clusterContext.metadataKeywords.map(k => `<span class="badge badge-neutral" style="margin-right:5px; color:black; background:#e0e0e0">${k}</span>`).join('')}</div>` 
                                    : ''}
                            </div>
                        </div>

                        <!-- Expert Review & Reasoning -->
                        <div class="evidence-section">
                            <h3>Expert Review & Logic</h3>
                            
                            <div style="margin-bottom: 15px;">
                                <span class="meta-label">Final Conclusion:</span>
                                <p>${cellType.conclusion || basic.justification}</p>
                            </div>

                            ${renderExpertAssessments(review.fullOutput?.expertAssessments)}
                            
                            ${fullOut.refinedCandidatesWithReasoning ? renderRefinements(fullOut.refinedCandidatesWithReasoning) : ''}
                        </div>
                    </div>
                </div>
                `;
            });

            container.innerHTML = html;
        }

        function renderExpertAssessments(assessments) {
            if (!assessments || assessments.length === 0) return '';
            let html = '<h4>Expert Assessments</h4>';
            assessments.forEach(expert => {
                html += `
                <div class="expert-quote">
                    <strong>${expert.expertType}:</strong> 
                    ${expert.conclusion || expert.biologicalInterpretation}
                    ${expert.alternatives ? `<br><span style="font-size:0.9em; color:#666">Alternative: ${expert.alternatives}</span>` : ''}
                </div>`;
            });
            return html;
        }

        function renderRefinements(candidates) {
            if (!candidates || candidates.length === 0) return '';
            let html = '<div style="margin-top:15px"><span class="meta-label">Candidate Refinement:</span><ul style="margin-left:20px; color:#555">';
            candidates.forEach(cand => {
                if(cand.label && cand.reasoning) {
                    html += `<li><strong>${cand.label}:</strong> ${cand.reasoning}</li>`;
                }
            });
            html += '</ul></div>';
            return html;
        }

        // Accordion Toggle
        window.toggleCard = function(header) {
            header.classList.toggle('active');
            const body = header.nextElementSibling;
            if (body.classList.contains('show')) {
                body.classList.remove('show');
            } else {
                body.classList.add('show');
            }
        }

        // Initialize
        init();
    </script>
</body>
</html>"""

json_data = json.dumps(data, ensure_ascii=False)
final_html = html_template.replace("DATA_PLACEHOLDER", json_data)

with open("report.html", "w", encoding="utf-8") as f:
    f.write(final_html)

print("Report saved as report.html")