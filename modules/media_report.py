import matplotlib.pyplot as plt
import textwrap
import ast
import numpy as np

LETTER_WIDTH_INCH = 8.5
LETTER_HEIGHT_INCH = 11
LEFT_MARGIN = 0.05
RIGHT_MARGIN = 0.95
TOP_MARGIN = 0.95
BOTTOM_MARGIN = 0.05
LINE_HEIGHT = 0.03  # fraction per line
MAX_CHARS_PER_LINE = 90

def wrap_lines(text_lines, width=MAX_CHARS_PER_LINE):
    wrapped = []
    for line in text_lines:
        wrapped.extend(textwrap.wrap(line, width=width))
    return wrapped

def paginate_lines(title, text_lines, bold_flags=None):
    if bold_flags is None:
        bold_flags = [False]*len(text_lines)

    figures = []
    wrapped_lines = []
    wrapped_bold_flags = []

    for i, line in enumerate(text_lines):
        line_wrapped = textwrap.wrap(line, width=MAX_CHARS_PER_LINE)
        wrapped_lines.extend(line_wrapped)
        wrapped_bold_flags.extend([bold_flags[i]]*len(line_wrapped))

    max_lines_per_page = int((TOP_MARGIN - BOTTOM_MARGIN) / LINE_HEIGHT)
    n_pages = (len(wrapped_lines) - 1) // max_lines_per_page + 1

    for i in range(n_pages):
        fig, ax = plt.subplots(figsize=(LETTER_WIDTH_INCH, LETTER_HEIGHT_INCH))
        ax.axis('off')
        start = i*max_lines_per_page
        end = min((i+1)*max_lines_per_page, len(wrapped_lines))
        y = TOP_MARGIN
        if i == 0 and title:
            fig.text(LEFT_MARGIN, y, title, fontsize=14, weight='bold')
            y -= LINE_HEIGHT
        for j in range(start, end):
            weight = 'bold' if wrapped_bold_flags[j] else 'normal'
            fig.text(LEFT_MARGIN, y, wrapped_lines[j], fontsize=10, weight=weight)
            y -= LINE_HEIGHT
        figures.append(fig)
    return figures

def generate_facet_figs(neighbors, knn_percentile, lof_percentile, df_sorted, prediction):
    facet_figs = []

    # kNN metrics
    if knn_percentile < 80:
        knn_label = "In-distribution (typical global distance)"
    elif knn_percentile < 95:
        knn_label = "Extrapolative (globally distant)"
    else:
        knn_label = "Out-of-distribution (very distant)"

    if lof_percentile < 80:
        lof_label = "Normal local density"
    elif lof_percentile < 95:
        lof_label = "Locally sparse"
    else:
        lof_label = "Strong local outlier"

    if knn_percentile > 95 and lof_percentile > 95:
        combined_label = "Strong outlier: globally distant and locally isolated"
    elif knn_percentile > 95 and lof_percentile <= 80:
        combined_label = "Distinct but coherent metabolic niche"
    elif knn_percentile <= 80 and lof_percentile > 95:
        combined_label = "Boundary case between metabolic clusters"
    else:
        combined_label = "Consistent with training distribution"

    knn_lines = [
        f"kNN distance percentile: {knn_percentile:.1f}% → {knn_label}",
        f"LOF percentile: {lof_percentile:.1f}% → {lof_label}",
        f"Overall assessment: {combined_label}"
    ]

    # kNN metrics and neighbors table
    neighbors_table = neighbors.copy()
    fig_table, ax_table = plt.subplots(figsize=(LETTER_WIDTH_INCH, LETTER_HEIGHT_INCH))
    ax_table.axis('off')
    table_data = [neighbors_table.columns.tolist()] + neighbors_table.values.tolist()
    table_artist = ax_table.table(cellText=table_data, cellLoc='left', loc='center')
    table_artist.auto_set_font_size(False)
    table_artist.set_fontsize(10)
    table_artist.scale(1, 1.5)
    for (row, col), cell in table_artist.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold')

    y_start = TOP_MARGIN
    fig_table.text(LEFT_MARGIN, y_start, "kNN Metrics", fontsize=14, weight='bold')
    y = y_start - LINE_HEIGHT
    for line in knn_lines:
        fig_table.text(LEFT_MARGIN, y, line, fontsize=10)
        y -= LINE_HEIGHT

    facet_figs.append(fig_table)

    # Media usage of neighboring taxa:
    df_display = df_sorted[['media_id', 'components', 'occurrence']].copy()
    df_display['components'] = df_display['components'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    media_lines = []
    bold_flags = []

    for _, row in df_display.iterrows():
        media_id = row['media_id']
        components_str = ", ".join(row['components'])
        occurrence_pct = row['occurrence'] * 100
        line = f"{media_id} (occurrence: {occurrence_pct:.1f}% of neighbors): {components_str}"
        media_lines.append(line)
        bold_flags.append(False)

    # Paginate and add media lines
    facet_figs.extend(paginate_lines("Media usage of neighboring taxa:", media_lines, bold_flags))

    # Shared ingredients:
    df_pred = prediction.copy()
    df_pred["components"] = df_pred["components"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )
    df_exploded = df_pred.explode("components")
    n_media = df_pred["media_id"].nunique()
    component_counts = df_exploded.groupby("components")["media_id"].nunique()

    shared_two = component_counts[component_counts >= 2].index.tolist()
    shared_all = component_counts[component_counts == n_media].index.tolist()

    fig, ax = plt.subplots(figsize=(LETTER_WIDTH_INCH, LETTER_HEIGHT_INCH))
    ax.axis("off")

    y = TOP_MARGIN
    fig.text(LEFT_MARGIN, y, "Shared ingredients:", fontsize=14, weight="bold")  # section header
    y -= LINE_HEIGHT * 1.5  # extra space after section header

    fig.text(
        LEFT_MARGIN, y,
        "Components shared by at least two media:",
        fontsize=10, weight="bold"
    )
    y -= LINE_HEIGHT  # first block

    fig.text(
        LEFT_MARGIN + 0.02, y,
        ", ".join(shared_two) if shared_two else "None",
        fontsize=10
    )
    y -= LINE_HEIGHT * 1.5  # space between blocks

    fig.text(
        LEFT_MARGIN, y,
        "Components shared by all media:",
        fontsize=10, weight="bold"
    )
    y -= LINE_HEIGHT  # second block

    fig.text(
        LEFT_MARGIN + 0.02, y,
        ", ".join(shared_all) if shared_all else "None",
        fontsize=10
    )
    y -= LINE_HEIGHT * 1.5  # extra space before info line

    # Add the italic info line at the end of this figure
    info_line = "For more information on media preparation and conditions, visit https://mediadive.dsmz.de/media."
    fig.text(LEFT_MARGIN, y, info_line, fontsize=10, style='italic')

    facet_figs.append(fig)

    return facet_figs


def cofactor_enrichment_barh(facet_figs, counts, baseline="Training", title="Metal cofactor enrichment compared to neighbors"):
    # Pivot counts for easier plotting
    df_pivot = counts.pivot(index="CofactorFinal", columns="Set", values="ra").fillna(0)
    
    cofactors = df_pivot.index.tolist()
    
    # Compute deviations from baseline
    baseline_values = df_pivot.get(baseline, 0)
    test_values = df_pivot.get("Test", 0)
    enrichment = test_values - baseline_values  # positive = test uses more

    y_pos = np.arange(len(cofactors))

    # Create figure
    fig, ax = plt.subplots(figsize=(LETTER_WIDTH_INCH, LETTER_HEIGHT_INCH))
    
    colors = ['steelblue' if x >= 0 else 'tomato' for x in enrichment]
    bars = ax.barh(y_pos, enrichment, color=colors)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(cofactors)
    ax.axvline(0, color='black', linewidth=1)  # baseline/zero line
    ax.set_xlabel(f"Deviation from {baseline} (%)")
    ax.set_title(title, fontsize=14, pad=20)
    
    # Annotate bars with percentages
    for bar, val in zip(bars, enrichment):
        width = bar.get_width()
        x_pos = width - 0.01 if val >= 0 else width + 0.01
        ha = 'right' if val >= 0 else 'left'
        ax.text(x_pos, bar.get_y() + bar.get_height()/2, f"{val*100:+.1f}%", 
                va='center', ha=ha, color='white', fontsize=10, fontweight='bold')

    ax.invert_yaxis()  # largest on top
    facet_figs.append(fig)


