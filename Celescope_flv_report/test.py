import pandas as pd
import json
import streamlit as st
from PIL import Image
import plotly.graph_objs as go
import plotly.offline as pltoff
import plotly.express as px
import numpy as np
import collections
import math

BarcodeRankPlotSegment = collections.namedtuple('BarcodeRankPlotSegment', ['start', 'end', 'cell_density', 'legend'])
BC_PLOT_COLORS = ['#dddddd', '#d1d8dc', '#c6d3dc', '#bacfdb', '#aecada', '#a3c5d9', '#97c0d9', '#8cbbd8', '#80b7d7',
                  '#74b2d7', '#6aadd6', '#66abd4', '#62a8d2', '#5ea5d1', '#59a2cf', '#559fce', '#519ccc', '#4d99ca',
                  '#4997c9', '#4594c7', '#4191c5', '#3d8dc4', '#3a8ac2', '#3787c0', '#3383be', '#3080bd', '#2c7cbb',
                  '#2979b9', '#2676b7', '#2272b6', '#1f6eb3', '#1d6ab0', '#1a65ac', '#1861a9', '#155ca6', '#1358a2',
                  '#10539f', '#0e4f9b', '#0b4a98', '#094695', '#09438f', '#0a4189', '#0c3f83', '#0d3d7c', '#0e3b76',
                  '#103970', '#11366a', '#123463', '#14325d', '#153057']

def segment_log_plot_by_length(y_data, x_start, x_end):
    """
    Given the extends of the mixed region [x_start, x_end), compute
    the x-indices that would divide the plot into segments of a
    prescribed length (in pixel coordinates) with a minimum number
    of barcodes in each segment
    """

    if x_end <= x_start:
        return []

    SEGMENT_NORMALIZED_MAX_LEN = 0.02
    MIN_X_SPAN = 20

    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))

    this_segment_len = 0.0
    segment_idx = [x_start]

    for i in range(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        this_segment_len += np.linalg.norm([dx, dy])
        if this_segment_len >= SEGMENT_NORMALIZED_MAX_LEN and i > (segment_idx[-1] + MIN_X_SPAN):
            segment_idx.append(i+1)
            this_segment_len = 0.0

    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)

    return segment_idx


def get_plot_segment(start_index, end_index, sorted_bc, cell_barcodes, legend=False):
    """
    Helper function to build a plot segment.
    """
    assert end_index > start_index
    num_cells = sum([1 for i in range(start_index, end_index) if sorted_bc[i] in cell_barcodes])
    density = float(num_cells)/float(end_index-start_index)
    return BarcodeRankPlotSegment(start=start_index, end=end_index, cell_density=density, legend=legend)


def counter_barcode_rank_plot_data(count_data_path):
    """
    get cell density for each plot_segments
    :param count_data_path:
    :return: sorted_counts, plot_segments, cell_nums
    """
    count_data = pd.read_csv(count_data_path, index_col=0, sep='\t')
    cell_bc = np.array(count_data[count_data['mark'] == 'CB'].index)
    sorted_bc = np.array(count_data.index)
    sorted_counts = np.array(count_data['UMI'])
    cell_nums = len(cell_bc)
    total_bc = len(sorted_bc)
    # find the first barcode which is not a cell
    first_non_cell = total_bc
    for i, bc in enumerate(sorted_bc):
        if bc not in cell_bc:
            first_non_cell = i
            break

    # find the last barcode which is a cell
    last_cell = 0
    for i in reversed(range(total_bc)):
        if sorted_bc[i] in cell_bc:
            last_cell = i
            break

    ranges = [0, first_non_cell, last_cell+1, total_bc]
    plot_segments = []
    plot_segments.append(BarcodeRankPlotSegment(start=0, end=ranges[1], cell_density=1.0, legend=True))
    plot_segments.append(BarcodeRankPlotSegment(start=ranges[2], end=ranges[3], cell_density=0.0, legend=True))

    mixed_segments = segment_log_plot_by_length(sorted_counts, ranges[1], ranges[2])
    for i in range(len(mixed_segments) - 1):
        plot_segments.append(
            get_plot_segment(mixed_segments[i], mixed_segments[i + 1], sorted_bc, cell_bc, legend=False))

    return sorted_counts, plot_segments, cell_nums


def get_plot_data(plot_segments, counts):
    plot_data = []
    for segment in plot_segments:
        plot_data.append(build_plot_data_dict(segment, counts))

    return plot_data


def convert_numpy_array_to_line_chart(array, ntype):
    array = np.sort(array)[::-1]

    rows = []
    previous_count = None
    for (index,), count in np.ndenumerate(array):
        if index == 0 or index == len(array)-1:
            rows.append([index, ntype(count)])
        elif previous_count != count:
            previous_index = rows[-1][0]
            if previous_index != index - 1:
                rows.append([index - 1, ntype(previous_count)])
            rows.append([index, ntype(count)])
        previous_count = count
    return rows


def build_plot_data_dict(plot_segment, counts):
    """
    Construct the data for a plot segment by appropriately slicing the
    counts
    Inputs:
    - plot_segment: BarcodeRankPlotSegment containing [start, end)
        of the segment, the cell density and legend visibility option
    - counts: Reverse sorted UMI counts for all barcodes.
    """

    start = max(0, plot_segment.start - 1)  # -1 for continuity between two charts
    end = plot_segment.end
    plot_rows = convert_numpy_array_to_line_chart(counts[start:end], int)
    name = 'Cells' if plot_segment.cell_density > 0 else 'Background'

    # Setup the tooltip
    if plot_segment.cell_density > 0.:
        n_barcodes = plot_segment.end - plot_segment.start
        n_cells = int(round(plot_segment.cell_density * n_barcodes))
        hover = "{:.0f}% Cells<br>({}/{})".format(100 * plot_segment.cell_density, n_cells, n_barcodes)
    else:
        hover = "Background"

    data_dict = {
        "x": [],
        "y": [],
        "name": name,
        "hoverinfo": "text",
        "text": hover,
        "type": "scattergl",
        "mode": "lines",
        "line": {
            "width": 3,
            "color": BC_PLOT_CMAP(plot_segment.cell_density),
        },
        "showlegend": plot_segment.legend,
    }
    offset = 1 + start  # it's a log-log plot, hence the 1
    for index, count in plot_rows:
        data_dict["x"].append(index + offset)
        data_dict["y"].append(count)

    # Handle case where the data is empty
    if len(data_dict["x"]) == 0:
        data_dict["x"].append(0)
        data_dict["y"].append(0)

    return data_dict

def BC_PLOT_CMAP(density):
    """
    Colormap utility fn to map a number to one of the colors in the gradient
    color scheme defined above
    Input
    - density : A real number in the range [0,1]
    """
    assert density >= 0.
    assert density <= 1.
    levels = len(BC_PLOT_COLORS)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return BC_PLOT_COLORS[ind]


def plot_barcode_rank(count_file_path):
    sorted_counts, plot_segments, _cell_nums = counter_barcode_rank_plot_data(count_file_path)
    plot_data = get_plot_data(plot_segments, sorted_counts)

    plotly_data = [go.Scatter(x=dat['x'], y=dat['y'], name=dat['name'], mode=dat['mode'], showlegend=dat['showlegend'],
                              marker={'color': dat['line']['color']}, line=dat['line'], text=dat['text']) for dat in
                   plot_data]

    layout = go.Layout(width=470, height=313,
                       title={"text": "Barcode rank", "font": {"color": "black"}, "x": 0.5},
                       xaxis={"type": "log", "title": "Barcode", "titlefont": {"color": "black"},
                              "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                       yaxis={"type": "log", "title": "UMI counts", "titlefont": {"color": "black"},
                              "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                       margin=dict(l=50, r=0, t=30, b=30),
                       plot_bgcolor="#FFFFFF")

    config = dict({"displayModeBar": True,
                   "staticPlot": False,
                   "showAxisDragHandles": False,
                   "modeBarButtons": [["toImage", "resetScale2d"]],
                   "scrollZoom": False,
                   "displaylogo": False})

    fig = go.Figure(data=plotly_data, layout=layout)

    chart = pltoff.plot(fig, include_plotlyjs=True, output_type='div', config=config)

    return chart


st.title("Celescope-1.14.0 Report")
page_text = """

CeleScope is a collection of bioinfomatics analysis pipelines to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

celescope flv_CR for Single cell full Length TCR/BCR data generated with sCircleTM Kits. It performs preprocessing, barcode conversion, assemble, annotation and filtering.

GitHub:

https://github.com/singleron-RD/CeleScope
"""
image = Image.open("C:/Users/admin/Desktop/百英/D37/logo2.png")
st.sidebar.image(image)
st.sidebar.markdown(page_text)


with open("C:/Users/admin/Desktop/百英/D37/metrics.json") as f:
    load_dict = json.load(f)

st.header("Sample")
with st.expander("See explanation"):
    st.write("""
        Chemistry : For more information, see 
        
        https://github.com/singleron-RD/CeleScope/blob/master/docs/chemistry.md.
    """)
sample_summary = load_dict["sample_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")

st.header("Demultiplexing")
with st.expander("See explanation"):
    st.write("""
        Raw Reads : Total reads from FASTQ files.
        
        Valid Reads : Reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads).
        
        Q30 of Barcodes : Percent of barcode base pairs with quality scores over Q30.
        
        Q30 of UMIs : Percent of UMI base pairs with quality scores over Q30.
    """)
barcode_summary = load_dict["barcode_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")

st.header("Mapping")
with st.expander("See explanation"):
    st.write("""
        Reads Mapped To Any V(D)J Gene : Fraction of reads that partially or wholly map to any germline V(D)J gene segment.
        
        Reads Mapped To IGH : Fraction of reads that map partially or wholly to a germline IGH gene segment.
        
        Reads Mapped To IGL : Fraction of reads that map partially or wholly to a germline IGL gene segment.
        
        Reads Mapped To IGK : Fraction of reads that map partially or wholly to a germline IGK gene segment.
    """)
sample_summary = load_dict["mapping_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")

st.header("Cells")
with st.expander("See explanation"):
    st.write("""
        Estimated Number of Cells : The number of barcodes estimated to be associated with cells that express targeted V(D)J transcripts.

        Fraction Reads in Cells : Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes.

        Mean Reads per Cell : Number of input read pairs divided by the estimated number of cells.

        Mean Used Reads per Cell : Mean number of read pairs used in assembly per cell-associated barcode.

        Median Used IGH UMIs per Cell : Median number of UMIs assigned to a IGH contig per cell.

        Median Used IGL UMIs per Cell : Median number of UMIs assigned to a IGL contig per cell.

        Median Used IGK UMIs per Cell : Median number of UMIs assigned to a IGK contig per cell.
    """)
sample_summary = load_dict["cells_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")

count_file = "C:/Users/admin/Desktop/百英/D37/count.txt"
chart=plot_barcode_rank(count_file)
st.plotly_chart(chart)

st.header("Annotation")
with st.expander("See explanation"):
    st.write("""
        Cells With Productive V-J Spanning Pair : Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair.A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.

        Cells With Productive V-J Spanning (IGK, IGH) Pair : Fraction of cell-associated barcodes with at least one productive contig for each chain of the (IGK, IGH) receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.

        Cells With Productive V-J Spanning (IGL, IGH) Pair : Fraction of cell-associated barcodes with at least one productive contig for each chain of the (IGL, IGH) receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.

        Cells With IGH Contig : Fraction of cell-associated barcodes with at least one IGH contig annotated as a full or partial V(D)J gene.

        Cells With CDR3-annotated IGH Contig : Fraction of cell-associated barcodes with at least one IGH contig where a CDR3 was detected.

        Cells With V-J Spanning IGH Contig : Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for IGH.

        Cells With Productive IGH Contig : Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for IGH, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.

        Cells With IGL Contig : Fraction of cell-associated barcodes with at least one IGL contig annotated as a full or partial V(D)J gene.

        Cells With CDR3-annotated IGL Contig : Fraction of cell-associated barcodes with at least one IGL contig where a CDR3 was detected.

        Cells With V-J Spanning IGL Contig : Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for IGL.

        Cells With Productive IGL Contig : Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for IGL, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.

        Cells With IGK Contig : Fraction of cell-associated barcodes with at least one IGK contig annotated as a full or partial V(D)J gene.

        Cells With CDR3-annotated IGK Contig : Fraction of cell-associated barcodes with at least one IGK contig where a CDR3 was detected.

        Cells With V-J Spanning IGK Contig : Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for IGK.

        Cells With Productive IGK Contig : Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for IGK, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.
    """)
sample_summary = load_dict["annotation_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")

st.header("Clonotypes")
clonotype_file = pd.read_csv("C:/Users/admin/Desktop/百英/D37/clonotypes.csv")
clonotype_file['ClonotypeID'] = clonotype_file['clonotype_id'].apply(lambda x: x.strip('clonetype'))
clonotype_file['Frequency'] = clonotype_file['frequency']
clonotype_file['Proportion'] = clonotype_file['proportion'].apply(lambda x: f'{round(x * 100, 2)}%')
clonotype_file['CDR3_aa'] = clonotype_file['cdr3s_aa'].apply(lambda x: x.replace(';', ' ; '))

df_table = clonotype_file[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
st.write(df_table)

clonotype_file['ClonotypeID'] = clonotype_file['ClonotypeID'].astype("int")
clonotype_file.sort_values(by=['ClonotypeID'], inplace=True)

df = pd.DataFrame({"Clonotype ID": [str(i) for i in list(clonotype_file.head(10).ClonotypeID)],
                   "Fraction of Cells": clonotype_file.head(10).proportion.tolist()})
fig = px.bar(df, x="Clonotype ID", y="Fraction of Cells", title="Top 10 Clonotype Frequencies", color_discrete_sequence=["#9EE6CF"], labels={'x': 'Clonotype ID', 'y': 'Fraction of Cells'})
fig.update_layout(xaxis = dict(dtick = 1))
st.plotly_chart(fig, theme="streamlit", use_container_width=True)


