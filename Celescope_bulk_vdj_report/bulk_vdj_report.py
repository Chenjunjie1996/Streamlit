import pandas as pd
import json
import streamlit as st
from PIL import Image
import plotly.express as px



st.title("Celescope-1.14.0 Report")

page_text = """

CeleScope is a collection of bioinfomatics analysis pipelines to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

celescope vdj for Single-cell Immune Repertoire data generated with GEXSCOPETM IR kits. It performs preprocessing, UMI consensus, vdj sequence alignment, UMI filtering and clonetypes counting.

GitHub:

https://github.com/singleron-RD/CeleScope
"""
image = Image.open("C:/Users/admin/Desktop/百英/D37/logo2.png")
st.sidebar.image(image)
st.sidebar.markdown(page_text)

with open("C:/Users/admin/Desktop/temp/metrics.json") as f:
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

st.header("Trimming")
with st.expander("See explanation"):
    st.write("""
        Raw Reads : Total reads from FASTQ files.

        Valid Reads : Reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads).

        Q30 of Barcodes : Percent of barcode base pairs with quality scores over Q30.

        Q30 of UMIs : Percent of UMI base pairs with quality scores over Q30.
    """)
sample_summary = load_dict["cutadapt_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")


st.header("Consensus")
with st.expander("See explanation"):
    st.write("""
        UMI Counts : Total UMI from FASTQ files.

        Mean UMI Length : Mean of all UMI length.

        Ambiguous Base Counts : Number of bases that do not pass consensus threshold.

        Filter UMI Counts : Filtered UMI from Consensus UMI fasta.
    """)
sample_summary = load_dict["consensus_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")

st.header("Mapping")
with st.expander("See explanation"):
    st.write("""
        Species : Human or Mouse.

        UMIs Mapped To Any VDJ Gene : UMIs Mapped to any germline VDJ gene segments.

        UMIs with CDR3 : UMIs with CDR3 sequence.

        UMIs with Correct CDR3 : UMIs with CDR3 might have stop codon and these Reads are classified as incorrect.

        UMIs Mapped Confidently To VJ Gene : UMIs with productive rearrangement mapped to VJ gene pairs and with correct CDR3.

        UMIs Mapped To TRA : UMIs mapped confidently to TRA.

        UMIs Mapped To TRB : UMIs mapped confidently to TRB.
    """)
sample_summary = load_dict["mapping_vdj_summary"]
for k, v in sample_summary.items():
    st.write(f"{k} : {v}")


st.header("Clonotypes")
clonotype_file = pd.read_csv("C:/Users/admin/Desktop/temp/Hum_0209PBMC_wx3w_T3lib_clonetypes.csv")
# clonotype_file['ClonotypeID'] = clonotype_file['clonotype_id'].apply(lambda x: x.strip('clonetype'))
# clonotype_file['Frequency'] = clonotype_file['frequency']
# clonotype_file['Proportion'] = clonotype_file['proportion'].apply(lambda x: f'{round(x * 100, 2)}%')
# clonotype_file['CDR3_aa'] = clonotype_file['cdr3s_aa'].apply(lambda x: x.replace(';', ' ; '))

df_table = clonotype_file[['ClonotypeID', 'aaSeqCDR3', 'Frequency', 'Proportion']]
st.write(df_table)

# clonotype_file['ClonotypeID'] = clonotype_file['ClonotypeID'].astype("int")
# clonotype_file.sort_values(by=['ClonotypeID'], inplace=True)

df = pd.DataFrame({"Clonotype ID": [str(i) for i in list(clonotype_file.head(10).ClonotypeID)],
                   "Fraction of Cells": clonotype_file.head(10).Proportion.tolist()})
fig = px.bar(df, x="Clonotype ID", y="Fraction of Cells", title="Top 10 Clonotype Frequencies", color_discrete_sequence=["#9EE6CF"], labels={'x': 'Clonotype ID', 'y': 'Fraction of Cells'})
fig.update_layout(xaxis = dict(dtick = 1))
st.plotly_chart(fig, theme="streamlit", use_container_width=True)