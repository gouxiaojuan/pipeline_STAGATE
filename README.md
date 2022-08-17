# Compute SRT spatial domain
This is a process of concatenating STAGATE
## Work Flow
This process uses the method of finding spatial domain in STAGATE. 
For the specific reason and process, see [STAGATE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8976049/)
## Requirements
This script runs based on python, and the required python packages are as follows：
* argparse
* numpy
* pandas
* scanpy
* sys
* STAGATE
* matplotlib
* sklearn

For the version requirements of the required packages, see [requirement.txt](https://github.com/gouxiaojuan/pipeline_STAGATE/blob/main/requirement.txt)<br>
For specific software installation, please refer to the [official tutorial of STARGATE](https://stagate.readthedocs.io/en/latest/index.html)<br>
## Use the script
To use this script you can calculate spatial domain of your own data, you need to provide 
a spatial expression matrix with <strong> row names for gene names, columns for cell names </strong>
and a file with spatial location information 
1. Example of a spatial representation matrix: [RNA_counts.csv](https://github.com/gouxiaojuan/pipeline_STAGATE/blob/main/sample/RNA_counts.zip)<br>
2. Example of a file with spatial location information: [position.csv](https://github.com/gouxiaojuan/pipeline_STAGATE/blob/main/sample/position.csv)<br>

If you have these two files ready, you can get the SVG using the following command：<br>
`$ python pipeline_STAGATE.py --counts_df=RNA_counts.csv --location_df=position.csv`

With this script you can get a [spatial_domain.pdf](https://github.com/gouxiaojuan/pipeline_STAGATE/blob/main/sample/spatial_domain.pdf)