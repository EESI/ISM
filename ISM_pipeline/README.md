# ISM pipeline on [Proteus](https://proteusmaster.urcf.drexel.edu/urcfwiki/index.php/Main_Page)

#### How to run:
1. Sequence data (as well as metadata) should be downloaded from [GISAID](https://www.gisaid.org/) website (nextfasta and nextmeta file).
2. One step must be completed before the execution of the next step. Run `pipeline_0_filter.sh` then `pipeline_1_mafft.sh` and finally `pipeline_2_run_ISM.py`.
3. Modification of those bash scripts are needed (change the email address, specific file path etc), read the comments in the scripts carefully.

#### Output explanation:
1. JSON files are saved to create the pie charts/time series plots in the [jupyter-notebook-report](https://github.com/EESI/ISM/blob/master/ISM-report-20200515-with_error_correction.ipynb) (e.g., `region_time_series.json`, `state_pie_chart.json` and `region_pie_chart.json`).
2. It is required to acknowledge GISAID and labs involved. The acknowledgement table is in `acknowledgement_table.csv`
3. The ISM table along with patient information are saved to `IMS_17nt_with_correction.csv` for other customized visualization.
4. The gene annotation of ISM positions are saved in `ISM_gene_in_reference.csv`.
5. `web_logo.json` is for the weblogo figure, which gives an intuitive explanation of ISMs.

#### suggestions:
1. New sequences are submitted to GISAID database, your analysis/visualization should consider future updates. Therefore, don't hardcode region numbers, color numbers in your framework. Instead, you should dynamically figure it out from the data.
2. Color map is another big challenge. Think about how to color code ISMs in a meaningful way and use the same ISM color map globally.

