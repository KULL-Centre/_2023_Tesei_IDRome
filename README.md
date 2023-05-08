<a target="_blank" href="https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDRLab.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

# Analyses of conformational ensembles of the human IDRome

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks, and data for reproducing the results presented in the manuscript _Conformational ensembles of the human intrinsically disordered proteome: Bridging chain compaction with function and sequence conservation_.

The CSV file `IDRome_DB.csv` and the Excel Sheet `IDRome_DB.xlsx` list the sequence and various sequence- and conformational properties of all the 29,998 IDRs.

Simulations trajectories and time series of conformational properties are available for all the 29,998 IDRs at [sid.erda.dk/sharelink/AVZAJvJnCO](https://sid.erda.dk/sharelink/AVZAJvJnCO).

We also provide a [Notebook](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDRLab.ipynb) on [Google Colab](https://colab.research.google.com/) to generate conformational ensembles of user-supplied sequences and conditions using the [CALVADOS](https://github.com/KULL-Centre/CALVADOS) model. 

### Layout
- `seq_conf_prop.ipynb` reproduces Fig. 1, 3, S3, S4, S6, and S10
- `go_analysis.ipynb` reproduces Fig. 2
- `conservation_analysis.ipynb` reproduces Fig. 4
- `clinvar_fmug.ipynb` reproduces Fig. 5, S8, and S9
- `uniprot_domains.ipynb` reproduces Fig. S1, S2, and S5
- `svr_model.ipynb` reproduces Fig. S7
- `go_uniprot_calls.ipynb` performs API calls to obtain gene ontology terms from UniProt
- `calc_seq_prop.ipynb` computes sequence descriptors for all the 29,998 IDRs
- `md_simulations/` contains code and data related to single-chain simulations performed using the CALVADOS model and [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/) v2.9.3 installed with [mphowardlab/azplugins](https://github.com/mphowardlab/azplugins)

### Usage

To open the Notebooks, install [Miniconda](https://conda.io/miniconda.html) and make sure all required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate idrome
    jupyter-notebook
```
