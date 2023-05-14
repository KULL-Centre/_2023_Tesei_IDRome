[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDRLab.ipynb)
[![DOI:10.1101/2023.05.08.539815](http://img.shields.io/badge/DOI-10.1101/2023.05.08.539815-B31B1B.svg)](https://doi.org/10.1101/2023.05.08.539815)
[![Video](http://img.shields.io/badge/►-Video-FF0000.svg)](https://youtu.be/v7YqJVEswM0)

# Analyses of conformational ensembles of the human IDRome

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks, and data for reproducing the results presented in the manuscript _Conformational ensembles of the human intrinsically disordered proteome: Bridging chain compaction with function and sequence conservation_.

The CSV file `IDRome_DB.csv` and the Excel Sheet `IDRome_DB.xlsx` list the sequence and various sequence- and conformational properties of all the 29,998 IDRs.

Simulation trajectories and time series of conformational properties are available for all the 29,998 IDRs at [sid.erda.dk/sharelink/AVZAJvJnCO](https://sid.erda.dk/sharelink/AVZAJvJnCO).

We also provide a [Notebook](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDRLab.ipynb) on [Google Colab](https://colab.research.google.com/) to generate conformational ensembles of user-supplied sequences using the [CALVADOS](https://github.com/KULL-Centre/CALVADOS) model. 


[![Video](https://img.youtube.com/vi/v7YqJVEswM0/default.jpg)](https://youtu.be/v7YqJVEswM0)

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

#### Commands to install [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/) v2.9.3 with [mphowardlab/azplugins](https://github.com/mphowardlab/azplugins) v0.11.0

```bash
    curl -LO https://github.com/glotzerlab/hoomd-blue/releases/download/v2.9.3/hoomd-v2.9.3.tar.gz
    tar xvfz hoomd-v2.9.3.tar.gz
    git clone https://github.com/mphowardlab/azplugins.git
    cd azplugins
    git checkout tags/v0.11.0
    cd ..
    cd hoomd-v2.9.3
    mkdir build
    cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=<path to python> \
        -DENABLE_CUDA=ON -DENABLE_MPI=ON -DSINGLE_PRECISION=ON -DENABLE_TBB=OFF \
        -DCMAKE_CXX_COMPILER=<path to g++> -DCMAKE_C_COMPILER=<path to gcc>
    make -j4
    cd ../hoomd
    ln -s ../../azplugins/azplugins azplugins
    cd ../build && make install -j4
```
