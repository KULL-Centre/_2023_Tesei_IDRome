[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDRLab.ipynb)
[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDR_SVR_predictor.ipynb)
[![DOI:10.1101/2023.05.08.539815](http://img.shields.io/badge/DOI-10.1101/2023.05.08.539815-B31B1B.svg)](https://doi.org/10.1101/2023.05.08.539815)
[![Video](http://img.shields.io/badge/►-Video-FF0000.svg)](https://youtu.be/v7YqJVEswM0)

# Conformational ensembles of the human IDRome

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks, and data for reproducing the results presented in the manuscript _Conformational ensembles of the human intrinsically disordered proteome_.

The CSV file `IDRome_DB.csv` lists amino acid sequences, sequence features, and conformational properties of all the 28,058 IDRs.

Simulation trajectories and time series of conformational properties are available for all the IDRs at [sid.erda.dk/sharelink/AVZAJvJnCO](https://sid.erda.dk/sharelink/AVZAJvJnCO).

We also provide Notebooks on [Google Colab](https://colab.research.google.com/) to (i) generate conformational ensembles of user-supplied sequences using the [CALVADOS](https://github.com/KULL-Centre/CALVADOS) model and (ii) predict scaling exponents and conformational entropies per residue using the SVR models:
- [`IDRLab.ipynb`](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDRLab.ipynb)
- [`IDR_SVR_predictor.ipynb`](https://colab.research.google.com/github/KULL-Centre/_2023_Tesei_IDRome/blob/main/IDR_SVR_predictor.ipynb)


[![Video](https://img.youtube.com/vi/v7YqJVEswM0/default.jpg)](https://youtu.be/v7YqJVEswM0)

### Layout
- `seq_conf_prop.ipynb` reproduces Fig. 1, 3, and Extended Data Fig. 2, 5, 6e-t, and 7
- `go_analysis.ipynb` reproduces Fig. 2
- `conservation_analysis.ipynb` reproduces Fig. 4
- `clinvar_fmug.ipynb` reproduces Fig. 5 and Extended Data Fig. 9
- `uniprot_domains.ipynb` reproduces Extended Data Fig. 1
- `svr_models.ipynb` reproduces Extended Data Fig. 8
- `go_uniprot_calls.ipynb` performs API calls to obtain gene ontology terms from UniProt
- `calc_seq_prop.ipynb` and `calc_seq_prop_SPOT.ipynb` compute sequence descriptors and generate the `IDRome_DB.csv` and `IDRome_DB_SPOT.csv` files
- `CALVADOS_tests.ipynb` reproduces Extended Data Fig. 3
- `AF2_PAEs.ipynb` reproduces Extended Data Fig. 4
- `CD-CODE.ipynb` reproduces Extended Data Fig. 6a-d
- `md_simulations/` contains code and data related to single-chain simulations performed using the CALVADOS model and [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/) v2.9.3 installed with [mphowardlab/azplugins](https://github.com/mphowardlab/azplugins)
- `idr_selection/` contains code and data to generate the pLDDT-based and SPOT-based sets of IDRs
- `idr_orthologs/` contains code and data to generate the set of orthologs of human IDRs
- `svr_models/` contains scikit-learn SVR models generated in `svr_models.ipynb`
- `zscores/` contains code and data to calculate [NARDINI](https://github.com/mshinn23/nardini) z-scores
- `go_analyses/` contains input and output data related to the Gene Ontology analyses in `go_analysis.ipynb`
- `QCDPred/` contains code and data related to [QCD calculations](https://github.com/KULL-Centre/papers/tree/main/2022/degron-predict-Johansson-et-al)
- `clinvar_fmug_cdcode/` contains code and data related to the analysis of the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [FMUG](https://fmug.amaral.northwestern.edu/), and [CD-CODE](https://cd-code.org/) databases

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
