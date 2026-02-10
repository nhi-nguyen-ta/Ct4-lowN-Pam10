# Ct4-lowN-Pam10

This repository contains the programs to reproduce the analyses described in [manuscript-link].  

> This software was tested on [Linux](https://en.wikipedia.org/wiki/Linux), but it should also work on [Apple MacOS](https://en.wikipedia.org/wiki/MacOS) and Microsoft Windows, using [Windows Subsystem for Linux (WSL)](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

Given in input a set of bacterial proteomes, the programs will perform the following:
- Perform KEGG orthology annotation using Kofamscan  
- Generate heatmaps with annotated KEGG pathways

### Pre-required software
- [HMMER](http://hmmer.org) (v3.4)
- [Kofamascan](https://www.genome.jp/ftp/tools/kofam_scan) (v1.3.0)
- [Kofamscan profiles](https://www.genome.jp/ftp/db/kofam/) (we used the version from 2025-06-01)
- [Kofamscan ko_list](https://www.genome.jp/ftp/db/kofam/) (we used the version from 2025-06-01)

You can find a copy of the version of Kofamscan that we used in `.packages/kofam_scan-1.3.0.tar.gz`

NOTE: using a different release of the profiles DB should not give very different results, although we strongly suggest that the same release is used for maximum reproducibility.  

> IMPORTANT: although the license of the code in this repository is GPL-3, for the required third-party software the license of each software applies.  

### Python environment setup  
Before being able to execute the python program to perform the analysis a python environment needs to be created and required python packages installed.  
There are different ways to setup a virtual environment:
- Using Python pip and venv
- Using Conda/Anaconda
- Using [Astral uv](https://docs.astral.sh/uv/)  

Due to its speed and flexibility, we use uv from Astral.  

The required python version, and packages are detailed in the file `pyproject.toml`, in the main directory of this repository.

Please install `uv` following the following the instruction for your operating system (OS)

Follow these steps to setup the environment:
- `Clone the git repository`
  - git clone https://github.com/nhi-nguyen-ta/Ct4-lowN-Pam10.git
- `Enter the project directory`
  - cd Ct4-lowN-Pam10
- `Create the environment using uv`
  - uv sync
- `Activate the python environment (we assume you are using a bash shell)`
  - source .venv/bin/activate

If everything went well, your terminal should show `(ct4-lown-pam10)` indicating that the environment is active.


# Main commands to execute

#### Proteome annotation
 python ~/work-repos/genome_data_analysis/src/hiruma-proteomes.py -i ~/Desktop/hiruma-annotation-final-test/input/ -o ~/Desktop/hiruma-annotation-final-test -k ~/Desktop/test-kofamscan/kofam_scan-1.3.0/exec_annotation -p ~/Desktop/test-kofamscan/kofamscan/kofamdb/profiles/prokaryote.hal -ko ~/Desktop/test-kofamscan/kofamscan/kofamdb/ko_list --threads 10 --write-logs --file-extension fasta

#### generate heatmaps
python ~/work-repos/genome_data_analysis/src/hiruma-plot-kegg-annotation-heatmaps.py --brite-list ~/Desktop/hiruma-annotation-final-test/kofam-annotation/brite_counts/master.brite.list.tsv --output-dir ~/Desktop/hiruma-annotation-final-test/figures --strain-names-mapping ~/Desktop/hiruma-annotation-final-test/RAST_based_strain_names_with_extensions.tsv --threads 10 --write-logs --file-extension svg

