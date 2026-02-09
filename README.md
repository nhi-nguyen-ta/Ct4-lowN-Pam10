# Ct4-lowN-Pam10

This repository contains the programs to reporduce the analyses described in <manuscript-link>.  

> This software was tested on [Linux](https://en.wikipedia.org/wiki/Linux) but it should work also on [Apple MacOS](https://en.wikipedia.org/wiki/MacOS) and Microsoft Windows using [Windows Subsystem for Linux (WSL)](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

Given in input a set of bacterial proteomes, the programs will perform the following:
- Perfrom KEGG orthology annotation using Kofamscan  
- Generate heatmaps with annotated KEGG pathways

### Pre-required software
- [HMMER](http://hmmer.org) (v3.4)
- [Kofamascan](https://www.genome.jp/ftp/tools/kofam_scan) (v1.3.0)
- [Kofamscan profiles](https://www.genome.jp/ftp/db/kofam/) (We used the version from 2025-06-01)
- [Kofamscan ko_list](https://www.genome.jp/ftp/db/kofam/) (We used the version from 2025-06-01)

You can find a copy of the version of Kofamscan that we used in __.packages/kofam_scan-1.3.0.tar.gz__

NOTE: using a different release of the profiles DB should not give very different results, although we strongly suggest use the same release we use for maximum reproducibility.  

> IMPORTANT: although the license of the code in this repository is GPL-3, for the required software the license of each software applies.  

### Python environment setup  
Before being able to execute the python program to perform the analysis a python environment need to be prepared.  
There are different ways to setup a virtual environment:
- Using Python pip and venv
- Using Conda/Anaconda
- Using [Astral uv](https://docs.astral.sh/uv/)  
Due to its speed and flexibility, we use uv from Astral.  