# MutaNET

## Introduction
Mutations in genomic key elements can influence gene expression and function in various ways, and hence greatly contribute to the phenotype. MutaNET comes with a next generation sequencing (NGS) pipeline that calls mutations based on paired-end NGS reads, an automated analysis tool and various file converters and mergers. The mutation analysis feature considers the coding region, protein domains, regulation and transcription factor binding site information, and can be used to analyse the potential impact of mutations on genes of interest.

MutaNET was developed and implemented in 2017 and published in 2018:
>Hollander, M, Hamed, M, Helms, V, Neininger, K (2018). MutaNET: a tool for automated analysis of genomic mutations in gene regulatory networks. Bioinformatics, 34, 5:864-866. https://www.ncbi.nlm.nih.gov/pubmed/29087464

Descriptions, explanations and examples of functionality and data sources can be found in `user_manual.pdf`. Detailed installation instructions for Windows, Linux and Mac OS X can be found in `installation_guide.pdf`. 

## Quickstart
### Starting MutaNET
MutaNET comes as Python 3 source code as well as an executable for Windows.

**Executable (Windows):** To start MutaNET, double–click on the executable `MutaNET32.exe` or `MutaNET64.exe`, depending on whether you have a 32–bit or 64–bit Windows installation. If you are not sure, choose the 32–bit executable. Make sure that the executable remains in the same directory as the config.yaml file. Otherwise the user interface will not start.

**From Source:** Open a command prompt or terminal and execute the following command:
```
python3 source_folder_path/mutaNET.py
```
or on Windows depending on your Python installation:
```
python source_folder_path/mutaNET.py
```
`Source_folder_path` is the path to the folder containing the source code of MutaNET. This requires Python 3 to be installed.

### Example Data Sets
When starting MutaNET for the first time, the file paths for small example data sets for the NGSpipeline, mutation analysis and file converters are already loaded to allow quick testing.  Keep in mind that for the NGS pipeline extra programs need to be installed. The  settings  for  the  example  data  can  be  restored  any  time  by  clicking  on `Settings → Restore default settings`.  The data sets can be found in the `example_data` folder, if you want to have a look at the file formats.

