# Create a custom genome reference for MODEL-AD

This repository provides code to concatenate human genes onto the mouse reference genome (GRCm39.112) to create a
custom genome reference that can be used for all MODEL-AD models. This allows for RNA sequencing alignment to human
genes for mice with transgenes. By default, the human genes APOE, APP, CLU, MAPT, and PSEN1 are added to the mouse
reference by this script.

# Running the script

The script can be run using all default values with Python:
```
python createCustomReference.py
```

This script also takes optional arguments to change various inputs for the script. Run the above command with `--help` to see the full list. Some key arguments:
* `--authToken <token>`: Input your Synapse auth token if your computer does not have your Synapse credentials saved.
* `-o <name>`: Rename your output genome to something other than the default ("universal_MODEL_AD_reference")
* `-a <file>`: Input a different gene configuration CSV than the default ("added_genes_config.csv"), which will allow you to add/remove which genes are added to the reference. Format of the new CSV must be identical to "added_genes_config.csv" in this repository.

# Reference sequence origins

## Mouse reference genome

The base mouse genome FASTA and GTF files were downloaded from Ensembl (GRCm39.112):
* FASTA: [Mus_musculus.GRCm39.dna.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)
* GTF: [Mus_musculus.GRCm39.112.gtf.gz](https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz)

and put on Synapse at 
[syn61779422](https://www.synapse.org/Synapse:syn61779422) and [syn61779397](https://www.synapse.org/Synapse:syn61779397) for re-use.


## Human reference sequences

Human gene sequences were downloaded from Ensembl (GRCh38.112) ([archive](https://may2024.archive.ensembl.org/)) as fasta files using the "Export Data" feature on each 
gene page:
* APOE - [ENSG00000130203](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000130203&db=core)
* APP - [ENSG00000142192](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000142192&db=core)
* CLU - [ENSG00000120885](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000120885&db=core)
* MAPT - [ENSG00000186868](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000186868&db=core)
* PSEN1 - [ENSG00000080815](https://may2024.archive.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000080815&db=core)

The output was "FASTA Sequence", the strand was "Forward strand", and we de-selected all options under "Options for 
FASTA sequence". All other fields were left as default.

Annotations for human genes were downloaded as a GTF file from Ensembl (GRCh38.112) ([Homo_sapiens.GRCh38.112.gtf.gz](https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz)) and the annotations for APOE, 
APP, CLU, MAPT, and PSEN1 were extracted into separate files using the command line:
```
cat Homo_sapiens.GRCh38.112.gtf | grep ENSG00000130203 > APOE_human.gtf
```
replacing the Ensembl ID and filename for each gene.
