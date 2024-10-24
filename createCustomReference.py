import subprocess, os, argparse, random
import csv, textwrap
import synapseclient
import pandas as pd

""""
# usage: python3 createCustomReference.py [--authToken ***]

Human gene sequences were downloaded from Ensembl (GRCh38.112) as fasta files using the "Export Data" feature on each 
gene page:
    APOE - https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000130203;r=19:44905791-44909393
    APP - https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000142192;r=21:25880535-26171128
    MAPT - https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000186868;r=17:45894527-46028334
    PSEN1 - https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000080815;r=14:73136418-73223691
The output was "FASTA Sequence", the strand was "Forward strand", and we de-selected all options under "Options for 
FASTA sequence". All other fields were left as default.

Annotations for human genes were downloaded as a GTF file from Ensembl (GRCh38.112) from
https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz and the annotations for APOE, 
APP, MAPT, and PSEN1 were extracted into separate files using the command line:
    cat Homo_sapiens.GRCh38.112.gtf | grep ENSG00000130203 > APOE_human.gtf
replacing the Ensembl ID and filename for each gene.

The base mouse genome fasta and GTF files were downloaded from Ensembl (GRCm39.112) from
https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
and https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz, and put on Synapse at 
syn61779422 and syn61779397 for re-use.
"""

# Constants - folders within this project folder where temporary files will be stored
MM_BASE_REFERENCE_FOLDER = "references"
EDITED_HUMAN_FILES_LOCATION = "edited_human_genes"

# Constants - defaults for command-line arguments
MM_BASE_FASTA_SYNID = "syn61779422"
MM_BASE_GTF_SYNID = "syn61779397"
ADDED_GENES_CONFIG_FILE = "added_genes_config.csv"
DEFAULT_OUTPUT_GENOME_NAME = "universal_MODEL_AD_reference"

# TODO ensure human gene fasta/gtf files end in \n before concatenation


def create_custom_reference(args):
    os.makedirs(MM_BASE_REFERENCE_FOLDER, exist_ok=True)
    os.makedirs(EDITED_HUMAN_FILES_LOCATION, exist_ok=True)

    # Download mouse genome reference from Synapse
    syn = synapseclient.Synapse()
    syn.login(authToken=args.authToken)

    print("Downloading reference FASTA and GTF files...")
    ref_fasta_gz = syn.get(
        args.base_fasta_synid, downloadLocation=MM_BASE_REFERENCE_FOLDER
    )
    if ref_fasta_gz["path"].endswith(".gz"):
        subprocess.run("gunzip --keep --force " + ref_fasta_gz["path"], shell=True)
        original_reference_location = ref_fasta_gz["path"][:-3]  # Chop off ".gz"
    else:
        original_reference_location = ref_fasta_gz["path"]

    ref_gtf_gz = syn.get(args.base_gtf_synid, downloadLocation=MM_BASE_REFERENCE_FOLDER)
    if ref_gtf_gz["path"].endswith(".gz"):
        subprocess.run("gunzip --keep --force " + ref_gtf_gz["path"], shell=True)
        original_gtf_location = ref_gtf_gz["path"][:-3]  # Chop off ".gz"
    else:
        original_gtf_location = ref_gtf_gz["path"]

    added_genes_df = pd.read_csv(args.added_genes_config_file)

    print("Adding the following genes to the reference genome:")
    print(added_genes_df)

    new_reference_name = args.output_genome_name + ".fa"
    new_gtf_name = args.output_genome_name + ".gtf"

    arguments_fa = ["cat", original_reference_location]
    arguments_gtf = ["cat", original_gtf_location]

    for i in range(0, added_genes_df.shape[0]):
        (edited_fa_filename, edited_gtf_filename) = modify_human_sequence(
            fa_filename=added_genes_df["fasta_file"][i],
            gtf_filename=added_genes_df["gtf_file"][i],
            chromosome=added_genes_df["chromosome"][i],
            edited_folder=EDITED_HUMAN_FILES_LOCATION,
            scramble_genome=args.scramble_genome,
        )

        arguments_fa = arguments_fa + [edited_fa_filename]
        arguments_gtf = arguments_gtf + [edited_gtf_filename]

    arguments_fa = arguments_fa + [">", new_reference_name]
    arguments_gtf = arguments_gtf + [">", new_gtf_name]

    subprocess.run(" ".join(arguments_fa), shell=True)
    subprocess.run(" ".join(arguments_gtf), shell=True)

    print("Creating Ensembl ID to gene name map...")
    df = pd.read_table(
        new_gtf_name, sep="\t", header=None, comment="#", dtype={0: "object"}
    )
    df = df[df[2] == "gene"]
    fields = df[8].str.split("; ")
    gene_ids = fields.apply(
        lambda row: [
            item.replace("gene_id ", "").replace('"', "")
            for item in row
            if "gene_id" in item
        ]
    )
    gene_names = fields.apply(
        lambda row: [
            item.replace("gene_name ", "").replace('"', "")
            for item in row
            if "gene_name" in item
        ]
    )

    id_map = pd.DataFrame({"ensembl_gene_id": gene_ids, "gene_symbol": gene_names})
    id_map["ensembl_gene_id"] = id_map["ensembl_gene_id"].apply(
        lambda row: row[0] if len(row) > 0 else ""
    )
    id_map["gene_symbol"] = id_map["gene_symbol"].apply(
        lambda row: row[0] if len(row) > 0 else ""
    )
    id_map = id_map.drop_duplicates()

    symbol_map_filename = args.output_genome_name + "_symbol_map.csv"
    id_map.to_csv(symbol_map_filename, index=False)

    print("gzipping fasta file...")
    subprocess.run("gzip --keep " + new_reference_name, shell=True)

    print("gzipping gtf file...")
    subprocess.run("gzip --keep " + new_gtf_name, shell=True)

    print(
        "Completed reference files:\n"
        + "\tReference FASTA: " + new_reference_name + ".gz\n" 
        + "\tReference GTF:   " + new_gtf_name + ".gz\n" 
        + "\tSymbol map:      " + symbol_map_filename
    )


def modify_human_sequence(
    fa_filename, gtf_filename, chromosome, edited_folder, scramble_genome
):
    gtf_df = pd.read_table(gtf_filename, sep="\t", header=None)
    gtf_df[0] = chromosome

    with open(fa_filename, "r") as old_fasta:
        fasta_lines = old_fasta.readlines()

        fasta_fields = fasta_lines[0].split(":")
        fasta_fields[0] = ">" + str(chromosome) + " dna"
        fasta_fields[3] = str(chromosome)

        orig_start_pos = int(fasta_fields[4])
        fasta_fields[4] = "1"
        fasta_fields[5] = str(int(fasta_fields[5]) - orig_start_pos + 1)

    gtf_df[3] = gtf_df[3] - orig_start_pos + 1
    gtf_df[4] = gtf_df[4] - orig_start_pos + 1

    # For benchmarking purposes only
    if scramble_genome:
        fasta_string = list("".join(fasta_lines[1:]).replace("\n", ""))

        random.seed(chromosome)
        random.shuffle(fasta_string)

        # Split back into lines of length 60 (as in original fasta file) and re-add \n characters to each line
        fasta_list = textwrap.wrap("".join(fasta_string), width=60)
        fasta_list = [line + "\n" for line in fasta_list]
    # Normal non-benchmarking case
    else:
        fasta_list = fasta_lines[1:]

    fasta_list = [":".join(fasta_fields)] + fasta_list

    edited_fa_filename = os.path.join(edited_folder, os.path.basename(fa_filename))
    edited_gtf_filename = os.path.join(edited_folder, os.path.basename(gtf_filename))

    with open(edited_fa_filename, "w") as new_fasta:
        new_fasta.writelines(fasta_list)

    gtf_df.to_csv(
        edited_gtf_filename, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE
    )

    return (edited_fa_filename, edited_gtf_filename)


def get_input_args():
    parser = argparse.ArgumentParser(
        description="Create a custom reference genome for Model AD that contains sequences for human transgenes."
    )

    parser.add_argument(
        "-a",
        "--added_genes_config_file",
        required=False,
        metavar="FILE",
        default=ADDED_GENES_CONFIG_FILE,
        help="Path (local or absolute) to the config CSV defining which genes to add and where they should be added.",
    )
    parser.add_argument(
        "-f",
        "--base_fasta_synid",
        required=False,
        metavar="FA_SYNID",
        default=MM_BASE_FASTA_SYNID,
        help="Synapse ID of the base mouse genome fasta to build on. The file at this location must be a .fa.gz file.",
    )
    parser.add_argument(
        "-g",
        "--base_gtf_synid",
        required=False,
        metavar="GTF_SYNID",
        default=MM_BASE_GTF_SYNID,
        help="Synapse ID of the base mouse genome GTF to build on. The file at this location must be a .gtf.gz file.",
    )
    parser.add_argument(
        "-o",
        "--output_genome_name",
        required=False,
        metavar="NAME",
        default=DEFAULT_OUTPUT_GENOME_NAME,
        help="Name of the output genome (without any file extension)",
    )
    parser.add_argument(
        "-s",
        "--scramble_genome",
        required=False,
        default=False,
        action="store_true",
        help="For benchmarking purposes only. If this flag is added, the human gene sequences will be randomly "
        + "shuffled before addition to the reference genome.",
    )
    parser.add_argument("-t", "--authToken", required=False, help="Synapse auth token")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = get_input_args()
    create_custom_reference(args)
