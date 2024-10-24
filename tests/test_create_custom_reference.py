import os, argparse
import pandas as pd
from unittest.mock import patch
import pytest
import synapseclient
import createCustomReference

REFERENCE_FA = os.path.join("tests", "test_input", "test_reference_genome.fa")
REFERENCE_GTF = os.path.join("tests", "test_input", "test_reference_genome.gtf")
TMP_FOLDER = os.path.join("tests", "tmp")


def mock_login(*args, **kwargs):
    return None

def mock_get(*args, **kwargs):
    if args[1] == "syn1111":
       return {"path": REFERENCE_FA}
    else:
        return {"path": REFERENCE_GTF}


@patch.object(synapseclient.Synapse, "login", mock_login)
@patch.object(synapseclient.Synapse, "get", mock_get)
def test_create_custom_reference():
    os.makedirs(TMP_FOLDER, exist_ok = True)
    
    args = argparse.ArgumentParser()
    args.added_genes_config_file = os.path.join("tests", "test_input", "test_config_file.csv")
    args.base_fasta_synid = "syn1111"
    args.base_gtf_synid = "syn1112"
    args.output_genome_name = os.path.join(TMP_FOLDER, "test_genome")
    args.authToken = "token123"
    args.scramble_genome = False

    createCustomReference.create_custom_reference(args)
        
    with open(args.output_genome_name + ".fa", "r") as fa_file:
        output_genome_fa = fa_file.read()
    
    with open(os.path.join("tests", "test_output", "test_reference_output.fa"), "r") as fa_file:
        true_genome_fa = fa_file.read()
        
    with open(args.output_genome_name + ".gtf", "r") as gtf_file:
        output_genome_gtf = gtf_file.read()
        
    with open(os.path.join("tests", "test_output", "test_reference_output.gtf"), "r") as gtf_file:
        true_genome_gtf = gtf_file.read()

    assert output_genome_fa == true_genome_fa
    assert output_genome_gtf == true_genome_gtf