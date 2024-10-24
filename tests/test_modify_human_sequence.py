import os
import pandas as pd
import pytest
import createCustomReference

# Test files
FA_INPUT = os.path.join("tests", "test_input", "test_gene.fa")
GTF_INPUT = os.path.join("tests", "test_input", "test_gene.gtf")
FA_OUTPUT = os.path.join("tests", "test_output", "test_gene_output.fa")
FA_OUTPUT_SCRAMBLED = os.path.join("tests", "test_output", "test_gene_scrambled_output.fa")
GTF_OUTPUT = os.path.join("tests", "test_output", "test_gene_output.gtf")

TMP_FOLDER = os.path.join("tests", "tmp")

@pytest.mark.parametrize("scramble_genome", (True, False))
def test_modify_human_sequence(scramble_genome):
    chromosome = 25
    os.makedirs(TMP_FOLDER, exist_ok = True)
    
    new_fa_file, new_gtf_file = createCustomReference.modify_human_sequence(
        fa_filename = FA_INPUT, 
        gtf_filename = GTF_INPUT, 
        chromosome = chromosome, 
        edited_folder = TMP_FOLDER,
        scramble_genome = scramble_genome
    )
    
    with open(new_fa_file, "r") as fa_file:
        new_fa = fa_file.read()
    
    if (scramble_genome):
        fa_output_file = FA_OUTPUT_SCRAMBLED
    else:
        fa_output_file = FA_OUTPUT
        
    with open(fa_output_file, "r") as fa_file:
        true_fa = fa_file.read()
        
    new_gtf = pd.read_table(new_gtf_file, sep="\t", header=None)
    true_gtf = pd.read_table(GTF_OUTPUT, sep="\t", header=None)
    
    assert new_fa == true_fa
    pd.testing.assert_frame_equal(new_gtf ,true_gtf)
    