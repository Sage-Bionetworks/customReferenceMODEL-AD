import os, argparse, shutil
from unittest.mock import patch
from io import StringIO
import pytest
import synapseclient
import createCustomReference

TMP_FOLDER = os.path.join("tests", "tmp")
REFERENCE_FA = os.path.join("tests", "test_input", "test_reference_genome.fa")
REFERENCE_GTF = os.path.join("tests", "test_input", "test_reference_genome.gtf")
EDITED_FILE_PREFIX = os.path.join("tests", "test_input", "edited_")

# For modify_human_sequence testing
ADD_FA_INPUT = os.path.join("tests", "test_input", "test_gene.fa")
ADD_GTF_INPUT = os.path.join("tests", "test_input", "test_gene.gtf")
FA_OUTPUT = os.path.join("tests", "test_output", "test_gene_output.fa")
FA_OUTPUT_SCRAMBLED = os.path.join(
    "tests", "test_output", "test_gene_scrambled_output.fa"
)
GTF_OUTPUT = os.path.join("tests", "test_output", "test_gene_output.gtf")


class MockSynapse:
    def login(self, *args, **kwargs):
        return None

    def get(self, *args, **kwargs):
        # For this test we use file paths as Synapse IDs, so args[0] will be a path instead of an ID
        return synapseclient.File(path=args[0], parent="syn1111")


def mock_modify_human_sequence(*args, **kwargs):
    # Both "fa_synid" and "gtf_synid" are passed in as file paths instead of Synapse IDs in this test.
    # Assumes that the pre-made edited files exist in tests/test_input and are named "edited_<original_filename>"
    edited_fa_filename = EDITED_FILE_PREFIX + os.path.basename(kwargs["fa_synid"])
    edited_gtf_filename = EDITED_FILE_PREFIX + os.path.basename(kwargs["gtf_synid"])
    return (edited_fa_filename, edited_gtf_filename)


class TestCreateCustomReference:
    # Replaces the CSV config file
    genes_config_string = (
        "gene_symbol,fasta_synid,gtf_synid,chromosome\n"
        + "TEST1,tests/test_input/test_gene.fa,tests/test_input/test_gene.gtf,25\n"
        + "TEST2,tests/test_input/test_gene2.fa,tests/test_input/test_gene2.gtf,20"
    )

    @pytest.fixture(scope="function", autouse=True)
    def setup_method(self):
        self.args = argparse.ArgumentParser()
        self.args.added_genes_config_file = StringIO(self.genes_config_string)
        self.args.base_fasta_synid = REFERENCE_FA
        self.args.base_gtf_synid = REFERENCE_GTF
        self.args.output_genome_name = "test_genome"
        self.args.output_folder_path = TMP_FOLDER
        self.args.authToken = "token123"
        self.args.scramble_genome = False

        # Clear files from previous tests
        if os.path.exists(TMP_FOLDER):
            shutil.rmtree(TMP_FOLDER)

        os.makedirs(TMP_FOLDER, exist_ok=True)

    @patch.object(synapseclient, "Synapse", MockSynapse)
    @patch.object(
        createCustomReference, "modify_human_sequence", mock_modify_human_sequence
    )
    def test_create_custom_reference(self):
        createCustomReference.create_custom_reference(self.args)

        output_genome = os.path.join(
            self.args.output_folder_path, self.args.output_genome_name
        )
        expected_output = os.path.join("tests", "test_output", "test_reference_output")
        with open(output_genome + ".fa", "r") as fa_file:
            output_genome_fa = fa_file.read()

        with open(expected_output + ".fa", "r") as fa_file:
            true_genome_fa = fa_file.read()

        with open(output_genome + ".gtf", "r") as gtf_file:
            output_genome_gtf = gtf_file.read()

        with open(expected_output + ".gtf", "r") as gtf_file:
            true_genome_gtf = gtf_file.read()

        assert output_genome_fa == true_genome_fa
        assert output_genome_gtf == true_genome_gtf

    @pytest.mark.parametrize("scramble_genome", (True, False))
    def test_modify_human_sequence(self, scramble_genome):
        chromosome = 25
        syn = MockSynapse()
        new_fa_file, new_gtf_file = createCustomReference.modify_human_sequence(
            syn=syn,
            fa_synid=ADD_FA_INPUT,
            gtf_synid=ADD_GTF_INPUT,
            chromosome=chromosome,
            edited_folder=TMP_FOLDER,
            scramble_genome=scramble_genome,
        )

        with open(new_fa_file, "r") as fa_file:
            new_fa = fa_file.read()

        if scramble_genome:
            fa_output_file = FA_OUTPUT_SCRAMBLED
        else:
            fa_output_file = FA_OUTPUT

        with open(fa_output_file, "r") as fa_file:
            true_fa = fa_file.read()

        with open(new_gtf_file, "r") as gtf_file:
            output_gtf = gtf_file.read()

        with open(GTF_OUTPUT, "r") as gtf_file:
            true_gtf = gtf_file.read()

        assert new_fa == true_fa
        assert output_gtf == true_gtf
