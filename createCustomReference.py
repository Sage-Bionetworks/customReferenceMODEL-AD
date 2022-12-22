import subprocess
import synapseclient, pandas, os
from synapseclient import File
import argparse

# usage: python3 createCustomReference.py --username *** --password ***
# probably needs testing outside of experimental testing

parser = argparse.ArgumentParser(description='Create custom reference genome for Model AD')

parser.add_argument('--username', help='Enter Synapse username')
parser.add_argument('--password', help='enter Synapse password')

args = parser.parse_args()

# probably not hard code
parentDirectory = 'syn50115986'

# download data from Synapses

syn = synapseclient.Synapse()
syn.login(args.username, args.password)

children = syn.getChildren(parent = parentDirectory)
for x in children:
	syn.get(x['id'], downloadLocation = "references")

# https://github.com/gencorefacility/reform

# python3 reform.py 
#   --chrom=<chrom> \
#   --position=<pos> \ 
#   --in_fasta=<in_fasta> \
#   --in_gff=<in_gff> \
#   --ref_fasta=<ref_fasta> \
#   --ref_gff=<ref_gff>

# Chromosome/Mouse Genome Reference Information
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27
# GRCm39
# Genome Reference Consortium Mouse Build 39

# APOE: 19
# APP: 11
# PS1: 14

addedGenes = ['APOE_Human.fa', 'APP_human.fa', 'PS1_human.fa']
addedGTF = ['APOE_Human.mod.gtf', 'APP_Human.mod.gtf', 'PS1_human.gtf']
onChromosome = ['NC_000085.7', 'NC_000077.7', 'NC_000080.7']

# new reference genome name
newReferenceName = "universal_MODEL_AD_Reference.fa"
# new gtf file name
newGTFName = "universal_MODEL_AD_Reference.gtf"

def addGenetoReference(geneFasta, geneGTF, chromosome, ref_fasta, ref_GTF):

    pre, ext = os.path.splitext(geneFasta)
    print("Running " + pre + "...")

    # affix to end of chromosome without replacement
    chrom = chromosome
    position = '-1'
    in_fasta = geneFasta
    in_gff = "Human_Genes_GTF/" + geneGTF
    ref_fasta = ref_fasta
    ref_gff = ref_GTF

    arguments = ['python3', './reform/reform.py', '--chrom', chrom, '--position', position, '--in_fasta', in_fasta, '--in_gff', in_gff, '--ref_fasta', ref_fasta, '--ref_gff', ref_gff]
    p = subprocess.run(arguments, capture_output=True)

    # print out error messages
    print( 'exit status:', p.returncode )
    print( 'stdout:', p.stdout.decode() )
    print( 'stderr:', p.stderr.decode() )

    # rename output files
    oldReferenceName = os.path.basename(ref_fasta)
    pre, ext = os.path.splitext(oldReferenceName)
    os.rename(pre + "_reformed.fa", newReferenceName)

    oldGTFName = os.path.basename(ref_GTF)
    pre, ext = os.path.splitext(oldReferenceName)
    os.rename(pre + "_reformed.gtf", newGTFName)

originalReferenceLocation = './references/GCF_000001635.27_GRCm39_genomic.fna'
originalGTFLocation = './references/GCF_000001635.27_GRCm39_genomic.gtf'

# loop through each gene

for i in range(0, 3):
    if(i == 0):
        addGenetoReference(geneFasta=addedGenes[i], geneGTF=addedGTF[i], chromosome=onChromosome[i], ref_fasta=originalReferenceLocation, ref_GTF = originalGTFLocation)
    else:
        addGenetoReference(geneFasta=addedGenes[i], geneGTF=addedGTF[i], chromosome=onChromosome[i], ref_fasta=newReferenceName, ref_GTF = newGTFName)

# upload to Synapse

file = File(path=newReferenceName, parent=parentDirectory)
file = syn.store(file)

file = File(path=newGTFName, parent=parentDirectory)
file = syn.store(file)
