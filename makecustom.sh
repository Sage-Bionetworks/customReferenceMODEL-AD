#Step1: Split human reference genome into chromosmes
mkdir human_chr_fastafiles
cd human_chr_fastafiles
csplit -s -z /path_to_human_reference_genome_file/Homo_sapiens.GRCh38.dna.primary_assembly.fa '/>/' '{*}'
for i in xx* ; do \
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  mv "$i" "$n.fa" ; \
 done

 #Step2:  create index for chromosome  with gene of interest. For e.g. human APOE gene is on chromosome 19
 samtools faidx 19.fa
 #Step3: extract human APOE gene sequences using location bases on version/build of human genome and store in a file
 samtools faidx 19.fa 19:44905754-44909393 >APOE_human.fa
 # rename header of extracted sequence
 sed 's/>.*/>21/' APOE_human.fa > APOE_Human.fa

 # Step4: Extract gene information of gene of interest from gtf file of human

 #Step5:change start and end f gene to inserts in mouse genome
 awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4-44905753"\t"$5-44905753"\t"$6"\t"$7"\t"$8"\t"$9}' APOE_Human.gtf | sed 's/^19/21/g' > APOE_Human.mod.gtf

#Step6a: Append APOE_Human.fa to mouse reference genome_file
cat APOE_Human.fa mouse_ref.fa > mouse_custom_ref.fa
#step6b: Append APOE_Human.mod.gtf to mouse reference gtf file
 cat APOE_Human.mod.gtf mouse_gtf.fa > mouse_custom_gtf.fa
