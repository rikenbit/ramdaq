# nf-core/ramdaq: Usage: Using provided annotations

## Download provided annotation (mouse)

```bash
wget -r -np -nH --cut-dirs=1 -l 1 -R "index.html*" https://bioinformatics.riken.jp/ramdaq_nfcore/ramdaq_annotation/mouse/
tar -xvf ramdaq_annotation/mouse/*.tar.gz -C ramdaq_annotation/mouse
```

## Download provided annotation (human)

```bash
wget -r -np -nH --cut-dirs=1 -l 1 -R "index.html*" https://bioinformatics.riken.jp/ramdaq_nfcore/ramdaq_annotation/human/
tar -xvf ramdaq_annotation/human/*.tar.gz -C ramdaq_annotation/human
```

## Run pipeline (mouse)

The below command processes a single-end and stranded FASTQ file using the GRCm38 mouse reference annotation:

```bash
nextflow run ramdaq.nf/main.nf -profile docker --genome GRCm38 --local_annot_dir ramdaq_annotation/mouse --single_end --stranded --outdir results_test --reads 'https://bioinformatics.riken.jp/ramdaq_nfcore/ramdaq_test_data/mouse/stranded_SE/SRR7993829_1.100K.fastq.gz'
```

The parameters are:

- `-profile docker`: Indicates the use of Docker.
- `--genome GRCm38`: Specifies reference genome version.
- `--local_annot_dir`: Specifes the location of local path to a directiory containing annotation files.
- `--outdir results_test`: Specifies the output directory.
- `--single_end`: Indicates that the data is single-end.
  - For paired-end data, just remove this parameter.
- `--stranded`: Indicates that the data is stranded.
  - For unstranded data, just remove this parameter.
- `--reads`: Specifies the path of FASTQ files.
  - Single-end data example: `--reads fastq_files/*.R1.fastq.gz`
  - Paired-end data example: `--reads fastq_files/*.{R1,R2}.fastq.gz`


## Run pipeline (human)

The below command processes a paired-end and unstranded FASTQ file using the GRCh38 human reference annotation:

```bash
wget https://bioinformatics.riken.jp/ramdaq_nfcore/ramdaq_test_data/human/unstranded_PE/SRR12594145_1.100K.fastq.gz
wget https://bioinformatics.riken.jp/ramdaq_nfcore/ramdaq_test_data/human/unstranded_PE/SRR12594145_2.100K.fastq.gz

nextflow run ramdaq.nf/main.nf -profile docker --genome GRCh38 --local_annot_dir ramdaq_annotation/human --outdir results_test --reads 'SRR12594145_{1,2}.100K.fastq.gz'
```

The parameters are:

- `-profile docker`: Indicates the use of Docker.
- `--genome GRCh38`: Specifies reference genome version.
- `--local_annot_dir`: Specifes the location of local path to a directiory containing annotation files.
- `--outdir results_test`: Specifies the output directory.
- `--reads`: Specifies the path of FASTQ files.
  - Single-end data example: `--reads fastq_files/*.R1.fastq.gz`
  - Paired-end data example: `--reads fastq_files/*.{R1,R2}.fastq.gz`

## Specification of mitocondorial gene, rrna, histone gene GTF

- All three annotations are based on the primary assembly annotation available at https://www.gencodegenes.org/.

- Mitochondrial GTFs were created by extracting the lines matching "chromosome M" from the primary assembly annotations.

- Histone GTFs were created by listing the gene names of histones concerning the following database and extracting the corresponding lines from primary assembly annotation.
  - human histone genes list : https://www.genenames.org/data/genegroup/#!/group/864
  - mouse histone genes list : http://www.informatics.jax.org/mgihome/nomen/gene_name_initiative.shtml
- Also, "H[0-9]" string search for gene names in GENCODE annotation found a histone that was not in the database above, so the following was added to GTFs.
  - "H1f0" "H4c11" "H4c12" "H2az2" "H1-10-AS1" "H2AZ1-DT" "H2AZP2" "H2AZP7" "H2AZP5" "H2AZP4" "H2AZP3" "H2AZ2P1" "H2AZP1" "H2AZP6" "H2ACP2" "H2ACP1" (human)
  - "H1f0" "H4c14" "H4c11" "H4c12" "H2az2" (mouse)

- The rRNA GTF is based on annotations obtained from the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables).
  - mouse rRNA annotation
    - Set mouse genome and Dec.2011(GRCm38/mm10) assembly in UCSC Table Browser.
    - Select Variation and Repeats group, RepeatMasker track, and specify "rRNA" in repClass column of Filter.
    - Output "RepeatMasker rRNA" GTF.
    - Next, select Genes and Gene Predictions group, NCBI RefSeq track, UCSC RefSeq(refGene) table, and set "NR_046233" in the name column of the filter.
    - Output "Rn45s 45S pre-ribosomal RNA" GTF.
    - Merge two GTFs.
  - human rRNA annotation
    - Set human genome and Dec.2013(GRCh38/hg38) assembly in UCSC Table Browser.
    - Select Repeats group, RepeatMasker track, and specify "rRNA" in repClass column of Filter.
    - Output "RepeatMasker rRNA" GTF.
    - Next, select Genes and Gene Predictions group, NCBI RefSeq track, UCSC RefSeq(refGene) table, and set "NR_046235" in the name column of the filter.
    - Output "45S pre-ribosomal RNA" GTF.
    - Merge two GTFs.
  - Since the format of chromosome names is different between UCSC and Gencode, rename the chromosome names of "scaffolds" chromosomes annotated by UCSC to "GenBank.Accn" format in the corresponding table below.
    - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13
      - Please refer to GenBank-Accn column and UCSC-style-name column.