# ramdaq: Usage: Using provided annotations

## Download provided annotation

```bash
nextflow run rikenbit/ramdaq \
  --dl_references \
  -profile <docker/singularity> \
  --annot_ver <ANNOT_VER> \
  --outdir ramdaq_dl_references \
  --reportdir pipeline_info_$(date +%Y%m%d_%H%M)
```
The parameters are:

- `--dl_references`: Enables the one-time download of reference assets (genomes, indexes, annotations). Use this before running analyses so later runs can reuse the downloaded assets.
- `-profile`: Selects the runtime profile.
- `--annot_ver`: Annotation version key to download (human or mouse). Run `nextflow run rikenbit/ramdaq --help` to list available keys.
  - Approximate free disk space required per annotation version:
    - Mouse (GRCm39_vM37): ~7.0 GB
    - Human (GRCh38_v48): ~9.5 GB
    
```bash
  --annot_ver [str] Name of human or mouse reference to download : GRCm39_vM37, GRCm39_vM35, GRCm39_vM33, GRCm39_vM32, GRCm39_vM30_chr19, GRCm39_vM30, GRCm39_vM29, GRCm39_vM28, GRCm39_vM27, GRCm39_vM26, GRCm38_vM25, NCBIM37_vM1, GRCh38_v48, GRCh38_v46, GRCh38_v44, GRCh38_v43, GRCh38_v41_chr19, GRCh38_v41, GRCh38_v40, GRCh38_v39, GRCh38_v38, GRCh38_v37, GRCh38_v35, GRCh37_v19
```

- `--outdir results_test`: Output directory for the downloaded references. Choose a reusable location. Species subdirectories `human/` and `mouse/` are created automatically.
- `--reportdir` : Directory for pipeline reports (timeline/report/trace/dag). Use a unique name to avoid collisions when downloading multiple versions (e.g., pipeline_info_$(date +%Y%m%d_%H%M)). If a same-named report directory already exists, Nextflow may refuse to run to prevent overwriting.

If a download is interrupted, rerun with `-resume` to continue from where it left off.


## Run pipeline (mouse)

The below command processes a single-end and stranded FASTQ file using the GRCm39_vM26 mouse reference annotation:

```bash
nextflow run rikenbit/ramdaq \
  -profile docker \
  --genome GRCm39_vM37 \
  --local_annot_dir ramdaq_dl_references/mouse \
  --single_end \
  --stranded fr-firststrand \
  --outdir results_test \
  --reads 'https://bioinformatics.riken.jp/ramdaq/ramdaq_test_data/mouse/stranded_SE/SRR7993829_1.100K.fastq.gz'
```

The parameters are:

- `-profile docker`: Indicates the use of Docker.
- `--genome`: Specifies reference genome version.
  - In the current ramdaq, human and mouse references (the source is [gencode](https://www.gencodegenes.org/)) can be used for local annotation.
  - The latest versions and specified options are as follows.
    - Human
      - GRCh38
        - gencode.v48.primary_assembly.annotation: `--genome GRCh38_v48`
    - Mouse
      - GRCm39
        - gencode.vM37.primary_assembly.annotation: `--genome GRCm39_vM37`
- `--local_annot_dir`: Specifes the location of local path to a directiory containing annotation files.
- `--outdir results_test`: Specifies the output directory.
- `--single_end`: Indicates that the data is single-end. For paired-end data, just remove this parameter.
- `--stranded fr-firststrand`: Indicates that the data is stranded where the first read corresponds to the reverse complemented counterpart of a transcript. For unstranded data, just remove this parameter.
- `--reads`: Specifies the path of FASTQ files.
  - Single-end data example: `--reads fastq_files/*.R1.fastq.gz`
  - Paired-end data example: `--reads fastq_files/*.{R1,R2}.fastq.gz`

## Run pipeline (human)

The below command processes a paired-end and unstranded FASTQ file using the GRCh38_v37 human reference annotation:

```bash
wget https://bioinformatics.riken.jp/ramdaq/ramdaq_test_data/human/unstranded_PE/SRR12594145_1.100K.fastq.gz
wget https://bioinformatics.riken.jp/ramdaq/ramdaq_test_data/human/unstranded_PE/SRR12594145_2.100K.fastq.gz

nextflow run rikenbit/ramdaq \
  -profile docker \
  --genome GRCh38_v37 \
  --local_annot_dir ramdaq_annotation/human \
  --outdir results_test \
  --reads 'SRR12594145_{1,2}.100K.fastq.gz'
```

The parameters are:

- `-profile docker`: Indicates the use of Docker.
- `--genome`: Specifies reference genome version.
- `--local_annot_dir`: Specifes the location of local path to a directiory containing annotation files.
- `--outdir`: Specifies the output directory.
- `--reads`: Specifies the path of FASTQ files.
  - Single-end data example: `--reads fastq_files/*.R1.fastq.gz`
  - Paired-end data example: `--reads fastq_files/*.{R1,R2}.fastq.gz`

## Specification of mitocondorial gene, rRNA gene, histone gene GTF

- GTF annotation files for mitocondorial genes (`*.annotation.mt.gtf`), rRNA genes (`*RNA45S_merge.gtf`), and histone genes (`*.annotation.histone.gtf`) are based on the GENCODE annotation of the primary assembly available at <https://www.gencodegenes.org/>.

- Mitochondrial GTFs were created by extracting the lines matching "chrM" from the primary assembly annotations.

- Histone GTFs were created by listing the gene names of histones from the two sources and extracting the corresponding lines from primary assembly annotation.
  - (1) Databases of human/mouse genes:
    - Human histone genes list: <https://www.genenames.org/data/genegroup/#!/group/864>
    - Mouse histone genes list: <http://www.informatics.jax.org/mgihome/nomen/gene_name_initiative.shtml>
  - (2) Rule-based serach: "H[0-9]" string search for gene names in GENCODE annotation found histone genes that were not in the database above, so the following genes were added to GTFs.
    - "H1f0" "H4c11" "H4c12" "H2az2" "H1-10-AS1" "H2AZ1-DT" "H2AZP2" "H2AZP7" "H2AZP5" "H2AZP4" "H2AZP3" "H2AZ2P1" "H2AZP1" "H2AZP6" "H2ACP2" "H2ACP1" (human)
    - "H1f0" "H4c14" "H4c11" "H4c12" "H2az2" (mouse)

- The rRNA GTF for each species was created by generating two GTF files on the UCSC Table Browser (<http://genome.ucsc.edu/cgi-bin/hgTables>) and mering them.
  - Procedures to create the mouse rRNA annotation:
    - (1) Set mouse genome and Dec.2011(GRCm38/mm10) assembly in UCSC Table Browser.
    - (2) First GTF: Select Variation and Repeats group, RepeatMasker track, and specify "rRNA" in repClass column of Filter. Output "RepeatMasker rRNA" GTF.
    - (3) Second GTF: Select Genes and Gene Predictions group, NCBI RefSeq track, UCSC RefSeq(refGene) table, and set "NR_046233" in the name column of the filter. Output "Rn45s 45S pre-ribosomal RNA" GTF.
    - (4) Merge the two GTFs.
  - Procedures to create the human rRNA annotation:
    - (1) Set human genome and Dec.2013(GRCh38/hg38) assembly in UCSC Table Browser.
    - (2) First GTF: Select Repeats group, RepeatMasker track, and specify "rRNA" in repClass column of Filter. Output "RepeatMasker rRNA" GTF.
    - (3) Second GTF: select Genes and Gene Predictions group, NCBI RefSeq track, UCSC RefSeq(refGene) table, and set "NR_046235" in the name column of the filter. Output "45S pre-ribosomal RNA" GTF.
    - (4) Merge the two GTFs.
  - Since the chromosome names aere different between UCSC and GENCODE, we renamed the chromosome names of "scaffolds" chromosomes annotated by UCSC to "GenBank.Accn" format in the corresponding table below:
    - "GCA_000001405.28_GRCh38.p13_assembly_report.txt" at [ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13)
      - Please refer to GenBank-Accn column and UCSC-style-name column.
