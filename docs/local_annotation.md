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
