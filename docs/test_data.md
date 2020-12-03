# ramdaq: Usage: Using test data

Below commands automatically download annotation data and test FASTQ data and run the pipeline.

Currently, 8 combinations of annotations/data types are availabe:

- Species: mouse and human
- Sequencing: single-end and paired-end
- Strandness: unstranded and stranded (fr-firststrand)

## mouse, single-end

### mouse, single-end, unstranded

```bash
nextflow run ramdaq/main.nf -profile test_SE_UST_M,docker --outdir results_DLtest_SE_UST_M
```

### mouse, single-end, stranded

```bash
nextflow run ramdaq/main.nf -profile test_SE_ST_M,docker --outdir results_DLtest_SE_ST_M
```

## mouse, paired-end

### mouse, paired-end, unstranded

```bash
nextflow run ramdaq/main.nf -profile test_PE_UST_M,docker --outdir results_DLtest_PE_UST_M
```

### mouse, paired-end, stranded

```bash
nextflow run ramdaq/main.nf -profile test_PE_ST_M,docker --outdir results_DLtest_PE_ST_M
```

## human, single-end

### human, single-end, unstranded

```bash
nextflow run ramdaq/main.nf -profile test_SE_UST_H,docker --outdir results_DLtest_SE_UST_H
```

### human, single-end, stranded

```bash
nextflow run ramdaq/main.nf -profile test_SE_ST_H,docker --outdir results_DLtest_SE_ST_H
```

## human, paired-end

### human, paired-end, unstranded

```bash
nextflow run ramdaq/main.nf -profile test_PE_UST_H,docker --outdir results_DLtest_PE_UST_H
```

### human, paired-end, stranded

```bash
nextflow run ramdaq/main.nf -profile test_PE_ST_H,docker --outdir results_DLtest_PE_ST_H
```
