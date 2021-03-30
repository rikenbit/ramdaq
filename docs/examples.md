# Examples

- [Examples](#examples)
  - [Case 1: Unstranded, Paired-end, Mouse](#case-1-unstranded-paired-end-mouse)
  - [Case 2: Stranded, Paired-end, Mouse](#case-2-stranded-paired-end-mouse)
  - [Case 3: Unstranded, Single-end, Mouse](#case-3-unstranded-single-end-mouse)
  - [Case 4: Stranded, Single-end, Mouse](#case-4-stranded-single-end-mouse)
  - [Case 5: Unstranded, Paired-end, Human](#case-5-unstranded-paired-end-human)
  - [Case 6: Stranded, Paired-end, Human](#case-6-stranded-paired-end-human)
  - [Case 7: Unstranded, Single-end, Human](#case-7-unstranded-single-end-human)
  - [Case 8: Stranded, Single-end, Human](#case-8-stranded-single-end-human)

## Case 1: Unstranded, Paired-end, Mouse

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCm39_vM26 --local_annot_dir <path to annotation> \
  --reads '*_{R1,R2}_001.fastq.gz' --outdir <path to outdir>
```

## Case 2: Stranded, Paired-end, Mouse

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCm39_vM26 --local_annot_dir <path to annotation> \
  --stranded fr-firststrand \
  --reads '*_{R1,R2}_001.fastq.gz' --outdir <path to outdir>
```

## Case 3: Unstranded, Single-end, Mouse

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCm39_vM26 --local_annot_dir <path to annotation> \
  --single_end \
  --reads '*_{R1}_001.fastq.gz' --outdir <path to outdir>
```

## Case 4: Stranded, Single-end, Mouse

```bash
nextflow run rikenbit/ramdaq -profile docker
  --genome GRCm39_vM26  --local_annot_dir <path to annotation> \
  --stranded fr-firststrand --single_end \
  --reads '*_{R1}_001.fastq.gz' --outdir <path to outdir>
```

## Case 5: Unstranded, Paired-end, Human

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCh38_v37 --local_annot_dir <path to annotation> ＼
  --reads '*_{R1,R2}_001.fastq.gz' --outdir <path to outdir>
```

## Case 6: Stranded, Paired-end, Human

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCh38_v37 --local_annot_dir <path to annotation> \
  --stranded fr-firststrand \
  --reads '*_{R1,R2}_001.fastq.gz' --outdir <path to outdir>
```

## Case 7: Unstranded, Single-end, Human

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCh38_v37 --local_annot_dir <path to annotation> \
  --single_end \
  --reads '*_{R1}_001.fastq.gz' --outdir <path to outdir>
```

## Case 8: Stranded, Single-end, Human

```bash
nextflow run rikenbit/ramdaq -profile docker \
  --genome GRCh38_v37 --local_annot_dir <path to annotation> \
  --stranded fr-firststrand --single_end \
  --reads '*_{R1}_001.fastq.gz' --outdir <path to outdir>
```
