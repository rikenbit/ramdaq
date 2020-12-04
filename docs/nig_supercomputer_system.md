# Using ramdaq on the NIG Supercomputer System

[NIG Supercomputer System](https://sc.ddbj.nig.ac.jp/en?set_language=en) is the supercomputer system hosted by the National Institute of Genetics, Japan. On the NIG Supercomputer System, users can use [Singularity](https://sylabs.io/)([Usage in Japanese](https://sc.ddbj.nig.ac.jp/ja/guide/software/singularity)) and [Univa Grid Engine (UGE)](https://www.univa.com/products/) ([Usage in Japanese](https://sc.ddbj.nig.ac.jp/ja/guide/software/univa-grid-engine))

In this document, we offer the following instructions:

1. How to install nextflow on NIG Supercomputer System
2. How to run ramdaq on NIG Supercomputer System
   1. Run ramdaq on a node
   2. Run ramdaq using Univa Grid Engine (UGE)

## 1. Install nextflow on NIG Supercomputer System

1\. qlogin

```bash
qlogin -l s_vmem=20G -l mem_req=20G
```

2\. Type this for Java use

```bash
export _JAVA_OPTIONS="-XX:ParallelGCThreads=1 -Xmx1g"
```

3\. Check Java version (OK if version number isZ "1.8.y_z" or "8")

```bash
java -version
```

4\. Install nextflow

```bash
curl -s https://get.nextflow.io | bash
```

5\. Test nextflow

```bash
./nextflow run hello
```

6\. Add the file path of nextflow to $PATH

## 2-1. Run ramdaq on a node

1\. Login to a node by `qlogin`

```bash
qlogin -l s_vmem=20G -l mem_req=20G
```

2\. Type this for Java use

```bash
export _JAVA_OPTIONS="-XX:ParallelGCThreads=1 -Xmx1g"
```

3\. Run ramdaq (Below command run ramdaq on the test data)

```bash
nextflow run rikenbit/ramdaq -profile test,singularity
```

## 2-2. Run ramdaq by Univa Grid Engine (UGE)

1\. Prepare a job script (See the example below)

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -l s_vmem=20G -l mem_req=20G
#$ -cwd

export _JAVA_OPTIONS="-XX:ParallelGCThreads=1 -Xmx1g"

nextflow run rikenbit/ramdaq -profile test,singularity
```

2\. Submit a job by `qsub` command

```bash
qsub -pe def_slot 4 -l short qsub_test.sh
```

> **Note**: We empirically recommend to set the number of threads (cores) to be **at least `max_cpus` + 1** when running ramdaq via UGE.
>
> For example, the test config set `max_cpus = 2` and we add the option `-pe def_slot 3` on `qsub`.


## Acknowlegements

- [Yasuhiro Tanizawa](https://github.com/nigyta)
- [Manabu Ishii](https://github.com/manabuishii)
- [Tomoya Tanjo](https://github.com/tom-tan)
- [Workflow Meetup](https://github.com/workflow-meetup-jp/workflow-meetup)