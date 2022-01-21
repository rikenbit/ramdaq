#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO: Add additional regexes for new tools in process get_software_versions
regexes = {
    'ramdaq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Fastq-mcf': ['v_fastqmcf.txt', r"Version: (\S+)"],
    'Hisat2': ['v_hisat2.txt', r"/opt/conda/envs/ramdaq-1.0dev/bin/hisat2-align-s version (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Bam2Wig': ['v_bam2wig.txt', r"bam2wig.py (\S+)"],
    'Bamtools': ['v_bamtools.txt', r"bamtools (\S+)"],
    'read_distribution': ['v_read_distribution.txt', r"read_distribution.py (\S+)"],
    'infer_experiment': ['v_infer_experiment.txt', r"infer_experiment.py (\S+)"],
    'inner_distance': ['v_inner_distance.txt', r"inner_distance.py (\S+)"],
    'junction_annotation': ['v_junction_annotation.txt', r"junction_annotation.py (\S+)"],
    'FeatureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'RSEM': ['v_rsem.txt', r"Current version: RSEM v(\S+)"],
    'R': ['v_R.txt', r"R version (\S+)"],
    'edgeR': ['v_edgeR.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['ramdaq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Fastq-mcf'] = '<span style="color:#999999;\">N/A</span>'
results['Hisat2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['Bam2Wig'] = '<span style="color:#999999;\">N/A</span>'
results['Bamtools'] = '<span style="color:#999999;\">N/A</span>'
results['read_distribution'] = '<span style="color:#999999;\">N/A</span>'
results['infer_experiment'] = '<span style="color:#999999;\">N/A</span>'
results['inner_distance'] = '<span style="color:#999999;\">N/A</span>'
results['junction_annotation'] = '<span style="color:#999999;\">N/A</span>'
results['FeatureCounts'] = '<span style="color:#999999;\">N/A</span>'
results['RSEM'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['edgeR'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'ramdaq Software Versions'
section_href: 'https://github.com/rikenbit/ramdaq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
