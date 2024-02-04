Workflow consists of three tools pan2prep.py, pan2index.py and pan2gene.py.
Can extract gene data out of pangenome in GFA format.

Commands to run:
python3 pan2prep.py -p <GFA-pangenome>
python3 pan2index.py -p <GFA-pangenome> -a <gene annotation> -j -r 100 -t 6
python3 pan2gene.py -p <GFA-pangenome> -v <vg path> -a <gene annotation> -j

used pangenome file: hprc-v1.0-pggb.gfa
used annotation file: gencode.v38.annotation.gtf

Note that pan2 only runs on gfa files like this containing grch38 as reference genome.
