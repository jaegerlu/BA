#!/usr/bin/env python3
import concurrent.futures
import sys
# sys.path.append('/home/j/jaegerl/.local/lib/python3.11/site-packages') #fixen iwann

import bisect
import subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
# from StringIO import StringIO
from Bio import AlignIO
from pymsaviz import MsaViz
import threading
import os  #


def read_genes(genes_path):
    genes = []  # [id,start_bp,end_bp,strand]
    with open(genes_path, 'r') as file:
        for line in file:
            genes.append(line.rstrip('\n').split('\t'))
    return genes


def read_gencode_annotation(annotation_path, just_protein_coding, gene_type_path):  # just_protein_coding = bool
    genes = []  # [id,start,end,strand,chrom]
    with open(annotation_path, 'r') as a, open(gene_type_path, 'w') as o:
        for line in a:
            if line.startswith('chr'):
                split = line.split('\t')
                if split[2] == 'gene':
                    infos = split[8].split(';')
                    for info in infos:
                        if info.startswith('gene_id'):
                            gene_id = info.split('"')[1]  # oder gene_name?
                        if info.startswith(' gene_type'):
                            gene_type = info.split('"')[1]
                            o.write(
                                gene_id + '\t' + gene_type + '\n')  ###### entferne nach einmaligem aufruf + in open()
                            if just_protein_coding and gene_type != 'protein_coding':
                                break
                    else:
                        chrom = split[0].split('r')[1]
                        start = split[3]
                        end = split[4]
                        strand = split[6]
                        if strand == '+':
                            strand = '1'
                        else:
                            strand = '-1'
                        genes.append([gene_id, start, end, strand, chrom])
    return genes


def read_gfa(gfa_path, chrom_names, segs_links_index_path, path_index_path, chrom_index_path):
    header = ''
    newChrom = True
    firstChrom = True
    segments = {}
    links = {}
    paths = {}
    ref_paths_indexes = {}  # chrom : [start,end]
    akt_chrom = '0'

    with open(gfa_path, 'r') as g, open(segs_links_index_path, 'w') as ind, open(
            path_index_path, 'w') as path_ind, open(chrom_index_path, 'w') as chr_ind:

        line = g.readline()
        tell = 0

        chrom_links = {}
        while line:
            if line.startswith('H'):
                header = line.rstrip('\n')
                ind.write(line)
            if line.startswith('S'):
                if newChrom:  ################
                    newChrom = False
                    akt_chrom = chrom_names.pop(0)
                    if not firstChrom:
                        paths[akt_chrom] = chrom_paths
                        chr_ind.write(akt_chrom + '\t' + str(tell) + '\n')
                        for key, value in chrom_links.items():
                            value_str = ','.join(value)
                            write_line = key + '\t' + value_str + '\n'
                            ind.write(write_line)
                        links.update(chrom_links)
                        chrom_links = {}
                    chrom_paths = {}

                split = line.split('\t')
                id = split[1]
                l = len(split[2].strip('\n'))  # length
                start = tell + len(id) + 3
                end = start + l
                segments[id] = [start, end]

                write_line = id + '\t' + str(start) + '\t' + str(end) + '\n'
                ind.write(write_line)  # + l

            if line.startswith('L'):
                if not newChrom:
                    newChrom = True
                split = line.split('\t')
                first = split[1]
                second = split[3]
                if first in chrom_links:
                    if second not in chrom_links[first]:
                        chrom_links[first].add(second)
                else:
                    chrom_links[first] = {second}
                # print(chrom_links)

            if line.startswith('P'):
                split = line.split('\t')
                id = split[1]
                start = tell + len(id) + 3
                end = start + len(split[2])
                # print(id)
                if id.startswith('grch38'):
                    ref_paths_indexes[akt_chrom] = [start, end]

                chrom_paths[id] = [start, end]  #

                write_line = id + '\t' + str(start) + '\t' + str(end) + '\t' + akt_chrom + '\n'
                path_ind.write(write_line)

            tell = g.tell()
            line = g.readline()

        paths[akt_chrom] = chrom_paths
        links.update(chrom_links)
        for key, value in chrom_links.items():
            value_str = ','.join(value)
            write_line = key + '\t' + value_str + '\n'
            ind.write(write_line)
        chr_ind.write(akt_chrom + '\t' + str(tell))

    return [header, segments, links, paths, ref_paths_indexes]


def read_Refpathseg_path_hash(hash_file):
    seg_contained_by_path = {}
    with open(hash_file, 'r') as h:
        for line in h:
            line = line.rstrip('\n').split('\t')
            seg_contained_by_path[line[0]] = set(line[1].split(','))
    return seg_contained_by_path


def read_index_file(index_file):
    header = ''
    segments = {}
    links = {}

    with open(index_file, 'r') as ind:
        for line in ind:
            split = line.rstrip('\n').split('\t')
            if split[0] == 'H':
                header = line.rstrip('\n')
            elif len(split) == 3:
                segments[split[0]] = [int(split[1]), int(split[2])]
            elif len(split) == 2:
                seg_to = split[1].split(',')
                links[split[0]] = seg_to

    return [header, segments, links]


def read_path_ind(path_ind_file):
    paths = {}  # chrom1 : {id: start,end}
    chrom_paths = {}
    aktchrom = 'chrom1'
    ref_paths_indexes = {}  # chrom : [start,end]
    with open(path_ind_file, 'r') as f:
        for line in f:
            split = line.rstrip('\n').split('\t')
            id = split[0]
            start = int(split[1])
            end = int(split[2])
            chrom = split[3]
            if chrom != aktchrom:
                paths[aktchrom] = chrom_paths
                chrom_paths = {}
                aktchrom = chrom
            if id.startswith('grch38'):
                ref_paths_indexes[chrom] = [int(start), int(end)]
            chrom_paths[id] = [int(start), int(end)]
        paths[aktchrom] = chrom_paths
    return [paths, ref_paths_indexes]


def get_ref_paths(gfa_file, ref_paths_indexes):
    ref_paths = {}
    with open(gfa_file, 'rb') as g:
        for chrom, pos in ref_paths_indexes.items():
            g.seek(pos[0])
            path = g.read(pos[1] - pos[0]).decode('utf8').split(',')  # schon gesplittet
            ref_paths[chrom] = path
    return ref_paths


def create_seg_to_paths_map(gfa_path, ref_segs, paths, seg_hash_path):
    seg_contained_by_path = {}  # segid: [path1,path2]
    lock = threading.Lock()
    path_sets = {}
    with open(gfa_path, 'rb') as g:
        for chrom, chrom_paths in paths.items():
            path_chrom_sets = {}
            for pid, pval in chrom_paths.items():
                g.seek(int(pval[0]))
                path = g.read(int(pval[1]) - int(pval[0])).decode('utf8').rstrip('\n')
                path_set = set()
                for pseg in path.split(','):
                    path_set.add(pseg[:-1])
                path_chrom_sets[pid] = path_set
            path_sets[chrom] = path_chrom_sets
    print('Paths to set ready')
    # ref_segs = ref_segs[:10000]  #####################für test
    print(ref_segs)

    def find_paths_for_segment(segment, chrom):
        #segment = segment[:-1]
        seg_in_path = set()
        for pid, pval in path_sets[chrom].items():
            if segment in pval:
                seg_in_path.add(pid)

        with lock:
            seg_contained_by_path[segment] = seg_in_path
        print(segment + ' fertig')

    for chrom in path_sets.keys():
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for ref_seg in ref_segs:
                executor.submit(find_paths_for_segment, ref_seg, chrom)
            # executor.map(find_paths_for_segment, ref_segs, chrom)

    with open(seg_hash_path, 'w') as r:
        for seg, pathss in seg_contained_by_path.items():
            r.write(seg + '\t' + ','.join(pathss) + '\n')

    return [seg_contained_by_path, path_sets]  # path_sets


def get_ref_seg_lens(ref_path,
                     segments):  # trägt alle daten zum ref path ein und gibt sortierte längen für bisect zurück
    # print(ref_path)
    ref_len_rev = {}
    ref_lens = {}
    #ref_dirs = {}
    for chrom, path in ref_path.items():
        kum_len = 0
        ref_lens_chrom = {}
        ref_len_rev_chrom = {}
        #ref_dirs_chrom = {}
        for ref in path:
            s_id = ref[:-1]
            #s_dir = ref[-1]
            length = segments[s_id][1] - segments[s_id][0]
            ref_lens_chrom[s_id] = kum_len
            ref_len_rev_chrom[kum_len] = s_id
            kum_len += length
        ref_lens[chrom] = ref_lens_chrom
        ref_len_rev[chrom] = ref_len_rev_chrom

    sorted_lens = {}
    # print(ref_len_rev.keys())
    for chrom, refrev in ref_len_rev.items():
        sorted_lens[chrom] = sorted(refrev.keys())

    return [sorted_lens, ref_lens, ref_len_rev]


def reverse_path(path):
    out = ''
    dirs = {'+':'-','-':'+'}
    for p in path.split(','):
        seg = p[:-1]
        dir = p[-1]
        out = seg+dirs[dir]+','+out
    return out[:-1]


def create_gene_index_file(gfa_file, paths, genes, sorted_lens, ref_lens, ref_len_rev, gene_index_path, seg_hash_path,
                           e):
    complete_genes = []
    segs_need_to_be_hashed = []
    appeared_chroms = []
    for gene in genes:
        [id, start, end, strand, chrom] = gene

        if int(start) > int(end):  # switche, bei opposite strand manchmal end kleiner als start
            temp = start
            start = end
            end = temp

        if chrom != '1':  ####################test für firstChrom
            continue

        new_index = bisect.bisect_right(sorted_lens[chrom], int(start))  # finde start und end segment
        if new_index:
            start_seg = ref_len_rev[chrom][sorted_lens[chrom][new_index - 1]]
        else:
            print('Fehler bei new_index')

        new_index = bisect.bisect_right(sorted_lens[chrom], int(end))
        if new_index:
            end_seg = ref_len_rev[chrom][sorted_lens[chrom][new_index - 1]]
        else:
            print('Fehler bei new_index')

        if strand == '1':
            segs_need_to_be_hashed.append(start_seg)
        else:
            segs_need_to_be_hashed.append(end_seg)

        start_g = int(start) - ref_lens[chrom][
            start_seg] - 1  # startposition des Gens innerhalb des Startsegments #### -1 wsh weg?
        end_g = int(end) - ref_lens[chrom][end_seg]  # +1 evtl weg bei ende von gen excludet/includet

        complete_genes.append([id, start, end, strand, chrom, start_seg, end_seg, start_g, end_g])
        if chrom not in appeared_chroms:
            appeared_chroms.append(chrom)

    if e:
        seg_contained_by_path = read_Refpathseg_path_hash(seg_hash_path)
    else:
        [seg_contained_by_path, path_sets] = create_seg_to_paths_map(gfa_file, segs_need_to_be_hashed, paths,
                                                                     seg_hash_path)
    print('hashing of segments ready')


    gene_segs = {}
    lock = threading.Lock()

    if os.path.exists(gene_index_path):
        # Lösche die Datei
        os.remove(gene_index_path)


    def create_gene_index_per_gene(id, strand, chrom, start_seg, end_seg, start_g, end_g, g):
        out = '>' + '\t'.join(str(x) for x in gene) + '\n'

        if start_seg in set(segs_need_to_be_hashed):
            print(start_seg+' in hashed')
        if  end_seg in set(segs_need_to_be_hashed):
            print(start_seg + ' in hashed')

        todo = [start_seg]
        completed = set()

        skip_last_seg = False
        if end_g == 0 and strand == '1':
            skip_last_seg = True
        if start_g == 0 and strand == '-1':
            skip_last_seg = True

        while todo:
            # print(todo)
            seg = todo.pop(0)

            if seg == end_seg:
                completed.add(seg)
                continue

            if seg in links:
                for n in links[seg]:
                    if int(n) <= int(
                            end_seg) and n not in completed and n not in todo:  ### int(id) Vergleich!
                        todo.append(n)
            completed.add(seg)

        if skip_last_seg and strand == '1':
            completed.remove(end_seg)
        elif skip_last_seg and strand == '-1':
            completed.remove(start_seg)

        already_contained_samples = set()
        if strand == '1':
            beginning_seg_dir = start_seg # + '+' ####plus weg?
        else:
            gene_segs[id] = completed  #########muss auch noch für index einlesung erstellt werden. neeee
            print('>' + '\t'.join(str(x) for x in gene) + '\n')
            beginning_seg_dir = end_seg # + '-'

        if beginning_seg_dir in seg_contained_by_path:
            print(seg_contained_by_path[beginning_seg_dir])
        else:
            print('Fehler, nicht enthalten')
        sample_paths = []
        for p_id in seg_contained_by_path[beginning_seg_dir]:
            if p_id == '':
                print(seg_contained_by_path[beginning_seg_dir])
                continue
            p_data = paths[chrom][p_id]

            p_id_split = p_id.split('#')
            sample_id = p_id_split[0] + '#' + p_id_split[1]
            # print(p_id + '\t' + str(smallseg) + '\t' + str(bigseg))
            if sample_id not in already_contained_samples:
                g.seek(int(p_data[0]))
                p_path = g.read(int(p_data[1]) - int(p_data[0])).decode('utf8').rstrip('\n').split(',')
#########schreibe dir raus, ob das startseg vom refpath + oder - ist. Falls anderer path andersum ist, drehe path um
                out_path = ''
                first_dir = p_path[0][-1]
                for p_seg in p_path:

                    if p_seg[:-1] in completed:
                        '''
                        dir = p_seg[-1]
                        if (strand == '1' and dir == '+') or (
                                strand == '-1' and dir == '-') or dir != first_dir:  # letzte prüfung, falls ein segment wirklich umgedreht ist
                            out_path = out_path + p_seg + ','
                        elif strand == '1' and dir == '-':  # überprüfe ob es pfade gibt, die sowohl + als auch - enthalten: dann entscheide nach 'allgemeiner Richtung und drehe um'
                            out_path = p_seg[:-1] + '+,' + out_path
                        elif strand == '-1' and dir == '+':
                            out_path = p_seg[:-1] + '-,' + out_path
                        else:
                            print('Fehler bei überprüfung von direction und strand')
                        '''
                        out_path = out_path + p_seg + ','

                out_path = out_path[:-1]
                if skip_last_seg and strand == '1':
                    out_path = ','.join(out_path.split(',')[:-1])
                if skip_last_seg and strand == '-1':
                    out_path = ','.join(out_path.split(',')[:-1])

                #if out_path.split(',')[0][:-1] == end_seg:      # der read kommt von dem opposite strang hmm ne warte
                #    out_path = reverse_path(out_path)

                sample_paths.append([sample_id,out_path])
                if sample_id.startswith('grch'):
                    ref_dir = out_path.split(',')[0][-1]
                already_contained_samples.add(sample_id)

        for s_path in sample_paths: #drehe reads um, wenn anderer strang von ref, genom oder strand = -1
            first_dir = s_path[1].split(',')[0][-1]
            if (first_dir == ref_dir and strand == '1') or (first_dir!=ref_dir and strand == '-1'):
                out += s_path[0]+'\t'+s_path[1]+'\n'
            else:
                print('reverse_path wird aufgerufen')
                out += s_path[0]+'\t'+reverse_path(s_path[1])+'\n'
        with lock:
            # g_ind.write(out)
            with open(gene_index_path, 'a') as g_ind:
                g_ind.write(out)
            print(str(os.path.getsize(gene_index_path))+' Bytes')
            # index_file_output.append(out)
            # gene_segs[id] = completed  #########muss auch noch für index einlesung erstellt werden. neeee
            # print('>' + '\t'.join(str(x) for x in gene) + '\n')
            #print(out)

    with open(gfa_file, 'rb') as g:  # open('GGE_genes.index', 'w') as g_ind,
        for gene in complete_genes:
            [id, start, end, strand, chrom, start_seg, end_seg, start_g, end_g] = gene

            with concurrent.futures.ThreadPoolExecutor() as executor:
                executor.submit(create_gene_index_per_gene, id, strand, chrom, start_seg, end_seg, start_g, end_g, g)
        # for gene in index_file_output:
        # g_ind.write(gene)


def reverse(seq):
    revcomp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    revcomp = ""
    for nt in reversed(seq):
        if nt == 'N':
            revcomp += 'N'
        elif nt == '-':
            revcomp += '-'
        else:
            revcomp += revcomp_dict[nt]
    return revcomp


def printrefSeq(gfa, strand):
    with open(gfa, 'r') as g:
        segments = {}  # id:seq
        for line in g:
            if line.startswith('S'):
                split = line.split('\t')
                segments[split[1]] = split[2][:-1]
            if line.startswith('P\tgrch'):
                split = line.split('\t')
                seq = ''
                for seg in split[2].split(','):
                    segid = seg[:-1]
                    dir = seg[-1]
                    if dir == '+':
                        seq = seq + segments[segid]
                    else:
                        seq = seq + reverse(segments[segid])
                if strand == '1':
                    print(seq)
                else:
                    print(reverse(seq))
                print('refSeq-len ' + str(len(seq)))
                break


def print_path(segments, path, gfa_file):
    out = ''
    with open(gfa_file, 'rb') as g:
        for p in path:
            p_seg = p[:-1]
            dir = p[-1]
            seg = segments[p_seg]
            g.seek(seg[0])
            seq = g.read(seg[1] - seg[0]).decode('utf8')
            if dir == '+':
                out += seq
            else:
                out += reverse(seq)
    print(out)


def readVCFfile(vcf_path,
                haplos):  # schaue falls gar keine varianten, haplos ist eine Liste aller haplotypesamples, eins por haplotyp

    ab_hier = False
    haplo_order = ['grch38']
    haplo_pos = []
    variants = []
    with open(vcf_path, 'r') as vcf:
        for line in vcf:
            if ab_hier:
                variants.append(line.rstrip('\n'))
            if line.startswith('#CHROM'):
                ab_hier = True
                split = line.rstrip('\n').split('\t')
                for j in range(9, len(split)):
                    for hap in haplos:
                        if hap == split[j]:
                            haplo_order.append(
                                hap)  ##########speichere auch die index position wo die nicht haplo samples eingerechnet sind
                            haplo_pos.append(j - 9)

    return [haplo_order, variants, haplo_pos]


def createMSA(gene_vcf_path, gene_ref_path, segments, path_seg_sets, pangenome, clustalw_path, haplos, start_g, end_g,
              strand):
    # erstelle haplos #siehe oben
    varianten_hash_id = {}
    varianten_hash_seg = {}  # seg:id
    msa = []  # [[*,seq]]
    # *: a für 'in allen haplos'
    #   b für 'nicht in allen sampeln' z.b. b:s,.,s,s,s,.,s #0 statt s
    #   o für 'variante aber len 1 und damit nicht neu aligniert'
    #   r für neu aligniert mit clustalw

    variant = set()  # segmente die durch vcf schon abgelaufen wurden

    [haplo_order, variants, haplo_pos] = readVCFfile(gene_vcf_path, haplos)

    ref_pos = 0
    n_haplo = len(haplos)

    gene_ref_path_split = gene_ref_path.split(',')
    start_seg = gene_ref_path_split[0][:-1]
    end_seg = gene_ref_path_split[-1][:-1]

    with open(pangenome, 'rb') as g:
        for line in variants:
            split = line.split('\t')
            var_id = split[2]
            # ref = split[3]
            # alt = split[4]
            print(split)
            startseg_data = split[7].split('>')
            if len(startseg_data) > 1:
                startseg = split[7].split('>')[1]
            else:
                startseg = split[7].split('<')[1]
            samples = split[9:]

            varianten_hash_id[var_id] = split  # für nested snarls. speichere nür nötiges. Gerade: alles lol
            varianten_hash_seg[startseg] = var_id

        print(gene_ref_path)  ##

        for ref_seg in gene_ref_path_split:
            ref_seg = ref_seg[:-1]
            if ref_seg in variant:
                variant.remove(ref_seg)
                continue

            in_sample = set()
            for haplo in haplos:
                if ref_seg in path_seg_sets[haplo]:
                    in_sample.add(haplo)
            seg_index = segments[ref_seg]
            g.seek(seg_index[0])
            seq = g.read(seg_index[1] - seg_index[0]).decode('utf-8')
            seq_len = len(seq)

            if ref_seg == start_seg:
                seq = seq[-(seq_len - start_g):]
            if ref_seg == end_seg:
                seq = seq[:-(seq_len - end_g)]

            if len(in_sample) == n_haplo:
                msa.append(['a', seq])
            else:
                b = 'b:' + give_sample_status(samples, haplo_pos)
                msa.append([b, seq])

            # print(varianten_hash_seg)
            # print('\n\n\n\n\n\n\n\n')
            # print(varianten_hash_id)
            if ref_seg in varianten_hash_seg:
                split = varianten_hash_id[varianten_hash_seg[ref_seg]]
                infos = split[7].split(';')
                ref = split[3]
                alt = split[4].split(',')
                samples = split[9:]
                for info in infos:
                    if info.startswith('AT'):
                        alleles = info[4:].split(',')
                        for allel in alleles:
                            variant.update(allel.split('>')[1:-1])
                    if info.startswith('LV'):
                        LV = info.split('=')[1]
                        # if LV != 0:        nested snarl
                allaltlen1 = True
                if len(ref) == 1:
                    for a in alt:
                        if not len(a) == 1:
                            allaltlen1 = False
                    if allaltlen1:
                        msa.append(['o:' + give_sample_status(samples, haplo_pos), [ref, ','.join(alt)]])
                if len(ref) > 1 or not allaltlen1:  # neues Alignement
                    flag = 'r:' + give_sample_status(samples, haplo_pos)
                    msa.append([flag, call_clustalw(ref, alt, gene_vcf_path, clustalw_path)])

    ##create msa_output

    msa_path = gene_vcf_path[:-3] + 'pdf'
    temp_path = gene_vcf_path[:-3] + 'temp'
    haplos_seqs = {}
    markers = []
    aktlen = 0
    for haplo in haplo_order:
        haplos_seqs[haplo] = ''
    for block in msa:
        flag = block[0]
        block_seqs = block[1]
        if flag == 'a':
            for id, seq in haplos_seqs.items():
                haplos_seqs[id] = seq + block_seqs
            aktlen += len(block_seqs)
        elif flag.startswith('b'):
            block_len = len(block_seqs)
            block_flag = flag[2:].split(',')
            for i, haplo in enumerate(haplo_order):
                if block_flag[i] == 0:
                    haplos_seqs[haplo] = haplos_seqs[haplo] + block_seqs
                else:
                    haplos_seqs[haplo] = haplos_seqs[haplo] + ('-' * block_len)
            aktlen += block_len
        elif flag.startswith('o'):
            markers.append(['o', aktlen + 1])
            block_flag = flag[2:].split(',')
            for i, haplo in enumerate(haplo_order):
                if block_flag[i] == '.':
                    haplos_seqs[haplo] = haplos_seqs[haplo] + '-'
                else:
                    haplos_seqs[haplo] = haplos_seqs[haplo] + block_seqs[int(block_flag[i])]
            aktlen += 1
        elif flag.startswith('r'):
            block_flag = flag[2:].split(',')
            block_seqs = [record.seq for record in block_seqs]
            block_len = len(block_seqs[0])
            for i, haplo in enumerate(haplo_order):
                if block_flag[i] == '.':
                    haplos_seqs[haplo] = haplos_seqs[haplo] + '-' * block_len
                else:
                    haplos_seqs[haplo] = haplos_seqs[haplo] + block_seqs[int(block_flag[i])]
            markers.append(['r', [aktlen, aktlen + block_len]])
            aktlen += block_len

    with open(temp_path, 'w') as temp:
        for haplo in haplo_order:
            temp.write('>' + haplo + '\n')
            #print(haplos_seqs[haplo])
            temp.write(str(haplos_seqs[haplo]) + '\n')

            # print(haplos_seqs[haplo])

    #print('until now')
    mv = MsaViz(temp_path, wrap_length=60, color_scheme='None')
    os.remove(temp_path)
    #print('until 1')
    # add markers
    for mark in markers:
        flag = mark[0]
        pos = mark[1]
        if flag == 'o':
            mv.add_markers([pos], color='orange', marker='o')
        # if flag == 'r':
        # mv.add_text_annotation((pos[0], pos[1]), text_color='red', range_color='red', text='aligned')
    #print('until 2')
    mv.savefig(msa_path)
    #print('until then')

    # *: a für 'in allen haplos'
    #   b für 'nicht in allen sampeln' z.b. b:s,.,s,s,s,.,s #0 statt s
    #   o für 'variante aber len 1 und damit nicht neu aligniert'
    #   r für neu aligniert mit clustalw


def give_sample_status(samples, haplo_pos):  # z.b. 0,2,1,1,1,0,.,1,0 nur für haplos
    output = ['0']
    for pos in haplo_pos:
        output.append(samples[pos])
    return ','.join(output)


def call_clustalw(ref, alt, gene_vcf_path, clustalw_path):
    fasta_file = gene_vcf_path[:-3] + 'fas'
    alignment_file = gene_vcf_path[:-3] + 'aln'

    with open(fasta_file, 'w') as f:
        f.write('>0\n' + ref + '\n')
        for count, seq in enumerate(alt):
            f.write('>' + str(count + 1) + '\n' + seq + '\n')

    clustal_command = '/'.join(clustalw_path.split('/')[:-1]) + '/./' + clustalw_path.split('/')[-1]
    clustalw_arguments = ['-infile=' + fasta_file, '-outfile=' + alignment_file, '-GAPOPEN=1', '-GAPEXT=0','-numiters','1']
    process = subprocess.Popen([clustal_command] + clustalw_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()

    if process.returncode != 0:
        print("Calling ClustalW failed")
    # assert process != 0, "Calling ClustalW failed"

    print(ref)
    print(alt)
    print(process.returncode)

    alignment = AlignIO.read(alignment_file, 'clustal')
    print(str(alignment))
    os.remove(fasta_file)
    os.remove(alignment_file)
    return alignment


def create_xg(genegfa, vg_path, strand):
    command = vg_path[:-2] + './vg'
    xg_file = genegfa[:-3] + 'xg'
    xg_arguments = ['index', genegfa, '-x', xg_file, '-t', '32', '-p']
    process = subprocess.Popen([command] + xg_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stderr.decode().strip())

    vcf_file = genegfa[:-3] + 'vcf'
    vcf_arguments = ['deconstruct', '-t', '32', '-a', '-e', '-v', '-P', 'grch38', genegfa]
    with open(vcf_file, 'w') as vcf:
        process = subprocess.Popen([command] + vcf_arguments, stdout=vcf, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stderr.decode().strip())

    fasta_file = genegfa[:-3] + 'fasta'
    fasta_arguments = ['paths', '-x', xg_file, '-F']
    process = subprocess.Popen([command] + fasta_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    # print(stdout.decode().strip())
    # print(stderr.decode().strip())
    fasta = stdout.decode().strip().split('>')[1:]
    fasta_len = len(fasta)
    fasta_dic = {}  # ATCCG: [1#2, 2#1, 26#2]
    # print(fasta)
    for fa in fasta:
        fa_split = fa.split('\n')
        fa_id = fa_split[0]
        fa_seq = ''.join(fa_split[1:])

        if fa_seq in fasta_dic:
            fasta_dic[fa_seq].append(fa_id)
        else:
            fasta_dic[fa_seq] = [fa_id]
    records = []
    for fa_seq, fa_ids in fasta_dic.items():
        sequence = Seq(fa_seq)

        haplo_sample = fa_ids[0]
        for id in fa_ids:
            if id.startswith('grch38'):
                haplo_sample = id
        record = SeqIO.SeqRecord(sequence, id=str(len(fa_ids)) + '/' + str(fasta_len) + ':' + haplo_sample,
                                 description='')  # id=str(len(fa_ids))+'/'+str(fasta_len)+':'+','.join(fa_ids)
        records.append(record)
    with open(fasta_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    haplos = []
    for rec in records:
        if not rec.id.split(':')[1].startswith('grch38'):
            haplos.append(rec.id.split(':')[1])

    return [haplos, vcf_file]


### mini gfa erstellung fehlt noch

def read_gene_index(file_path):
    paths = []
    first = True
    gene_index_data = []
    with open(file_path, 'r') as g:
        for line in g:
            if line.startswith('>'):
                if first:
                    first = False
                else:
                    gene_index_data.append([id, strand, chrom, start_seg, end_seg, start_g, end_g, paths])
                    paths = []

                [id, start, end, strand, chrom, start_seg, end_seg, start_g, end_g] = line.rstrip('\n').lstrip(
                    '>').split('\t')
                start_g = int(start_g)
                end_g = int(end_g)
            else:
                if line:
                    paths.append(line.rstrip('\n'))
    return gene_index_data


def create_mini_gfa(pangenome, header, segments, links, id, startseg, endseg, strand, start_g, end_g, directory,
                    g_paths):
    with open(pangenome, 'rb') as g, open(directory + id + '.gfa', 'w') as o:
        segs = set()
        o.write(header + '\n')
        for g_path in g_paths:  # get all segs
            g_path = g_path.rstrip('\n')
            for g_p_seg in g_path.split('\t')[1].split(','):
                g_p_seg = g_p_seg[:-1]
                segs.add(g_p_seg)
        for seg in segs:  # write all segs
            seg_data = segments[seg]
            g.seek(seg_data[0])
            seq = g.read(seg_data[1] - seg_data[0]).decode('utf8')
            seq_len = len(seq)
            if seg == startseg:
                if strand == '1':
                    seq = seq[-(seq_len - start_g):]
                    # print(seq)
                else:
                    seq = seq[:-(seq_len - start_g)]
                seq_len -= start_g
            if seg == endseg:
                if strand == '1':
                    seq = seq[:-(seq_len - end_g)]
                else:
                    seq = seq[-(seq_len - end_g):]
            o.write('S\t' + seg + '\t' + seq + '\n')
        for path in g_paths:    #write all paths
            o.write('P\t' + path + '\t*\n')
        for from_s in segs: #write all links
            if from_s in links:
                for to_s in links[from_s]:
                    if to_s in segs:
                        o.write('L\t' + from_s + '\t+\t' + to_s + '\t+\t*\n')
            else:
                print(from_s+' nicht in links enthalten')
        print(chrom+' '+id+'.gfa erstellt')

        return o.name


def search_ref_path_and_get_path_set_sets(g_paths):
    ref_path = ''
    path_set_sets = {}
    for g_path in g_paths:
        split = g_path.rstrip('\t').split('\t')
        if g_path.startswith('grch38'):
            ref_path = split[1]
        p_set = set()
        for seg in split[1].split(','):
            p_set.add(seg[:-1])
        path_set_sets[split[0]] = p_set
    return [ref_path, path_set_sets]


def create_mini_gfa_and_etc(pangenome, header, segments, links, id, startseg, endseg, starnd, start_g, end_g, directory,
                            g_paths, vg_path, clustalw_path):
    print('etc aufgerufen')
    minigfa_path = create_mini_gfa(pangenome, header, segments, links, id, startseg, endseg, starnd, start_g, end_g,
                                   directory, g_paths)
    [haplos, vcf_file] = create_xg(minigfa_path, vg_path, strand)
    [ref_path, path_set_sets] = search_ref_path_and_get_path_set_sets(g_paths)
    createMSA(vcf_file, ref_path, segments, path_set_sets, pangenome, clustalw_path, haplos, start_g, end_g, strand)
    print(id + ' creations fertig')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extracts genetic information out of pangenome')
    parser.add_argument('-p', '--gfa_path')
    parser.add_argument('-g', '--genes_path')
    parser.add_argument('-v', '--vg_path')
    parser.add_argument('-w', '--clustalw_path')
    parser.add_argument('-i', '--generate_index', action='store_true')
    parser.add_argument('-a', '--annotation_path')
    parser.add_argument('-e', '--designer_path', action='store_true')

    args = parser.parse_args()

    gfa_path = args.gfa_path
    genes_path = False
    if args.genes_path:
        genes_path = args.genes_path  # simuliert die eingabe
    annotation_path = False
    if args.annotation_path:
        annotation_path = args.annotation_path
    vg_path = args.vg_path
    clustalw_path = args.clustalw_path
    generate_index = args.generate_index
    designer_path = args.designer_path

    numbers = [str(i + 1) for i in range(22)]
    additional_elements = ['X', 'Y', 'M']
    chrom_names = numbers + additional_elements

    just_protein_coding = True

    gene_index_path = 'GGE_genes.index'
    segs_links_index_path = 'GGE_segs+links.index'
    path_index_path = 'GGE_path.index'
    chrom_index_path = 'GGE_chrom.index'
    gene_type_path = 'GGEpy_gene_type_per_gene'
    seg_hash_path = 'GGE_refpathseg_path_hash.txt'
    if genes_path:
        genes = read_genes(genes_path)
    elif annotation_path:
        genes = read_gencode_annotation(annotation_path,just_protein_coding, gene_type_path)
    print('genes eingelesen')

    if designer_path:
        [header, segments, links] = read_index_file(segs_links_index_path)
        print('segs und links eingelesen')
        [paths, ref_paths_indexes] = read_path_ind(path_index_path)
        print('paths eingelesen')
        ref_paths = get_ref_paths(gfa_path, ref_paths_indexes)
        print('ref paths geholt')
        [sorted_lens, ref_lens, ref_len_rev] = get_ref_seg_lens(ref_paths, segments)
        print('ref seg lens berechnet')

        create_gene_index_file(gfa_path, paths, genes, sorted_lens, ref_lens, ref_len_rev, gene_index_path,
                               seg_hash_path, True)
        print('gene_index_file_created')

    elif generate_index:
        [header, segments, links, paths, ref_paths_indexes] = read_gfa(gfa_path, chrom_names, segs_links_index_path,
                                                                       path_index_path, chrom_index_path)
        print('gfa-file finished')
        ref_paths = get_ref_paths(gfa_path, ref_paths_indexes)
        print('ref path geholt')
        # [seg_contained_by_path, path_sets] = create_seg_to_paths_map(gfa_path,ref_paths,paths)
        # print('refseg_seg_maps fertig')
        [sorted_lens, ref_lens, ref_len_rev] = get_ref_seg_lens(ref_paths, segments)
        print('ref seg lens berechnet')
        create_gene_index_file(gfa_path, paths, genes, sorted_lens, ref_lens, ref_len_rev, gene_index_path,seg_hash_path, False)
        print('gene_index_file_created')
    else:
        [header, segments, links] = read_index_file(segs_links_index_path)
        print('segs und links eingelesen')
        [paths, ref_paths_indexes] = read_path_ind(path_index_path)
        print('paths eingelesen')
        ref_paths = get_ref_paths(gfa_path, ref_paths_indexes)
        print('ref paths geholt')
        seg_contained_by_path = read_Refpathseg_path_hash(seg_hash_path)
        print('ref_path_hash eingelesen')
        [sorted_lens, ref_lens, ref_len_rev] = get_ref_seg_lens(ref_paths, segments)
        print('ref seg lens berechnet')

    gene_index_data = read_gene_index(gene_index_path)
    print('gene_index eingelesen')

    chrom_dics = {}  # erstelle Chrom directories if they dont exist yet
    for chrom in chrom_names:
        chrom_dic = 'chrom' + chrom
        chrom_dics[chrom] = chrom_dic
        if not os.path.exists(chrom_dic):
            os.makedirs(chrom_dic)

    for g_i_data in gene_index_data:
        [id, strand, chrom, start_seg, end_seg, start_g, end_g, g_paths] = g_i_data
        directory = chrom_dics[chrom]+'/'

        create_mini_gfa_and_etc(gfa_path, header, segments, links, id, start_seg, end_seg, strand,
        start_g, end_g, directory, g_paths, vg_path, clustalw_path)###################test
        '''
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(create_mini_gfa_and_etc, gfa_path, header, segments, links, id, start_seg, end_seg, strand,
                            start_g, end_g, directory, g_paths, vg_path, clustalw_path)
        '''
