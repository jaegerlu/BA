#!/usr/bin/env python3
import os
import argparse
import pickle
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq



def get_genes(genes_path,annotation_path,chromosomes,just_protein_coding):
    if genes_path:
        return read_genes(genes_path, chromosomes)
    elif annotation_path:
        return read_gencode_annotation(annotation_path,just_protein_coding,chromosomes)
    else:
        return False

def read_gencode_annotation(annotation_path, just_protein_coding,chromosomes):  # just_protein_coding = bool
    genes = {}  # chrom: set[id,id]
    with open(annotation_path, 'r') as a:
        for line in a:
            if line.startswith('chr'):
                split = line.split('\t')
                if split[2] == 'gene':
                    infos = split[8].split(';')
                    for info in infos:
                        if info.startswith('gene_id'):
                            gene_id = info.split('"')[1]
                        if info.startswith(' gene_name'):
                            gene_name = info.split('"')[1]
                        if info.startswith(' gene_type'):
                            gene_type = info.split('"')[1]
                            if just_protein_coding and gene_type != 'protein_coding':
                                break
                    else:
                        chrom = split[0].split('r')[1]
                        if chrom not in chromosomes:
                            continue
                        if chrom not in genes:
                            genes[chrom] = set()
                        genes[chrom].add(gene_name + '_' + gene_id)
    return genes


def read_genes(genes_path, chromosomes):
    genes = {}  # chrom: set[id,id]
    with open(genes_path, 'r') as file:
        for line in file:
            split = line.rstrip('\n').split('\t')
            if split[4] in chromosomes:
                if split not in genes:
                    genes[split[4]] = set()
                genes[split[4]].add(split[0])
    return genes

def read_header(gfa_path):
    with open(gfa_path, 'r') as g:
        for line in g:
            header = line.rstrip('\n')
            return header


def read_index_file(index_file):  # works for links, segments, hashs
    with open(index_file, 'rb') as f:
        result = pickle.load(f)
    return result

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


def createMSA(gene_vcf_path, gene_ref_path, segments, path_seg_sets, pangenome, haplos, start_g, end_g,
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

    ref_pos = 0  ############################ was wenn variante schon am startelement? funktioniert es dann auch?
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
            #print(split)
            startseg_data = split[7].split('>')
            if len(startseg_data) > 1:
                startseg = split[7].split('>')[1]
            else:
                startseg = split[7].split('<')[1]
            samples = split[9:]

            varianten_hash_id[var_id] = split  # für nested snarls. speichere nür nötiges. Gerade: alles lol
            varianten_hash_seg[startseg] = var_id

        #print(gene_ref_path)  ##

        for n,ref_seg in enumerate(gene_ref_path_split):
            ref_seg = ref_seg[:-1]
            if ref_seg in variant:
                variant.remove(ref_seg)
                continue

            in_sample = set()
            for haplo in haplos:
                if ref_seg in path_seg_sets[haplo]:
                    in_sample.add(haplo)
            if ref_seg == '':
                continue
            seg_index = segments[ref_seg]
            g.seek(seg_index[0])
            seq = g.read(seg_index[1] - seg_index[0]).decode('utf-8')
            seq_len = len(seq)
            if strand == '-1':
                temp = start_seg
                start_seg = end_seg
                end_seg = temp
                temp = start_g
                start_g = end_g
                end_g = temp


            if ref_seg == start_seg and n == 0:
                if strand == '+':
                    seq = seq[-(seq_len - start_g):]
                else:
                    seq = seq[seq_len-start_g:]
                seq_len -= start_g
            if ref_seg == end_seg and n == len(gene_ref_path_split) - 1:
                if strand == '+':
                    seq = seq[:-(seq_len - end_g)]
                else:
                    seq = seq[:seq_len-end_g]

            if len(in_sample) == n_haplo:
                msa.append(['a', seq])
            #else:
                #b = 'b:' + give_sample_status(samples, haplo_pos)
                #msa.append([b, seq])

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
                        msa.append(['o:' + give_sample_status(samples, haplo_pos), [ref]+alt])
                if len(ref) > 1 or not allaltlen1:  # neues Alignement
                    flag = 'r:' + give_sample_status(samples, haplo_pos)
                    block = [''] * (len(alt) + 1)
                    block[0] = ref
                    for i in range(len(alt)):
                        block[i + 1] = '-' * len(ref)
                    for j in range(len(alt)):
                        for i in range(len(block)):
                            if j + 1 == i:
                                block[i] += alt[j]
                            else:
                                block[i] += '-' * len(alt[j])
                    msa.append([flag, block])
                    # msa.append([flag, call_clustalw(ref, alt, gene_vcf_path, clustalw_path)])

    ##create msa_output

    isar_path = gene_vcf_path[:-3] + 'msa'
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
            # block_seqs = [record.seq for record in block_seqs]
            block_len = len(block_seqs[0])
            for i, haplo in enumerate(haplo_order):
                if block_flag[i] == '.':
                    haplos_seqs[haplo] = haplos_seqs[haplo] + '-' * block_len
                else:
                    haplos_seqs[haplo] = haplos_seqs[haplo] + block_seqs[int(block_flag[i])]
            markers.append(['r', [aktlen, aktlen + block_len]])
            aktlen += block_len

    with open(isar_path, 'w') as isar:
        for haplo in haplo_order:
            isar.write('>' + haplo + '\n')
            # print(haplos_seqs[haplo])
            isar.write(str(haplos_seqs[haplo]) + '\n')

def give_sample_status(samples, haplo_pos):  # z.b. 0,2,1,1,1,0,.,1,0 nur für haplos
    output = ['0']
    for pos in haplo_pos:
        output.append(samples[pos])
    return ','.join(output)


def create_xg(genegfa, vg_path, strand):
    command = vg_path[:-2] + './vg'
    xg_file = genegfa[:-3] + 'xg'
    xg_arguments = ['index', genegfa, '-x', xg_file, '-t', '32', '-p']
    process = subprocess.Popen([command] + xg_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    #print(stderr.decode().strip())

    vcf_file = genegfa[:-3] + 'vcf'
    vcf_arguments = ['deconstruct', '-t', '32', '-a', '-e', '-v', '-P', 'grch38', genegfa]
    with open(vcf_file, 'w') as vcf:
        process = subprocess.Popen([command] + vcf_arguments, stdout=vcf, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    #print(stderr.decode().strip())

    fasta_file = genegfa[:-3] + 'fasta'
    fasta_arguments = ['paths', '-x', xg_file, '-F']
    process = subprocess.Popen([command] + fasta_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    # print(stdout.decode().strip())
    # print(stderr.decode().strip())
    fasta = stdout.decode().strip().split('>')[1:]
    fasta_len = len(fasta)
    fasta_dic = {}  # ATCCG: [1#2, 2#1, 26#2]
    haplo_tab = {}  # für tabellen ausgabe aller haplos hap1 = {[hap2,hap3]} samples vom selben haplotype

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

        haplo_tab[haplo_sample] = fa_ids

        record = SeqIO.SeqRecord(sequence, id=str(len(fa_ids)) + '/' + str(fasta_len) + ':' + haplo_sample,
                                 description='')  # id=str(len(fa_ids))+'/'+str(fasta_len)+':'+','.join(fa_ids)
        records.append(record)
    with open(fasta_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    haplos = []
    for rec in records:
        if not rec.id.split(':')[1].startswith('grch38'):
            haplos.append(rec.id.split(':')[1])

    # write tab file
    with open(genegfa[:-3] + 'haptab', 'w') as h:
        h.write('#example_samp\tother_samples_in_haplo\n')
        for hap, samples in haplo_tab.items():
            for i, s in enumerate(samples):
                if hap == s:
                    break
            samples.pop(i)

            h.write(hap + '\t' + ','.join(samples) + '\n')

    return [haplos, vcf_file]


def create_mini_gfa(pangenome, header, segments, id, startseg, endseg, strand, start_g, end_g, directory,
                    g_paths):
    with open(pangenome, 'rb') as g, open(directory + id + '.gfa', 'w') as o:
        segs = set()
        o.write(header + '\n')
        empty = []
        links = set()       # (10,12) für link from 10 to 12
        for g_path in g_paths:  # get all segs
            last = None
            g_path = g_path.rstrip('\n')
            for g_p_seg in g_path.split('\t')[1].split(','):
                if g_p_seg == '':
                    empty.append(g_path)
                    continue
                dir = g_p_seg[-1]
                g_p_seg = g_p_seg[:-1]
                segs.add(g_p_seg)
                if last is not None:
                    if dir == '+':
                        links.add((last,g_p_seg))
                    else:
                        links.add((g_p_seg,last))
                last = g_p_seg
        for p in empty:
            g_paths.remove(p)
        if strand == '-1':
            temp = startseg
            startseg = endseg
            endseg = temp
            temp = start_g
            start_g = end_g
            end_g = temp
        for i,seg in enumerate(segs):  # write all segs
            seg_data = segments[seg]
            g.seek(seg_data[0])
            seq = g.read(seg_data[1] - seg_data[0]).decode('utf8')
            seq_len = len(seq)
            if seg == startseg and i == 0:
                if strand == '+':
                    seq = seq[-(seq_len - start_g):]
                else:
                    seq = seq[seq_len-start_g:]
                seq_len -= start_g
            if seg == endseg and i == len(segs) - 1:
                if strand == '+':
                    seq = seq[:-(seq_len - end_g)]
                else:
                    seq = seq[:seq_len-end_g]

            o.write('S\t' + seg + '\t' + seq + '\n')
        for path in g_paths:  # write all paths
            o.write('P\t' + path + '\t*\n')
        for link in links:
            from_s = link[0]
            to_s = link[1]
            o.write('L\t' + from_s + '\t+\t' + to_s + '\t+\t*\n')

        #print(chrom + ' ' + id + '.gfa erstellt')

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


def create_mini_gfa_and_etc(pangenome, header, segments, id, startseg, endseg, starnd, start_g, end_g, directory,
                            g_paths, vg_path):
    #print('etc aufgerufen')
    minigfa_path = create_mini_gfa(pangenome, header, segments, id, startseg, endseg, starnd, start_g, end_g,
                                   directory, g_paths)
    [haplos, vcf_file] = create_xg(minigfa_path, vg_path, strand)
    [ref_path, path_set_sets] = search_ref_path_and_get_path_set_sets(g_paths)
    createMSA(vcf_file, ref_path, segments, path_set_sets, pangenome, haplos, start_g, end_g, strand)
    print(id + ' creations fertig')

def read_gene_index(file_path):
    paths = []
    with open(file_path, 'r') as g:
        for line in g:
            if line.startswith('>'):
                [id, start, end, strand, chrom, start_seg, end_seg, start_g, end_g] = line.rstrip('\n').lstrip(
                    '>').split('\t')
                start_g = int(start_g)
                end_g = int(end_g)
            else:
                if line:
                    paths.append(line.rstrip('\n'))
    return [id, strand, chrom, start_seg, end_seg, start_g, end_g, paths]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extracts genetic information out of pangenome')
    parser.add_argument('-p', '--gfa_path')
    parser.add_argument('-g', '--genes_path')
    parser.add_argument('-v', '--vg_path')
    parser.add_argument('-a', '--annotation_path')
    parser.add_argument('-n', '--chromosomes')
    parser.add_argument('-j', '--just_protein_coding', action='store_true')
    parser.add_argument('-c', '--config', help="File containing default argument values")

    args = parser.parse_args()

    if args.config:
        with open(args.config, 'r') as f:
            config_args = yaml.safe_load(f)
            parser.set_defaults(**config_args)
        args = parser.parse_args()

    gfa_path = args.gfa_path
    genes_path = False
    if args.genes_path:
        genes_path = args.genes_path  # simuliert die eingabe
    annotation_path = False
    if args.annotation_path:
        annotation_path = args.annotation_path

    numbers = [str(i + 1) for i in range(22)]
    additional_elements = ['X', 'Y', 'M']
    chrom_names = numbers + additional_elements

    chromosomes = set(chrom_names)
    if args.chromosomes:
        chromosomes = set(args.chromosomes.split(','))
    just_protein_coding = False
    if args.just_protein_coding:
        just_protein_coding = True

    vg_path = args.vg_path


    numbers = [str(i + 1) for i in range(22)]
    additional_elements = ['X', 'Y', 'M']
    chrom_names = numbers + additional_elements


    seg_index_path = 'segments/pan2gene_seg.pkl.chr'
    #link_index_path = 'links/pan2gene_link.pkl.chr'
    path_index_path = 'pan2gene_path.pkl'
    chrom_index_path = 'pan2gene_chrom.index'
    gene_index_temp = 'gene_index'

    stat_file = 'stat.info'



    #segments = read_index_file(seg_index_path)
    #print('segs eingelesen')
    #[paths, ref_paths_indexes] = read_path_ind(path_index_path, chrom_names)
    #print('paths eingelesen')
    genes = get_genes(genes_path,annotation_path,chromosomes,just_protein_coding)
    print('genes geholt')
    header = read_header(gfa_path)
    print('header eingelesen')
    #ref_paths = get_ref_paths(gfa_path, ref_paths_indexes)
    #print('ref paths geholt')
    # seg_contained_by_path = read_index_file(seg_hash_path)
    # print('ref_path_hash eingelesen')
    #[sorted_lens, ref_len_rev] = get_ref_seg_lens(ref_paths, segments)
    #print('ref seg lens berechnet')

    #gene_index_data = read_gene_index(gene_index_path)
    #print('gene_index eingelesen')
    #links = read_index_file(link_index_path)
    #print('Links eingelesen')

    chrom_dics = {}  # erstelle Chrom directories if they dont exist yet
    if not os.path.exists('out'):
        os.makedirs('out')
    for chrom in chrom_names:
        chrom_dic = 'out/chrom' + chrom
        chrom_dics[chrom] = chrom_dic
        if not os.path.exists(chrom_dic):
            os.makedirs(chrom_dic)

    for chrom in chrom_names:
        if chrom not in chromosomes:
            continue
        #lese segments und links ein
        segments = read_index_file(seg_index_path + chrom)
        #links = read_index_file(link_index_path + chrom)
        for file in os.listdir(gene_index_temp+'/chrom'+chrom):
            if genes:
                if file.split('.index')[0] not in genes[chrom]:
                    continue
            [id, strand, chrom, start_seg, end_seg, start_g, end_g, paths] = read_gene_index(os.path.join('gene_index/chrom'+chrom, file))
            create_mini_gfa_and_etc(gfa_path,header,segments,id,start_seg,end_seg,strand,start_g,end_g,chrom_dics[chrom]+'/',paths,vg_path)
            #print(id + ' fertig')
