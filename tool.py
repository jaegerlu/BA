#!/usr/bin/env python3
import sys
# sys.path.append('/home/j/jaegerl/.local/lib/python3.11/site-packages') #fixen iwann
# import gfapy
# from gfapy.sequence import rc
import bisect
import subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
#from StringIO import StringIO
from Bio import AlignIO
import os  #


def read_genes(genes_path):
    genes = []  # [id,start_bp,end_bp,strand]
    with open(genes_path, 'r') as file:
        for line in file:
            genes.append(line.split('\t'))
    return genes


def read_gfa(first_chrom_ref_seg, segments, links, paths, gfa_path, index_file, header):
    ref_path = []
    newChrom = True
    firstChrom = True

    with open(gfa_path, 'r') as g:
        if index_file:
            ind = open(index_file, 'w')
            #ind_index = 0
            #ind_indexes = [] #entählt die indexe des index files, immer das erste segment
        line = g.readline()
        tell = 0
        l = 0
        chrom_links = {}
        while line:
            if line.startswith('H'):
                header.append(line)
                if index_file:
                    ind.write(line)
                    #ind_index += len(line.encode('utf-8'))
                    #ind_indexes.append(ind_index)
            if line.startswith('S'):
                if newChrom:  ################
                    newChrom = False
                    if not firstChrom:
                        paths.append(chrom_paths)
                        if index_file:
                            for key, value in chrom_links.items():
                                value_str = ','.join(value)
                                write_line = key + '\t' + value_str + '\n'
                                #ind_index += len(write_line.encode('utf-8'))
                                ind.write(write_line)
                        links.update(chrom_links)
                        chrom_links = {}
                    chrom_paths = {}

                split = line.split('\t')
                id = split[1]
                l = len(split[2].strip('\n'))  # length
                start = tell + len(id) + 3
                end = start + l
                segments[id] = [str(start), str(end)]
                if index_file:
                    write_line = id + '\t' + str(start) + '\t' + str(end) + '\n'
                    #ind_index += len(write_line.encode('utf-8'))
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
                    print('ja')
                    ref_p = split[2].split(',')
                    ref_path = ref_path + ref_p
                    first_chrom_ref_seg.append(ref_p[0][:-1])
                    if index_file:
                        write_line = ref_p[0][:-1] + '\n'
                        #ind_index += len(write_line.encode('utf-8'))
                        #ind_indexes.append(ind_index)
                        ind.write(write_line)
                path_segs = split[2].split(',')
                smallest_s = -1
                biggest_s = 0
                first_s = True
                for p_s in path_segs:
                    p_s = int(p_s[:-1])
                    if first_s:
                        first_s = False
                        smallest_s = p_s
                    if p_s < smallest_s:
                        smallest_s = p_s
                    if p_s > biggest_s:
                        biggest_s = p_s
                chrom_paths[id] = [start, end, smallest_s, biggest_s]  #
                if index_file:
                    write_line = id + '\t' + str(start) + '\t' + str(end) + '\t' + str(smallest_s) + '\t' + str(
                            biggest_s) + '\n'
                    #ind_index += len(write_line.encode('utf-8'))
                    ind.write(write_line)

            tell = g.tell()
            line = g.readline()

        paths.append(chrom_paths)
        links.update(chrom_links)
    g.close()
    if index_file:
        for key, value in chrom_links.items():
            value_str = ','.join(value)
            write_line = key + '\t' + value_str + '\n'
            #ind_index += len(write_line.encode('utf-8'))
            ind.write(write_line)
        # ind.flush()
        ind.close()
        #return [ref_path, ind_indexes]
    #f = open(index_file, 'r')
    #print(f.read())
    #print(os.stat(index_file).st_size)
    return ref_path


def read_index_file(index_file, gfa_file, segments, links, paths, first_chrom_ref_seg, header):
    with open(index_file, 'r') as ind, open(gfa_file, 'rb') as gfa:
        chrom_path = {}
        ref_path = []
        ### ne entferne das ganze inindex und indx wieder
        for line in ind:
            split = line.rstrip('\n').split('\t')
            if split[0] == 'H':
                header.append(line)
            elif len(split) == 3:
                segments[split[0]] = [int(split[1]), int(split[2])]
                if chrom_path:
                    paths.append(chrom_path)
                    chrom_path = {}
            elif len(split) == 5:
                chrom_path[split[0]] = [int(split[1]), int(split[2]), int(split[3]), int(split[4])]
                print(split[0] + '\t' + split[3] + '\t' + split[4]) #################
                if split[0].startswith('grch38'):
                    gfa.seek(int(split[1]))
                    ref_path = ref_path + gfa.read(int(split[2]) - int(split[1])).decode('utf8').split(',')
            elif len(split) == 2:
                seg_to = split[1].split(',')
                links[split[0]] = seg_to
            elif len(split) == 1:
                first_chrom_ref_seg.append(split[0])

        paths.append(chrom_path)
    return ref_path


# print(ref_p[:50])
def get_ref_seg_lens(ref_path, ref_lens, ref_len_rev, kum_len,
                     segments, first_chrom_ref_seg,
                     chrom_start_pos):  # trägt alle daten zum ref path ein und gibt sortierte längen für bisect zurück
    # print(ref_path)
    for ref in ref_path:
        s_id = ref[:-1]
        dir = ref[-1]
        # print(s_id)
        length = int(segments[s_id][1]) - int(segments[s_id][0])
        if s_id in first_chrom_ref_seg:
            chrom_start_pos.append(kum_len)
        ref_lens[s_id] = kum_len
        ref_len_rev[kum_len] = s_id
        kum_len += length
        ref_dir[s_id] = dir
        # ref_path.append(s_id)
    '''
    if indindexes:
        f = open(gfa_path[:-3] + 'indx', 'w')
        for i in range(len(chrom_start_pos)):
            f.write(str(chrom_start_pos[i]) + '\t' + str(indindexes[i]) + '\n')
        print(len(indindexes))
        print(indindexes)
        print(len(chrom_start_pos))
        print(chrom_start_pos)
    '''

    return sorted(ref_len_rev.keys())  # sorted_lens


def createMiniGFA(gfa_path, gene, out, ref_len_rev, sorted_lens, ref_dir, segments, paths, links,
                  chrom_start_pos, header):  # return path to new gfa
    with open(gfa_path, 'rb') as g:

        # if index_file:
        # ix = open(index_file,'w')
        g_id = gene[0]
        start = int(gene[1])
        end = int(gene[2])
        # strand = int(gene[3])
        with open(out + g_id + '.gfa', 'w') as o:  # moment kann GFA.py nicht gfas schreiben? dann lieber so
            o.write(header[0])

            # print(sorted_lens)
            new_index = bisect.bisect_right(sorted_lens, start)  # finde start und end segment
            if new_index:
                start_seg = ref_len_rev[sorted_lens[new_index - 1]]
            else:
                print('Fehler bei new_index')

            new_index = bisect.bisect_right(sorted_lens, end)
            if new_index:
                end_seg = ref_len_rev[sorted_lens[new_index - 1]]
            else:
                print('Fehler bei new_index')
            # print('start_seg\t' + start_seg)
            # print('end_seg\t' + end_seg)

            start_g = start - ref_lens[start_seg] - 1  # startposition des Gens innerhalb des Startsegments
            end_g = end - ref_lens[end_seg]  # +1 evtl weg bei ende von gen excludet/includet

            # print('end_g\t' + str(end_g) + '\tend: ' + str(end) + '\tref_lens[end_seg]' + str(ref_lens[end_seg]))
            todo = [start_seg]
            links_g = set()
            completed = set()
            skip_last_seg = False
            if end_g == 0:
                skip_last_seg = True
            while todo:
                # print(todo)
                seg = todo.pop(0)
                g.seek(int(segments[seg][0]))
                seq = g.read(int(segments[seg][1]) - int(segments[seg][0])).decode('utf8')
                # print(seq)
                if seg == start_seg:
                    seq_len = len(seq)
                    if ref_dir[seg] == '+':
                        seq = seq[-(seq_len - start_g):]
                        # print(seq)
                    else:
                        seq = seq[:-(seq_len - start_g)]
                elif seg == end_seg:
                    if skip_last_seg:
                        continue
                    seq_len = len(seq)
                    if ref_dir[seg] == '+':
                        seq = seq[:-(seq_len - end_g)]
                    else:
                        seq = seq[-(seq_len - end_g):]

                # print(seg)

                o.write('S\t' + seg + '\t' + seq + '\n')

                if seg in links:
                    for n in links[seg]:
                        if int(n) <= int(end_seg) and n not in completed and n not in todo and seg != end_seg:
                            todo.append(n)

                completed.add(seg)

                # füge den link in links hinzu
                # schaue ob endseg in aufspaltung größere id hat als das andere... das andere brauchen wir nicht... obwohl dann liegt da auch ein Ende drauf

            # finde das pfade des chromosoms:
            new_index = bisect.bisect_right(chrom_start_pos, start) - 1  # finde start und end segment
            already_contained_samples = set()
            for p_id, p_data in paths[new_index].items():
                smallseg = p_data[2]
                bigseg = p_data[3]
                if int(start_seg) >= smallseg and int(end_seg) <= bigseg:  # pfad enthält gene komplett
                    p_id_split = p_id.split('#')
                    sample_id = p_id_split[0] + '#' + p_id_split[1]
                    print(p_id + '\t' + str(smallseg) + '\t' + str(bigseg))
                    if sample_id not in already_contained_samples:
                        g.seek(p_data[0])
                        p_path = g.read(p_data[1] - p_data[0]).decode('utf8').split(',')

                        out_path = ''
                        path_is_valid = False       #schaue ob der path wirklich min ein segment des Gens enthält, oder das Gen skippt(?)
                        for p_seg in p_path:
                            if skip_last_seg:
                                if int(start_seg) <= int(p_seg[:-1]) < int(end_seg):
                                    out_path = out_path + p_seg + ','
                            else:
                                if int(start_seg) <= int(p_seg[:-1]) <= int(end_seg):
                                    out_path = out_path + p_seg + ','
                            if p_seg[:-1] in completed:
                                path_is_valid = True
                        if path_is_valid:
                            out_path = out_path[:-1]
                            o.write('P\t' + sample_id + '\t' + out_path + '\t*\n')
                            already_contained_samples.add(sample_id)
                        #else: <irgendwas mit den skippenden Pfaden machen?>

            for seg_from in completed:  # schreibe links
                if seg_from == end_seg:
                    continue
                for seg_to in links[seg_from]:
                    if seg_to in completed:  #catcht größere Elemente als in gfa enthalten
                        o.write('L\t' + seg_from + '\t+\t' + seg_to + '\t+\t*\n')

        return o.name


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
                break


def create_xg(genegfa, vg_path, strand, muscle_path):
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
    #print(stdout.decode().strip())
    #print(stderr.decode().strip())
    fasta = stdout.decode().strip().split('>')[1:]
    fasta_len = len(fasta)
    fasta_dic = {}          #ATCCG: [1#2, 2#1, 26#2]
    #print(fasta)
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
        if strand == '1':
            sequence = Seq(fa_seq)
        else:
            sequence = Seq(reverse(fa_seq))

        haplo_sample = fa_ids[0]
        for id in fa_ids:
            if id.startswith('grch38'):
                haplo_sample = id
        record = SeqIO.SeqRecord(sequence, id = str(len(fa_ids))+'/'+str(fasta_len)+':'+haplo_sample, description='')  # id=str(len(fa_ids))+'/'+str(fasta_len)+':'+','.join(fa_ids)
        records.append(record)
    with open(fasta_file, 'w') as f:
        SeqIO.write(records,f,'fasta')

    msl_file = genegfa[:-3] + 'msa'
    muscle_split = muscle_path.split('/')
    if len(muscle_split) == 1:
        muscle_command = './' + muscle_path
    else:
        muscle_command = '/'.join(muscle_split[:-1]) + '/./' + muscle_split[-1]
    msl_arguments = ['-in',fasta_file,'-out',msl_file,'-clw','-diags','-maxiters','1']
    process = subprocess.Popen([muscle_command] + msl_arguments, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stderr.decode().strip())












if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extracts genetic information out of pangenome')
    parser.add_argument('-p', '--gfa_path')
    parser.add_argument('-g', '--genes_path')
    parser.add_argument('-v', '--vg_path')
    parser.add_argument('-m', '--muscle_path')
    parser.add_argument('-i', '--index_output', action='store_true')
    parser.add_argument('-j', '--index_input')
    args = parser.parse_args()

    gfa_path = args.gfa_path
    genes_path = args.genes_path  # simuliert die eingabe
    vg_path = args.vg_path
    muscle_path = args.muscle_path
    index_file = False
    if args.index_output:
        index_file = gfa_path[:-3] + 'index'  # vielleicht nicht dynamisch + optional, ob erstellt werden soll
    index_input = args.index_input
    out = 'Chromosom1/'

    genes = read_genes(genes_path)

    genes = genes[1:]  # bei ensembl Daten
    genes = genes[:3]
    print('gene Einlesung abgeschlossen :)')

    # gfa = gfapy.Gfa.from_file(gfa_path) #lese gfa-file ein
    # print('Gfa file eingelesen. Zeilen: ' + str(len(gfa.lines)))

    ref_path = []  # speichere hier alle 25 grch38 paths ein, zusammengefügt
    ref_lens = {}  # seg: kumlength
    ref_len_rev = {}  # kumlength: seg
    ref_dir = {}  # seg: +/-
    kum_len = 0
    first_chrom_ref_seg = []  # speichere erstes Segment vom ref path jedes Chromosoms
    chrom_start_pos = []  # im ref. Damit später nur für Gen ein Paths des Chromosoms gelesen werden müssen
    header = []  # nur ein element drinnen

    segments = {}  # id : [start_b, end_b]
    links = {}  # first : set(alle seconds)
    paths = []  # list für chrom: hash: id : [start_b,end_b,smallestseg,biggestseg]

    #indindexes = [] #list of all indexes from the start of a chromosome of the index-file

    if index_input:
        ref_path = read_index_file(index_input,gfa_path,segments,links,paths,first_chrom_ref_seg,header)
        print('index File Einlesung abgeschlossen')
    else:
        #if index_file:
            #[ref_path, indindexes] = read_gfa(first_chrom_ref_seg, segments, links, paths, gfa_path, index_file, header)
        #else:
        ref_path = read_gfa(first_chrom_ref_seg, segments, links, paths, gfa_path, index_file, header)
        print('gfa Einlesung abgeschlossen')

    #sorted_lens = []

    sorted_lens = get_ref_seg_lens(ref_path, ref_lens, ref_len_rev, kum_len,
                                   segments, first_chrom_ref_seg,
                                   chrom_start_pos)  # trägt alle daten zum ref path ein und gibt sortierte längen für bisect zurück
    print('Ref-Path-Zuordnung fertig')
    #f = open(index_file, 'r')
    #print(f.read())
    #print(index_file)
    #print(os.stat(index_file).st_size)

    # print(sorted_lens)
    # print(links)

    for gene in genes:
        gene_gfa = createMiniGFA(gfa_path, gene, out, ref_len_rev, sorted_lens, ref_dir, segments, paths, links,
                                 chrom_start_pos, header)
        print('Mini-GFA Creation abgeschlossen! Referenzpath:')
        printrefSeq(gene_gfa, gene[3])
        create_xg(gene_gfa, vg_path, gene[3], muscle_path)

    #print(os.stat(index_file).st_size)
