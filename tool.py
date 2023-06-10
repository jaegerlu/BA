#!/usr/bin/env python3
import sys
# sys.path.append('/home/j/jaegerl/.local/lib/python3.11/site-packages') #fixen iwann
# import gfapy
# from gfapy.sequence import rc
import bisect
import subprocess
import argparse


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

    with open(gfa_path, 'r') as g, open(index_file, 'w') as ind:
        line = g.readline()
        tell = 0
        l = 0
        while line:
            if line.startswith('H'):
                header.append(line)
            if line.startswith('S'):
                if newChrom:  ################
                    newChrom = False
                    if not firstChrom:
                        paths.append(chrom_paths)
                        ind.write(str(links))
                    chrom_paths = {}

                split = line.split('\t')
                id = split[1]
                l = len(split[2].strip('\n'))  # length
                start = tell + len(id) + 3
                end = start + l
                ind.write(id + '\t' + str(start) + '\t' + str(end) + '\n')  # + l
                segments[id] = [str(start), str(end)]

            if line.startswith('L'):
                if not newChrom:
                    newChrom = True
                split = line.split('\t')
                first = split[1]
                second = split[3]
                if first in links:
                    if second not in links[first]:
                        links[first].add(second)
                else:
                    links[first] = {second}

            if line.startswith('P'):
                split = line.split('\t')
                id = split[1]
                start = tell + len(id) + 3
                end = start + len(split[2])
                #print(id)
                if id.startswith('grch38'):
                    print('ja')
                    ref_p = split[2].split(',')
                    ref_path = ref_path + ref_p
                    first_chrom_ref_seg.append(ref_p[1])
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
                ind.write(
                    id + '\t' + str(start) + '\t' + str(end) + '\t' + str(smallest_s) + '\t' + str(biggest_s) + '\n')

            tell = g.tell()
            line = g.readline()

        paths.append(chrom_paths)
    g.close()
    return ref_path


# print(ref_p[:50])
def get_ref_seg_lens(ref_path, ref_lens, ref_len_rev, kum_len,
                     segments, first_chrom_ref_seg,
                     chrom_start_pos):  # trägt alle daten zum ref path ein und gibt sortierte längen für bisect zurück
    print(ref_path)
    for ref in ref_path:
        s_id = ref[:-1]
        dir = ref[-1]
        #print(s_id)
        length = int(segments[s_id][1]) - int(segments[s_id][0])
        if s_id in first_chrom_ref_seg:
            chrom_start_pos.append(kum_len)
        ref_lens[s_id] = kum_len
        ref_len_rev[kum_len] = s_id
        kum_len += length
        ref_dir[s_id] = dir
        # ref_path.append(s_id)

    return sorted(ref_len_rev.keys())  # sorted_lens


def createMiniGFA(gfa_path, gene, out, ref_len_rev, sorted_lens, ref_dir, segments, paths, links,
                  chrom_start_pos, header):  # return path to new gfa
    with open(gfa_path, 'rb') as g:

        if index_file:
            ix = open(index_file,'w')
        g_id = gene[0]
        start = int(gene[1])
        end = int(gene[2])
        # strand = int(gene[3])
        with open(out + g_id + '.gfa', 'w') as o:  # moment kann GFA.py nicht gfas schreiben? dann lieber so
            o.write(header[0])
            
            #print(sorted_lens)
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
            #print('start_seg\t' + start_seg)
            #print('end_seg\t' + end_seg)

            start_g = start - ref_lens[start_seg] - 1  # startposition des Gens innerhalb des Startsegments
            end_g = end - ref_lens[end_seg]  # +1 evtl weg bei ende von gen excludet/includet

            #print('end_g\t' + str(end_g) + '\tend: ' + str(end) + '\tref_lens[end_seg]' + str(ref_lens[end_seg]))
            todo = [start_seg]
            links_g = set()
            completed = set()
            while todo:
                seg = todo.pop(0)
                g.seek(int(segments[seg][0]))
                seq = g.read(int(segments[seg][1]) - int(segments[seg][0])).decode('utf8')
                #print(seq)
                if seg == start_seg:
                    seq_len = len(seq)
                    if ref_dir[seg] == '+':
                        seq = seq[-(seq_len - start_g):]
                        #print(seq)
                    else:
                        seq = seq[:-(seq_len - start_g)]
                elif seg == end_seg:
                    seq_len = len(seq)
                    if ref_dir[seg] == '+':
                        seq = seq[:-(seq_len - end_g)]
                    else:
                        seq = seq[-(seq_len - end_g):]

                #print(seg)
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
                    if sample_id not in already_contained_samples:
                        g.seek(p_data[0])
                        p_path = g.read(p_data[1] - p_data[0]).decode('utf8').split(',')
                        out_path = ''
                        for p_seg in p_path:
                            if int(start_seg) <= int(p_seg[:-1]) <= int(end_seg):
                                out_path = out_path + p_seg + ','
                        out_path = out_path[:-1]
                        o.write('P\t' + sample_id + '\t' + out_path + '\t*\n')
                        already_contained_samples.add(sample_id)

            for seg_from in completed:  # schreibe links
                if seg_from == end_seg:
                    continue
                for seg_to in links[seg_from]:
                    o.write('L\t' + seg_from + '\t+\t' + seg_to + '\t+\t*\n')

        return o.name


def reverse(seq):
    revcomp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    revcomp = ""
    for nt in reversed(seq):
        if nt == 'N':
            revcomp += 'N'
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

def create_xg(genegfa, vg_path):
    command = vg_path[:-2] + './vg'
    xg_file = genegfa[:-3] + 'xg'
    xg_arguments = ['index', genegfa, '-x', xg_file, '-t', '32', '-p']
    process = subprocess.Popen([command] + xg_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stderr.decode().strip())

    vcf_file = genegfa[:-3] + 'vcf'
    vcf_arguments = ['deconstruct', '-t', '32', '-a', '-e', '-P', 'grch38', genegfa, '>', vcf_file]
    process = subprocess.Popen([command] + vcf_arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stderr.decode().strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extracts genetic information out of pangenome')
    parser.add_argument('-p', '--gfa_path')
    parser.add_argument('-g', '--genes_path')
    parser.add_argument('-v', '--vg_path')
    parser.add_argument('-i', '--index_output', action='store_true')
    parser.add_argument('-j', '--index_input')
    args = parser.parse_args()

    gfa_path = args.gfa_path
    genes_path = args.genes_path  # simuliert die eingabe
    vg_path = args.vg_path
    index_file = False
    if args.index_output:
        index_file = gfa_path[:3] + 'index'  # vielleicht nicht dynamisch + optional, ob erstellt werden soll
    index_input = args.index_input
    out = 'Chromosom1/'

    genes = read_genes(genes_path)

    genes = genes[1:] # bei ensembl Daten
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
    header = [] # nur ein element drinnen

    segments = {}  # id : [start_b, end_b]
    links = {}  # first : set(alle seconds)
    paths = []  # list für chrom: hash: id : [start_b,end_b,smallestseg,biggestseg]

    ref_path = read_gfa(first_chrom_ref_seg, segments, links, paths, gfa_path, index_file, header)
    print('gfa Einlesung abgeschlossen')

    sorted_lens = []

    sorted_lens = get_ref_seg_lens(ref_path, ref_lens, ref_len_rev, kum_len,
                                   segments, first_chrom_ref_seg,
                                   chrom_start_pos)  # trägt alle daten zum ref path ein und gibt sortierte längen für bisect zurück
    print('Ref-Path-Zuordnung fertig')
    print(sorted_lens)

    for gene in genes:
        gene_gfa = createMiniGFA(gfa_path, gene, out, ref_len_rev, sorted_lens, ref_dir, segments, paths, links,
                                 chrom_start_pos, header)
        print('Mini-GFA Creation abgeschlossen! Referenzpath:')
        printrefSeq(gene_gfa,gene[3])
        create_xg(gene_gfa, vg_path)

