#!/usr/bin/env python3
from multiprocessing import set_start_method
try:
    set_start_method("spawn")
except RuntimeError:
    pass
from multiprocessing import Manager, Pool
import os
import argparse
import pickle
import bisect
import portion as P
import gc
import random


def read_genes(genes_path):
    genes = []  # [id,start_bp,end_bp,strand,chrom]
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
                            gene_id = info.split('"')[1]
                        if info.startswith(' gene_name'):
                            gene_name = info.split('"')[1]
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
                        genes.append([gene_name + '_' + gene_id, start, end, strand, chrom])
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


def read_path_ind(path_ind_file, chrom_names):
    with open(path_ind_file, 'rb') as f:
        paths = pickle.load(f)

    ref_paths_indexes = {}
    for chrom in chrom_names:
        ref_paths_indexes[chrom] = paths[chrom]['grch38#chr' + chrom]  # oder?
    return [paths, ref_paths_indexes]


def get_ref_paths(gfa_file, ref_paths_indexes):
    ref_paths = {}
    with open(gfa_file, 'rb') as g:
        for chrom, pos in ref_paths_indexes.items():
            g.seek(pos[0])
            path = g.read(pos[1] - pos[0]).decode('utf8').split(',')  # schon gesplittet
            ref_paths[chrom] = path
    return ref_paths


def reverse_path(path, is_split=False):
    out = ''
    dirs = {'+': '-', '-': '+'}
    if is_split:
        path_split = path
    else:
        path_split = path.split(',')
    for p in path_split:
        # if len(p) < 2:
        # print(path)
        # print(p)
        seg = p[:-1]
        dir = p[-1]
        out = seg + dirs[dir] + ',' + out
    return out[:-1]


def create_gene_index_file(gfa_file, paths, genes, sorted_lens, ref_len_rev, gene_index_path, ref_path,
                           promoter_length, n_threads, gene_index_temp, chromosomes, stats, max_val_size):
    complete_genes = {}  # chrom:  list
    # segs_need_to_be_hashed = {} #chrom : list[123,53]
    for gene in genes:
        [id, start, end, strand, chrom] = gene

        if int(start) > int(end):  # switche, bei opposite strand manchmal end kleiner als start
            temp = start
            start = end
            end = temp
            print(id + 'start > end ' + ('$%&' * 50))

        if strand == '1':
            start = int(start) - promoter_length
            if start < 0:
                start = 0
        else:
            end = int(end) + promoter_length

        new_indexs = bisect.bisect_right(sorted_lens[chrom], int(start)) - 1  # finde start und end segment
        if isinstance(new_indexs, int):
            start_seg_pos = sorted_lens[chrom][new_indexs]
            start_seg = ref_len_rev[chrom][start_seg_pos][0]
            start_ref_index = ref_len_rev[chrom][start_seg_pos][1]

        else:
            print('Fehler bei new_index')
            print(new_indexs)
            print(bisect.bisect_right(sorted_lens[chrom], int(start)))
            print(start)

        new_indexe = bisect.bisect_right(sorted_lens[chrom], int(end)) - 1
        if isinstance(new_indexe, int):
            end_seg_pos = sorted_lens[chrom][new_indexe]
            end_seg = ref_len_rev[chrom][end_seg_pos][0]
            end_ref_index = ref_len_rev[chrom][end_seg_pos][1]
        else:
            print('Fehler bei new_index')
            print(new_indexe)
            print(bisect.bisect_right(sorted_lens[chrom], int(end)))
            print(end)

        start_g = int(
            start) - start_seg_pos - 1  # startposition des Gens innerhalb des Startsegments #### -1 wsh weg? #bzw falls unter 0, ein segment navh vorne?
        if start_g == -1:
            new_indexs -= 1
            start_seg_pos = sorted_lens[chrom][new_indexs]
            start_seg = ref_len_rev[chrom][start_seg_pos][0]
            start_ref_index = ref_len_rev[chrom][start_seg_pos][1]
            start_g = int(start) - start_seg_pos - 1

        end_g = int(end) - end_seg_pos  # +1 evtl weg bei ende von gen excludet/includet
        if end_g == 0:
            new_indexe -= 1
            end_seg_pos = sorted_lens[chrom][new_indexe]
            end_seg = ref_len_rev[chrom][end_seg_pos][0]
            end_ref_index = ref_len_rev[chrom][end_seg_pos][1]
            end_g = int(end) - end_seg_pos

        if chrom not in complete_genes:
            complete_genes[chrom] = []
            # links = read_index_file(links_path + chrom)
            # links_rev = read_index_file(links_rev_path + chrom)
            # print(chrom)
            # print_RAM()

        complete_genes[chrom].append(
            [id, start, end, strand, chrom, start_seg, end_seg, start_g, end_g,
             start_ref_index, end_ref_index])
   # with open('complete_genes', 'wb') as ind:
    #    pickle.dump(complete_genes, ind)

    print('gene index step 1 fertig')
    links = 0
    links_rev = 0
    sorted_lens = 0
    ref_len_rev = 0
    genes = 0
    gc.collect()

    manager = Manager()
    lock = manager.Lock()
    if os.path.exists(gene_index_path):
        # Lösche die Datei
        os.remove(gene_index_path)

    if not os.path.exists(gene_index_temp):
        os.makedirs(gene_index_temp)

    for chrom, genes in complete_genes.items():
        if chrom not in chromosomes:
            continue
        print(chrom)
        print('n genes in chrom: ' + str(len(genes)))

        ref_sub_paths = {}

        if stats:
            stats = 'pan2index_stats_chr' + chrom + '.txt'
            if os.path.exists(stats):
                os.remove(stats)

        for gene in genes:
            [id, sstart, send, strand, chrom, start_seg, end_seg, start_g, end_g,
             start_ref_index, end_ref_index] = gene

            if not os.path.exists(gene_index_temp + '/chrom' + chrom):
                os.makedirs(gene_index_temp + '/chrom' + chrom)

            with open(gene_index_temp + '/chrom' + chrom + '/' + id + '.index',
                      'w') as f:  # evtl lass file mit start starten, damit es beim zusammenfügen dann auch in der richtigen reihenfolge
                f.write('>' + '\t'.join(str(x) for x in gene[:-2]))

            ref_sub_paths[id] = ref_path[chrom][start_ref_index: end_ref_index + 1]

        print('Segment Indexing Chromosome '+chrom+'...')
        pos = create_segment_path_index(paths[chrom], gfa_file, ref_sub_paths)
        if stats:
            pos_stats(pos, stats)

        with Pool(processes=n_threads) as pool:
            for i, gene in enumerate(genes):
                print(i)
                [gene_id, sstart, send, strand, chrom, start_seg, end_seg, start_g, end_g,
                 start_ref_index, end_ref_index] = gene
                bucket = generate_bucket(pos, ref_sub_paths[gene_id])
                outstring = str(i) + '/' + str(len(genes))
                pool.apply_async(search_for_gene_via_indeces,
                                 (lock, bucket, ref_sub_paths[gene_id], paths[chrom], gene_id, chrom, strand,
                                  gene_index_temp, gfa_file, max_val_size, outstring, stats, True,))

            pool.close()
            pool.join()
        pos = 0  # damit nicht zwei pos gleichzeitig geladen werden müssen


def pos_stats(pos, stats):
    max_segs = {}  # id:anzahl
    most_appear = []
    out = []
    for seg, appear in pos.items():
        segtemp = {}  # pathid: häufigkeit
        whole = len(appear)
        for tupel in appear:
            path = tupel[0]
            if path not in segtemp:
                segtemp[path] = 1
            else:
                segtemp[path] = segtemp[path] + 1
        if len(most_appear) < 3:
            most_appear.append(whole)
            max_segs[seg] = whole
        elif whole > min(most_appear):
            m = min(most_appear)
            most_appear.remove(m)
            most_appear.append(whole)
            for s, n in max_segs.items():
                if n == m:
                    break
            max_segs.pop(s)
            max_segs[seg] = whole
        all = sum(segtemp.values())
        out.append([len(segtemp.values()), all])
    with open(stats, 'a') as s:
        s.write('best\t' + str(max_segs) + '\n')
        print(max_segs)
    with open(stats + 'pos.pkl', 'wb') as ind:
        pickle.dump(out, ind)
    #######################################


def search_for_gene_via_indeces(lock, bucket, ref_sub_path, paths_chrom, gene_id, chrom, strand,
                                gene_index_temp,
                                gfa,
                                max_val_size, outstring, stats,
                                is_int=500):  # is_int True für ob man weiß ob die Segment ids im pangenom ints:

    results = {}
    print('bucket befüllt######################################################')
    ref_set = set()
    ref_set_rev = set()
    dirs = {'+': '-', '-': '+'}

    span = is_int
    if span == True:
        #print('f1')
        last_seg = int(ref_sub_path[0][:-1])
        span = 0
        for seg in ref_sub_path:
            #print('f2')
            ref_set.add(seg)
            #print('f3')
            dir = seg[-1]
            seg = seg[:-1]
            #print('f4')
            ref_set_rev.add(seg + dirs[dir])
            #print(seg)
            seg = int(seg)
            #print('f5')
            dif = abs(seg - last_seg)
            if dif > span:
                span = dif
            #print('f6')
            last_seg = seg
            #print('f7')
        if span < 10:
            span = 10
        if span > 500:
            span = 500
        #print('f8')

    else:
        for seg in ref_sub_path:
            ref_set.add(seg)
            dir = seg[-1]
            seg = seg[:-1]
            ref_set_rev.add(seg + dirs[dir])
    #print('ref sets erstellt')
    #print('flag 1')
    for sample_id, dir1 in bucket.items():
        bestscore = -float('inf')
        bestval = None  # [start,end]                                            # bzw [[],[]] falls contigübergrenzend
        # safe_other_tupels = None
        # results[sample_id] = None  # hier am Ende sample_id: [bestinterval]
        # print(sample_id)
        for contig_id, dir2 in dir1.items():
            # print(contig_id)
            if contig_id != sample_id:
                [cstart, cend] = paths_chrom[sample_id + '#' + contig_id]
            else:
                [cstart, cend] = paths_chrom[sample_id]
            conlen = cend - cstart
            val = []

            for i in dir2:
                val.append(check_tupel_borders([i - span, i + span], conlen))
            if len(val) > max_val_size:
                val = reduce_val(val, max_val_size)
            # print(len(val))
            merge = merge_overlapping_onelist(val)  # [[+/-,count,start,end], ...]
            if stats:
                with lock:
                    with open(stats, 'a') as stat:
                        stat.write(
                            'v+\t' + str(len(val)) + '\t' + str(len(merge)) + '\t' + str(len(ref_sub_path)) + '\n')
            # better = False
            for m in merge:
                [m_start, m_end] = [m[1], m[2]]  # check_tupel_borders([m[1],m[2]],conlen)
                size = m_end - m_start

                count = m[0]
                score = count / size

                if score > bestscore:
                    # better = True
                    bestscore = score
                    bestval = [contig_id, m_start, m_end]

        # jetzt irgendwie bestes aln richtig zuschneiden und überprüfen
        #print('flag 2')
        samp_result = ''
        if bestval is not None:  # and len(bestval) != 2:
            idk_flag = False
            [val_contig, val_start, val_end] = bestval
            # falls start > end:
            if val_start > val_end:
                temp = val_start
                val_start = val_end
                val_end = temp

            start_segs = []
            end_segs = []

            if sample_id == val_contig:
                [ps, pe] = paths_chrom[sample_id]
            else:
                [ps, pe] = paths_chrom[sample_id + '#' + val_contig]
            with lock:
                with open(gfa, 'rb') as g:  ########### with lock
                    g.seek(ps)
                    path = g.read(pe - ps).decode('utf8').split(',')
                if stats:
                    with open(stats, 'a') as stat:
                        stat.write('s\t' + str(bestscore) + '\n')
            subpath = path[val_start:val_end + 1]

            crosses_contig_i = False
            crosses_contig_j = False
            if val_start == 0:
                crosses_contig_i = True
            if val_end >= len(path) - 1:
                crosses_contig_j = True

            pos_count = 0
            neg_count = 0
            for s in subpath:  # check if contig is reversed
                if s in ref_set:
                    pos_count += 1
                if s in ref_set_rev:
                    neg_count += 1

            if pos_count >= neg_count:
                val_dir = '+'
                dir_set = ref_set
            else:
                val_dir = '-'
                dir_set = ref_set_rev

            i_ = None
            j_ = None

            if len(ref_sub_path) >= 6:
                start_segs = ref_sub_path[:3]
                end_segs = ref_sub_path[-3:]
            # elif len(ref_sub_path) > 3:
            #    start_segs = ref_sub_path[:2]
            #    end_segs = ref_sub_path[-2:]
            else:
                start_segs = [ref_sub_path[0]]
                end_segs = [ref_sub_path[-1]]
            if val_dir == '-':  # and strand == '1' or val_dir == '+' and strand == '-1':
                temp = start_segs
                start_segs = reverse_path(end_segs, True).split(',')
                end_segs = reverse_path(temp, True).split(',')

            if len(start_segs[0]) == 1:
                print(ref_sub_path)

            for i, seg in enumerate(subpath):
                # print(seg+'\t'+start_segs[0])
                if seg == start_segs[0]:
                    break
                if len(start_segs) > 1:
                    if seg == start_segs[1]:
                        i -= 1
                        break
                    elif seg == start_segs[2]:
                        i -= 2
                        break
                if i >= 2 * span:  # gehe so lange weiter bis du auf ein element triffst das teil ist.
                    # i = ?
                    i = -1
                    if crosses_contig_i:
                        i = 0
                    break

            if i == -1:  # gehe ducrh und nehme erstes passendes segment falls startsegment nicht gefunden
                i_, flag = enlarge_subpath(path, '-', start_segs, dir_set, val_start, len(ref_sub_path))
                if i_ is not None:
                    subpath = path[i_:val_end + 1]
                    # i_ = 0
                    if flag == True:
                        crosses_contig_i = True
                else:
                    i = span
                    idk_flag = True

            # jetzt von hinten
            # print('------------------------')
            for j in range(len(subpath) - 1, 0, -1):
                # print(seg+'\t'+end_segs[-1])
                seg = subpath[j]
                if seg == end_segs[-1]:
                    break
                if len(end_segs) > 1:
                    if seg == end_segs[-2]:
                        j += 1
                        break
                    elif seg == end_segs[-3]:
                        j += 2
                        break
                if j < len(subpath) - 2 * span:
                    j = -1
                    if crosses_contig_j:  # unterscheide zwischen i und j
                        j = len(subpath) - 1
                    # j = ?
                    break

            if j == -1:  # gehe ducrh und nehme erstes passendes segment falls startsegment nicht gefunden
                j_, flag = enlarge_subpath(path, '+', end_segs, dir_set, val_end, len(ref_sub_path))
                if j_ is not None:
                    if i_ is None:
                        subpath = path[val_start + i:j_ + 1]
                    else:
                        subpath = path[i_:j_ + 1]
                    if flag == True:
                        crosses_contig_j = True
                else:
                    j = len(subpath) - 1 - span
                    idk_flag = True
            # print('i: ' + str(i))
            # print('len - j: ' + str(len(subpath) - 1 - j))
            safej = str(len(subpath) - 1 - j)
            if i_ is not None and j_ is not None:
                res_path = subpath
            elif i_ is not None:
                res_path = subpath[:j + 1]
            elif j_ is not None:
                res_path = subpath
            else:
                res_path = subpath[i:j + 1]
            if strand == '-1' and val_dir == '+' or strand == '1' and val_dir == '-':
                samp_result += reverse_path(res_path, True)
            else:
                samp_result += ','.join(res_path)

            samp_result = '\t' + samp_result  # str(bestval)+'span:'+str(span)+' val_dir:'+val_dir+' i:'+str(i)+' len-j:'+safej+' i_='+str(i_)+' j_:'+str(j_) + '\t' +

            if crosses_contig_i:
                samp_result = '*' + samp_result
            if crosses_contig_j:
                samp_result = '*' + samp_result
            if idk_flag:
                samp_result = '?' + samp_result

            # fall für contiggrenzen überschreitend

            results[sample_id] = samp_result
    #print('flag 3')
    with open(gene_index_temp + '/chrom' + chrom + '/' + gene_id + '.index', 'a') as f:
        #print('flag 4')
        for samp, path in results.items():
            f.write('\n' + samp + path)
            # print('\n' + samp + '\t' + path)
        if strand == '-1':
            ref = reverse_path(ref_sub_path, True)
        else:
            ref = ','.join(ref_sub_path)
        f.write('\n' + 'grch38#chr' + chrom + '\t' + ref)

    print(gene_id + ' fertig: ' + outstring)


def create_segment_path_index(paths_chrom, gfa, ref_sub_paths):
    seg_positions = {}  # pro chrom: seg-id:[[path,index],]
    refset = set()
    # refdirset = set()
    for ref_path in ref_sub_paths.values():
        for seg in ref_path:
            # refdirset.add(seg)
            refset.add(seg[:-1])

    with open(gfa, 'rb') as g:
        for pid, pval in paths_chrom.items():
            #print(pid)
            if pid.startswith('grch38'):
                continue
            g.seek(pval[0])
            path = g.read(pval[1] - pval[0]).decode('utf8')
            for i, seg in enumerate(path.split(',')):
                dir = seg[-1]
                seg_ = seg[:-1]
                if seg_ in refset:
                    if seg_ not in seg_positions:
                        seg_positions[seg_] = []  # macht es Sinn zu unterscheiden?

                    seg_positions[seg_].append([pid, i])
    return seg_positions


def generate_bucket(pos, ref_sub_path):
    bucket = {}  # sample_id: contig_id:  [i1,i2,i3]
    # print(len(ref_sub_path))
    for seg in ref_sub_path:  # befüllung von bucket
        seg = seg[:-1]  # irgendwie muss dann noch überprüft werden ob es entgegengesetzt zum ref path ist

        if seg not in pos:
            continue
        if len(pos[seg]) > 1000:  # zu häufig auftretende segmente sind nicht hilfreich
            continue

        for po in pos[seg]:
            [p_id, index] = po
            id_split = p_id.split('#')
            sample_id = '#'.join(id_split[:2])
            if len(id_split) > 2:
                contig_id = '#'.join(id_split[2:])
            else:
                contig_id = sample_id  # für grch38 und c..m13
            if sample_id not in bucket:
                bucket[sample_id] = {}
            if contig_id not in bucket[sample_id]:
                bucket[sample_id][contig_id] = []
            bucket[sample_id][contig_id].append(index)
    return bucket


def reduce_val(val, max_val_size):
    while len(val) > max_val_size:
        index_to_remove = random.randint(0, len(val) - 1)
        val.pop(index_to_remove)
    return val


def enlarge_subpath(path, dir, border_segs, dir_set, startposition, gene_len):
    # print('enlarge called')
    if dir == '+':  # hier heißt dir nur, ob nach links oder rechts geschaut werden muss, also hier das Ende offen
        j = startposition
        safej = None
        while j <= len(path) - 1:
            if path[j] in border_segs:
                safej = j
                for p in range(1, 6):
                    if p + j > len(path) - 1:
                        break
                    if path[p + j] in border_segs:
                        safej = p + j
                break  # ende gefunden
            j += 1

        if safej is not None:
            return safej, None
        else:
            j -= 1
            look = int(gene_len * 0.1)
            if look > j - startposition:
                look = j - startposition
            count = 0
            for seg in path[-look:]:
                if seg in dir_set:
                    count += 1
            score = count / look
            # print(j)
            if len(path) > 2:
                if score >= 0.2 and (path[j] in dir_set or path[j - 1] in dir_set or path[j - 2] in dir_set):
                    return j, True
                else:
                    return None, False
            else:
                if score >= 0.2:
                    return j, True
                else:
                    return None, False

    if dir == '-':  # nach links gehen
        i = startposition
        safei = None
        while i >= 0:
            if path[i] in border_segs:
                safei = i
                for p in range(1, 6):
                    if i - p < 0:
                        break
                    if path[i - p] in border_segs:
                        safei = i - p
                break  # ende gefunden
            i -= 1

        if safei is not None:
            return safei, None
        else:
            i += 1
            look = int(gene_len * 0.1)
            if look > startposition - i:
                look = startposition - i
            count = 0
            for seg in path[:look]:
                if seg in dir_set:
                    count += 1
            score = count / look
            # print(i)
            if len(path) > 2:
                if score >= 0.2 and (path[i] in dir_set or path[i + 1] in dir_set or path[i + 2] in dir_set):
                    return i, True
                else:
                    return None, False
            else:
                if score >= 0.2:
                    return i, True
                else:
                    return None, False


def check_tupel_borders(tupel, length):
    if tupel[0] < 0:
        tupel[0] = 0
    if tupel[1] > length - 1:
        tupel[1] = length - 1
    return tupel


def merge_overlapping_onelist(interval_list):
    intervals = [P.closed(a, b) for a, b in interval_list]
    merge = P.Interval(*intervals)
    out = []
    for interval in merge:
        # Zähle, wie oft jedes Intervall in der Liste vorkommt
        count = sum(1 for i in interval_list if i[0] >= interval.lower and i[1] <= interval.upper)
        out.append([count, interval.lower, interval.upper])
    return out  # [[count,start,end],[c,s,e]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extracts genetic information out of pangenome')
    parser.add_argument('-p', '--gfa_path')
    parser.add_argument('-g', '--genes_path')
    parser.add_argument('-a', '--annotation_path')
    parser.add_argument('-r', '--promoter_length')
    parser.add_argument('-t', '--n_processes')
    parser.add_argument('-n', '--chromosomes')
    parser.add_argument('-j', '--just_protein_coding', action='store_true')
    parser.add_argument('-s', '--stats', action='store_true')
    parser.add_argument('-m', '--max_val_size')
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
    stats = False
    if args.stats:
        stats = True
    max_val_size = 10000
    if args.max_val_size:
        max_val_size = int(args.max_val_size)

    n_threads = 1
    if args.n_processes:
        n_threads = int(args.n_processes)

    promoter_length = int(args.promoter_length)

    numbers = [str(i + 1) for i in range(22)]
    additional_elements = ['X', 'Y', 'M']
    chrom_names = numbers + additional_elements

    chromosomes = chrom_names
    if args.chromosomes:
        chromosomes = args.chromosomes.split(',')
    just_protein_coding = False
    if args.just_protein_coding:
        just_protein_coding = True

    gene_index_path = 'pan2gene_genes.index'
    seg_index_path = 'segments/pan2gene_seg.pkl.chr'
    link_index_path = 'links/pan2gene_link.pkl.chr'
    # link_rev_index_path = 'links/pan2gene_link_rev.pkl.chr'
    path_index_path = 'pan2gene_path.pkl'
    chrom_index_path = 'pan2gene_chrom.index'
    gene_type_path = 'pan2gene_gene_type_per_gene'
    ref_lens_path = 'pan2gene_ref_lens.pkl'
    sorted_lens_path = 'pan2gene_sorted_lens.pkl'
    gene_index_temp = 'gene_index'
    # seg_hash_path = 'pan2gene_refpathseg_path_hash.pkl'
    stat_file = 'stat.info'

    if genes_path:
        genes = read_genes(genes_path)
    elif annotation_path:
        genes = read_gencode_annotation(annotation_path, just_protein_coding, gene_type_path)
    print('genes eingelesen')

    # segments = read_index_file(seg_index_path)
    # print('segs eingelesen')
    [paths, ref_paths_indexes] = read_path_ind(path_index_path, chrom_names)
    print('paths eingelesen')
    header = read_header(gfa_path)
    print('Header eingelesen')

    ref_paths = get_ref_paths(gfa_path, ref_paths_indexes)
    print('ref paths geholt')
    # [sorted_lens, ref_len_rev] = get_ref_seg_lens(ref_paths, segments)
    # print('ref seg lens berechnet')
    # segments = 0
    # with open(ref_lens_path, 'wb') as ind:
    #    pickle.dump(ref_len_rev, ind)
    # with open(sorted_lens_path, 'wb') as ind:
    #    pickle.dump(sorted_lens, ind)

    ref_len_rev = read_index_file(ref_lens_path)
    print('ref lens rev gelesen')
    sorted_lens = read_index_file(sorted_lens_path)
    print('sorted_lens gelesen')

    create_gene_index_file(gfa_path, paths, genes, sorted_lens, ref_len_rev, gene_index_path,
                           ref_paths,
                           promoter_length, n_threads, gene_index_temp, chromosomes, stats, max_val_size)
    print('gene_index_files_created')
