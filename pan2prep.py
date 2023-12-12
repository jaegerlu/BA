#!/usr/bin/env python3
import pickle
import argparse
import os

def read_gfa(gfa_path, chrom_names, seg_index_path, path_index_path, chrom_index_path, link_index_path):
    header = ''
    newChrom = True
    firstChrom = True
    segments = {}
    links = {}
    paths = {}
    ref_paths_indexes = {}  # chrom : [start,end]
    akt_chrom = '0'

    with open(gfa_path, 'r') as g, open(chrom_index_path, 'w') as chr_ind:

        line = g.readline()
        tell = 0

        if not os.path.exists(link_index_path.split('/')[0]):
            os.makedirs(link_index_path.split('/')[0])
        if not os.path.exists(seg_index_path.split('/')[0]):
            os.makedirs(seg_index_path.split('/')[0])

        while line:
            if line.startswith('H'):
                header = line.rstrip('\n')
                # ind.write(line)
            if line.startswith('S'):
                if newChrom:  ################
                    chr_ind.write(akt_chrom + '\t' + str(tell))
                    newChrom = False
                    if not firstChrom:
                        paths[akt_chrom] = chrom_paths
                        with open(seg_index_path + akt_chrom, 'wb') as ind:
                            pickle.dump(segments, ind)

                        with open(link_index_path + akt_chrom, 'wb') as ind:
                            pickle.dump(links, ind)
                        links = {}
                        segments = {}
                    akt_chrom = chrom_names.pop(0)

                    chrom_paths = {}

                split = line.split('\t')
                id = split[1]
                l = len(split[2].strip('\n'))  # length
                start = tell + len(id) + 3
                end = start + l
                segments[id] = [start, end]

                # write_line = id + '\t' + str(start) + '\t' + str(end) + '\n'
                # ind.write(write_line)  # + l

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
                # print(chrom_links)

            if line.startswith('P'):
                split = line.split('\t')
                id = split[1]
                start = tell + len(id) + 3
                end = start + len(split[2])
                # print(id)
                if id.startswith('grch38'):
                    ref_paths_indexes[akt_chrom] = [start, end]
                    firstChrom = False

                chrom_paths[id] = [start, end]  #

                # write_line = id + '\t' + str(start) + '\t' + str(end) + '\t' + akt_chrom + '\n'
                # path_ind.write(write_line)

            tell = g.tell()
            line = g.readline()

        paths[akt_chrom] = chrom_paths

        # for s, i in segments.items():
        #    if s in links:
        #        value_str = ','.join(links[s])
        #    else:
        #        value_str = ''
        #
        #    write_line = s + '\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + value_str + '\n'
        #    ind.write(write_line)

        '''
        for key, value in chrom_links.items():
            value_str = ','.join(value)
            write_line = key + '\t' + str(segments[key][0]) + '\t' + str(segments[key][1]) + '\t' + value_str + '\n'
            ind.write(write_line)

        for key, value in chrom_links.items():
            value_str = ','.join(value)
            write_line = key + '\t' + value_str + '\n'
            ind.write(write_line)
        '''

        chr_ind.write(akt_chrom + '\t' + str(tell))
        with open(seg_index_path + akt_chrom, 'wb') as ind:
            pickle.dump(segments, ind)
        with open(link_index_path + akt_chrom, 'wb') as ind:
            pickle.dump(links, ind)
        with open(path_index_path, 'wb') as ind:
            pickle.dump(paths, ind)

    return [paths, ref_paths_indexes]


def get_ref_paths(gfa_file, ref_paths_indexes):
    ref_paths = {}
    with open(gfa_file, 'rb') as g:
        for chrom, pos in ref_paths_indexes.items():
            g.seek(pos[0])
            path = g.read(pos[1] - pos[0]).decode('utf8').split(',')  # schon gesplittet
            ref_paths[chrom] = path
    return ref_paths

def get_ref_seg_lens(ref_path,
                     seg_index_path):  # trägt alle daten zum ref path ein und gibt sortierte längen für bisect zurück
    ref_len_rev = {}
    sorted_lens = {}
    for chrom, path in ref_path.items():
        segments = read_index_file(seg_index_path+chrom)
        kum_len = 0
        ref_len_rev_chrom = {}
        sorted_lens_chrom = []
        print('ref_lens '+chrom)
        for i, ref in enumerate(path):
            s_id = ref[:-1]
            length = segments[s_id][1] - segments[s_id][0]

            ref_len_rev_chrom[kum_len] = [s_id, i]

            sorted_lens_chrom.append(kum_len)

            kum_len += length
        ref_len_rev[chrom] = ref_len_rev_chrom
        sorted_lens[chrom] = sorted_lens_chrom


    return [sorted_lens, ref_len_rev]

def read_index_file(index_file):  # works for links, segments, hashs
    with open(index_file, 'rb') as f:
        result = pickle.load(f)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extracts genetic information out of pangenome')
    parser.add_argument('-p', '--gfa_path')
    parser.add_argument('-c', '--config', help="File containing default argument values")
    args = parser.parse_args()

    if args.config:
        with open(args.config, 'r') as f:
            config_args = yaml.safe_load(f)
            parser.set_defaults(**config_args)
        args = parser.parse_args()

    gfa_path = args.gfa_path
    numbers = [str(i + 1) for i in range(22)]
    additional_elements = ['X', 'Y', 'M']
    chrom_names = numbers + additional_elements




    seg_index_path = 'segments/pan2gene_seg.pkl.chr'
    link_index_path = 'links/pan2gene_link.pkl.chr'
    path_index_path = 'pan2gene_path.pkl'
    chrom_index_path = 'pan2gene_chrom.index'
    ref_lens_path = 'pan2gene_ref_lens.pkl'
    sorted_lens_path = 'pan2gene_sorted_lens.pkl'
    pos_path = 'pos/pos.pkl.chr'
    stat_file = 'stat.info'




    [paths, ref_paths_indexes] = read_gfa(gfa_path, chrom_names, seg_index_path,
                                                            path_index_path, chrom_index_path, link_index_path)
    print('gfa-file finished')
    ref_paths = get_ref_paths(gfa_path, ref_paths_indexes)
    print('ref path geholt')

    [sorted_lens, ref_lens] = get_ref_seg_lens(ref_paths, seg_index_path)
    print('ref seg lens berechnet')

    with open(ref_lens_path, 'wb') as ind:
        pickle.dump(ref_lens, ind)
    with open(sorted_lens_path, 'wb') as ind:
        pickle.dump(sorted_lens, ind)
