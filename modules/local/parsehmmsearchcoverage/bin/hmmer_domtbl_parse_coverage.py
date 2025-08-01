import argparse
import math
import sys
import re
import fileinput
import json
import re
from collections import defaultdict

parser = argparse.ArgumentParser(description='Parse table output from HMMer (`hmmsearch`) to calculate HMM model coverages from metagenomic reads.')
parser.add_argument('-i', "--input_fp", type=str,
                    default='-',
                    help="Input fasta/fastq filepath. Use '-' for STDIN (default).")
parser.add_argument('-o', "--output_fp", type=str,
                    default='-',
                    help="Output TSV filepath. Use '-' for STDOUT (default).")
parser.add_argument('-d', "--db", type=str,
                    required=True,
                    help="Query HMM database path.")
parser.add_argument('-s', "--stats_output_fp", type=str,
                    default='',
                    help="Output JSON filepath for recording mapping stats.")

args = parser.parse_args()


def get_db_metadata(db_path):
    with open(db_path, 'rt') as f:
        metadata = {}
        for line in f:
            if line[0] in {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}:
                suffix = 'HMM'
                if line[:len(suffix)] == suffix:
                    continue
                line_ = line.strip().split()
                k = line_[0]
                v = ' '.join(line_[1:])
                metadata[k] = v
            if line[:2]=='//':
                yield metadata
                metadata = {}
        else:
            yield metadata

cols = [
    'target',
    'query',
    'target_start',
    'target_end',
    'query_start',
    'query_end',
    'comp_score',
    'comp_bias',
    'cell_evalue',
    'cell_frac',
]

def extract_from_line(line_dict):
    d, b, e = int(line_dict['read_frame']), int(line_dict['read_frame_begin']), int(line_dict['read_frame_end'])
    if d<0:
        b, e = e, b

    return {
        'read_frame': (d,b,e),
        'query': line_dict['query'],
        'comp_score': float(line_dict['comp_score']),
        'cell_evalue': float(line_dict['cell_evalue']),
        'target_start': int(line_dict['target_start']),
        'target_end': int(line_dict['target_end']),
        'query_start': int(line_dict['query_start']),
        'query_end': int(line_dict['query_end'])
    }

if __name__ == '__main__':
    # Parse database
    hmm_lengths = {d['NAME']:int(d['LENG']) for d in get_db_metadata(args.db) if 'NAME' in d}

    # Parse file
    read_hits = defaultdict(list)
    for line in fileinput.input([] if args.input_fp=='-' else args.input_fp):
        if line[0] == '#':
            continue
        line_dict = dict(zip(cols, [v.strip() for v in line.strip().split()]))

        read_header_split = re.findall(r'^(.*?)_frame=(-?\d+)_begin=(\d+)_end=(\d+)\s*$',
                                       line_dict['target'])[0]
        line_dict['read_name'] = read_header_split[0]
        line_dict['read_frame'] = read_header_split[1]
        line_dict['read_frame_begin'] = read_header_split[2]
        line_dict['read_frame_end'] = read_header_split[3]

        read_hits[line_dict['read_name']].append(extract_from_line(line_dict))

    top_read_hits = {}
    for k,vs in read_hits.items():
        # greedy resolution of overlaps
        deoverlapped = []
        ali_coverage = set()
        for d in sorted(vs, key=lambda x:x['cell_evalue']):
            phase = int(d['read_frame'][0])
            direction = -1 if phase<0 else 1
            phase *= direction
            start, end = d['read_frame'][1:3]

            m = lambda x: (start-1)+direction*(x-1)*3 + (phase-1)
            nt_base_idxs = list(range(*list(sorted((m(d['target_start']), m(d['target_end']))))))

            if not any([i in ali_coverage for i in nt_base_idxs]):
                deoverlapped.append(d)
                for i in nt_base_idxs:
                    ali_coverage.add(i)

        top_read_hits[k] = list(deoverlapped)

    # get hmm coverage and read counts
    hmm_hits_coverage = {}
    hmm_hit_count = defaultdict(int)
    for k,vs in top_read_hits.items():
        for d in vs:
            hmm_hit_count[d['query']] += 1
            if not d['query'] in hmm_hits_coverage:
                hmm_hits_coverage[d['query']] = {i+1:0 for i in range(hmm_lengths[d['query']])}
            for i in range(d['query_start'], d['query_end']):
                hmm_hits_coverage[d['query']][i] += 1

    # Collect and write
    hmm_hits_coverage_stats = {}
    for k,d in hmm_hits_coverage.items():
        if not len(d)>0:
            continue
        depth = sum(list(d.values()))/len(d)
        breadth = sum([v>0 for _,v in d.items()])/len(d)

        depth_ = 709 if depth>709 else depth  # prevents and float overflow with math.exp
        expected = 1-(1/math.log2(1+math.exp(depth_)))

        hmm_hits_coverage_stats[k] = {
            'depth': depth,
            'breadth': breadth,
            'count': hmm_hit_count[k],
            'expected_breadth': expected,
            'ratio': breadth/expected,
        }

    outfile = sys.stdout if args.output_fp=='-' else open(args.output_fp, 'wt')
    outfile.write(f"# function\tread_count\tcoverage_depth\tcoverage_breadth\n")
    for k,d in sorted(hmm_hits_coverage_stats.items(), key=lambda x:-x[1]['depth']):
        outfile.write(f"{k}\t{d['count']}\t{d['depth']}\t{d['breadth']}\n")
    outfile.close()

    if args.stats_output_fp:
        stats = {
            'reads_mapped': len(top_read_hits),
            'hmm_count': len(hmm_hit_count),
            'read_hit_count': sum(list(hmm_hit_count.values()))
        }
        with open(args.stats_output_fp, 'wt') as f:
            json.dump(stats, f)


