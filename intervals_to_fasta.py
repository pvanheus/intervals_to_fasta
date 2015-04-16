#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import re
import operator
import Bio.SeqIO
import pvh.gff_utils

def write_gene(output_file, seq, start, end):
    #print(start, end)
    subseq = seq[start:end]
    # write genomic (1-based) coordinates of ig region start and end
    subseq.id = '{}_{}_{}'.format(id_from_gb(seq), start+1, end)
    subseq.description = '{}-{} '.format(start+1, end) + subseq.description
    Bio.SeqIO.write(subseq, output_file, 'fasta')

GB_RE = re.compile(r'gi\|\d+\|\w+\|(?P<id>[^|]+)\|')
def id_from_gb(gb_id):
    match = GB_RE.match(gb_id.id)
    if match == None:
        raise ValueError("{} doesn't look like a Genbank ID".format(gb_id.id))
    return match.group('id')

def add_to_coordinates(genome_id, strand, coordinates, new_coordinate):
    assert strand in ('+', '-'), "Unknown strand: {} expected + or -".format(strand)
    coordinates_by_strand = coordinates.get(genome_id, {'+': [], '-': []})
    coordinate_list = coordinates_by_strand[strand]
    coordinate_list.append(new_coordinate)
    coordinates_by_strand[strand] = coordinate_list
    coordinates[genome_id] = coordinates_by_strand

def intervals_to_fasta(fasta_input_file, gff_input_file, output_file, min_intergenic_length=1,
                       pad_start=0, pad_end=0):
    count = 0
    start = end = -1
    coordinates = dict()
    gene_locations = dict()
    for line in gff_input_file:
        count += 1
        if line.startswith('#'):
            continue
        line = line.rstrip()
        fields = line.split('\t')
        if len(fields) != 9:
            print("GFF3 format error on line {}:\n{}".format(count, line), file=sys.stderr)
            continue
        if fields[2] == 'gene':
            attributes = pvh.gff_utils.parse_gff_attributes(fields[8])
            gene_id = attributes["ID"]
            start = int(fields[3]) - 1 # want to stop 1 before the start of the gene, so this is end of ig region
            end = int(fields[4]) # don't need to specify end+1 because GFF coordinates are 1 based. ig region starts here
            strand = fields[6]
            # cache the gene location - for eukaryotes this will be a superset of the CDS
            gene_locations[gene_id] = (start, end, strand)
            # rotate coordinates one to the left
        elif fields[2] == 'CDS':
            genome_id = fields[0]
            attributes = pvh.gff_utils.parse_gff_attributes(fields[8])
            parent_id = attributes['Parent']
            (start, end, strand) = gene_locations[parent_id]
            coordinate = (start, end)
            if strand == '.':
                strand = '+' # assume positive strand if none specified
            add_to_coordinates(genome_id, strand, coordinates, coordinate)
    genome_dict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(fasta_input_file, 'fasta'), key_function=id_from_gb)            
    for genome_id in sorted(coordinates.keys()):
        for strand in coordinates[genome_id]:
            # sort on start coordinate
            coordinate_list = sorted(coordinates[genome_id][strand], key=operator.itemgetter(0))
            for i in range(len(coordinate_list)-1):
                start = coordinate_list[i][1] - pad_start
                # we need to check that start and end don't run off the end of the sequence because
                # even though our intervals are intergenic, the first or last gene feature
                # might be smaller than pad_start or pad_end
                if start < 0:
                    start = 0
                end = coordinate_list[i+1][0] + pad_end
                if end > (len(genome_dict[genome_id]) - 1):
                    end = (len(genome_dict[genome_id]) - 1)
                if (end - start) >= min_intergenic_length:
                    write_gene(output_file, genome_dict[genome_id], start, end)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract intervals from GFF/FASTA description of Genome')
    parser.add_argument('--pad_start', default=0, type=int, help="Numbers of extra bases to add at start of interval")
    parser.add_argument('--pad_end', default=0, type=int, help="Number of extra bases to add at end of interval")
    parser.add_argument('fasta_input_file', type=argparse.FileType(), help='FASTA format description of genome')
    parser.add_argument('gff_input_file', type=argparse.FileType(), help='GFF3 format description of genome')
    parser.add_argument('output_file', default=sys.stdout, nargs='?', type=argparse.FileType('w'), 
        help='Output file (FASTA format)')
    args = parser.parse_args()

    intervals_to_fasta(args.fasta_input_file, args.gff_input_file, args.output_file, args.pad_start, args.pad_end)