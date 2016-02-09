#!/usr/bin/env python

import getopt, sys, re, os, glob, csv
from classifier import tree, NGclassify, consts, datatypes, parse_mhcs
from bioinf.seqs import SeqList
import io_modules.csv
import io_modules.old_table
import io_modules.serialize
import os.path

# folder where to find data for haplogroup classification and functional annotation
data_file = os.path.dirname(sys.argv[0])

def usage():
    print """\nAssigns haplogroup to contigs and performs functional annotation
        Options:
        -i        Contig file [mtDNAassembly-Contigs.fasta]
        -m        MUSCLE executable PATH [/usr/local/bin/muscle]
        -b        basename for output files
        -s        file with most reliable haplogroup prediction
        """

def pickle_csv(csvfile, pickle_fname=None):
    tree_file = csv.reader(open(csvfile, 'rb'))
    if pickle_fname is None:
        pickle_fname = csvfile + '.pickle'
    aplo_list = io_modules.csv.parse_csv(tree_file)
    htree = tree.HaplogroupTree(aplo_list=aplo_list)
    pickle_file = open(pickle_fname, 'wb')
    pickle_file.write(htree.serialize())

def write_old_table(pickle_fname, out_fname):
    htree = tree.HaplogroupTree(pickle_data=open(pickle_fname, 'rb').read())
    fh = csv.writer(open(out_fname, 'wb'))
    for haplo_name in htree:
        io_modules.old_table.write_haplogroup(fh, '', htree[haplo_name])

def merge_tables(f, g, h):
    fgh = f + g + h
    mergedlist = []
    for jj in fgh:
        if jj not in mergedlist:
            mergedlist.append(jj)
    o = []
    o.append(["", "RSRS", "MHCS", "rCRS"])
    y = "yes"
    n = ""
    for i in mergedlist:
        if i in f and i in g and i in h:
            o.append([i.pprint(),y,y,y])
        elif i in f and i in g:
            o.append([i.pprint(),y,y,n])
        elif i in f and i in h:
            o.append([i.pprint(),y,n,y])
        elif i in g and i in h:
            o.append([i.pprint(),n,y,y])
        elif i in f:
            o.append([i.pprint(),y,n,n])
        elif i in g:
            o.append([i.pprint(),n,y,n])
        elif i in h:
            o.append([i.pprint(),n,n,y])
    return o

def align_sequence(muscle_exe, sequence, rif=None, ):
    """sequence is a datatypes.Sequence, rif"""
    if rif is None:
        rif = datatypes.Sequence('RSRS', consts.RCRS)
    seq_diff = NGclassify.SequenceDiff()
    seq_diff.gen_diff(muscle_exe, rif, datatypes.Sequence(sequence.name, str(sequence)))
    return seq_diff

def h_analysis(htrees, seq_diff, regions, mhcs_dict):
    a = NGclassify.Classify()
    #print "Classification of sequence %s" % seq_diff.obj.name
    for htree, name in htrees:
        print "Classification according to tree:", name
        a.classify_by_tree(htree, seq_diff, regions)
        #print "start is ", seq_diff.start
        #print "end is ", seq_diff.end
        #print "haplo_stats: ", a.haplo_stats
        print "genome_state is ", a.get_genome_state()
        (haplo_stats_sorted, haplo_best) = a.prediction_sorting()
        #print haplo_best
        #print "haplo_stats_sorted is:\n", haplo_stats_sorted
        print "="*20
        #print "haplo_best is: ", haplo_best
        #print "finding MHCS for sequence %s" % seq_diff.obj.name
        mhcss = a.get_mhcss(mhcs_dict)
        #print "MHCS ID for sequence %s is %s" % (seq_diff.obj.name, ','.join(list(mhcss)))
    print '-'*30
    return a

def load_sequences(fname):
    a = SeqList()
    a.load_file(fname)
    print "Loaded %d contig sequences" % len(a)
    return a

def write_output(class_obj, seq_diff, seq_diff_mhcs, seq_diff_rcrs, merged_tables, outfile):
    print "Writing results for sequence %s" % outfile
    class_obj.pprint(open(outfile + '.csv', 'w'))
    class_obj.pprint_sorted(open(outfile + '.sorted.csv', 'w'))
    merged_tables_file = open(outfile + '_merged_diff.csv', 'w')
    for row in merged_tables:
        merged_tables_file.write(','.join(row)+'\n')

def main_mt_hpred():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:m:b:s:")
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit()
    contig_file = 'mtDNAassembly-contigs.fasta'
    muscle_exe='/usr/local/bin/muscle'
    basename='mtDNAassembly-contigs'
    best_results_file = 'mt_classification_best_results.csv'
    for o,a in opts:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-i": contig_file = a
        elif o == "-m": muscle_exe = a
        elif o == "-b": basename = a
        elif o == "-s": best_results_file = a
        else:
            assert False, "Unhandled option."

    print "Your best results file is ", best_results_file
    # sample name
    f = os.path.abspath(contig_file)
    sample_name = contig_file.split('-')[0]
    
    # haplogroup tree parsing
    htrees = [(tree.HaplogroupTree(pickle_data=open(data_file + '/data/phylotree_r16.pickle', 'rb').read()), data_file + '/data/phylotree_r16.pickle')]
    # mhcs parsing
    mhcs_dict = parse_mhcs.parse2mhcs_dict(data_file + '/data/mhcs.tab')
    
    print "\nLoading contig sequences from file %s" % contig_file
    contig_array = load_sequences(contig_file)
    contig_array_seqdiff = [] # list of lists
    contig_total_seqdiff = [] # list of variants
    contig_array_mappings = []
    
    print "\nAligning Contigs to mtDNA reference genome...\n"
    
    # update each contig's SeqDiff
    for x,contig in enumerate(contig_array):
        if x == 0:
            contig_seq_diff = align_sequence(muscle_exe, contig)
            contig_seq_diff.find_segment() # avoid having long gaps at 5' and 3' (not actual gaps but due to the alignment)
            contig_seq_diff.regions.append([contig_seq_diff.start, contig_seq_diff.end])
        else:
            incoming_seqdiff = align_sequence(muscle_exe, contig)
            incoming_seqdiff.find_segment()
            contig_seq_diff.diff_list.extend(incoming_seqdiff.diff_list)
            contig_seq_diff.regions.append([incoming_seqdiff.start, incoming_seqdiff.end])

    print "\nSequence haplogroup assignment\n"
    seq_classify = h_analysis(htrees, contig_seq_diff, contig_seq_diff.regions, mhcs_dict)
    seq_classify.sample_name = sample_name

    print "Contig alignment to MHCS and rCRS"
    m = list(seq_classify.mhcss)[0]
    print "Aligning contigs to MHCS SeqDiff object"
    its_mhcs = datatypes.Sequence(m, mhcs_dict[m])
    for x, contig in enumerate(contig_array):
        if x == 0:
            contig_mhcs_seq_diff = align_sequence(muscle_exe, contig, its_mhcs)
            contig_mhcs_seq_diff.find_segment()
            contig_mhcs_seq_diff.regions.append([contig_seq_diff.start, contig_seq_diff.end])
        else:
            incoming_mhcs_seqdiff = align_sequence(muscle_exe, contig, its_mhcs)
            incoming_mhcs_seqdiff.find_segment()
            contig_mhcs_seq_diff.diff_list.extend(incoming_mhcs_seqdiff.diff_list)
            contig_mhcs_seq_diff.regions.append([incoming_mhcs_seqdiff.start, incoming_mhcs_seqdiff.end])

    print "rCRS SeqDiff object"
    rcrs = datatypes.Sequence('rCRS', consts.rcrs)
    for x, contig in enumerate(contig_array):
        if x == 0:
            contig_rcrs_seq_diff = align_sequence(muscle_exe, contig, rcrs)
            contig_rcrs_seq_diff.find_segment()
            contig_rcrs_seq_diff.regions.append([contig_seq_diff.start, contig_seq_diff.end])
        else:
            incoming_rcrs_seqdiff = align_sequence(muscle_exe, contig, rcrs)
            incoming_rcrs_seqdiff.find_segment()
            contig_rcrs_seq_diff.diff_list.extend(incoming_rcrs_seqdiff.diff_list)
            contig_rcrs_seq_diff.regions.append([incoming_rcrs_seqdiff.start, incoming_rcrs_seqdiff.end])

    # try gathering diff from reference sequences
    print "Merging seq_diffs..."
    mergedtables = merge_tables(contig_seq_diff.diff_list, contig_mhcs_seq_diff.diff_list, contig_rcrs_seq_diff.diff_list)

    # OUTPUTS
    write_output(seq_classify, contig_seq_diff.diff_list, contig_mhcs_seq_diff.diff_list, contig_rcrs_seq_diff.diff_list, mergedtables, basename)
    open(os.path.join('../', best_results_file), 'a').write(','.join([basename, ';'.join([i[0] for i in seq_classify.haplo_best.items()])])+'\n')
    
if __name__ == "__main__":
    try:
        main_mt_hpred()
    except:
        sys.stderr.write('Unable to compute haplogroup. Exit')
        sys.exit(1)

