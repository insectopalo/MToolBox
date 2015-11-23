#!/usr/bin/env python

#yay

import ast
import re, sys, glob, os, getopt


from operator import itemgetter
from itertools import repeat
from classifier import tree, NGclassify, consts, parse_mhcs
from classifier.datatypes2 import *
from collections import Counter

# folder where to find data for haplogroup classification and functional annotation
data_file = os.path.dirname(sys.argv[0])

#dizionario ambiguita'
dIUPAC={'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','C','G','T']}

# mt coding loci. This list is used to discriminate, among non-patho_table variants,
# those belonging to coding and non-coding regions to correctly annotate "syn" and "frameshift" mutations
mt_coding_loci = ['MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MT-ND4L', 'MT-ATP8', 'MT-ATP6', 'MT-CYB']

#sequenza di riferimento RSRS
RSRS = 'GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCCCATCCCATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTAAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAATCTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACNNCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGTTCACCCTCTAAATCACCACGATCAAAAGGGACAAGCATCAAGCACGCAACAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAAACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTCGCCAGAACACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATATCCCTCTAGAGGAGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATACCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAGACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTGAGAGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTCAAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGTCGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGAGCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAGTTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTCCAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTAACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCACTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATCACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATGACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAATACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACCATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACNTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTTGTTAAGATGGCAGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTCAGAGGTTCAATTCCTCTTCTTAACAACATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTCTAATCGCAATGGCATTCCTAATGCTTACCGAACGAAAAATTCTAGGCTATATACAACTACGCAAAGGCCCCAACGTTGTAGGCCCCTACGGGCTACTACAACCCTTCGCTGACGCCATAAAACTCTTCACCAAAGAGCCCCTAAAACCCGCCACATCTACCATCACCCTCTACATCACCGCCCCGACCTTAGCTCTCACCATCGCTCTTCTACTATGAACCCCCCTCCCCATACCCAACCCCCTGGTTAACCTCAACCTAGGCCTCCTATTTATTCTAGCCACCTCTAGCCTAGCCGTTTACTCAATCCTCTGATCAGGGTGAGCATCAAACTCAAACTACGCCCTGATCGGCGCACTGCGAGCAGTAGCCCAAACAATCTCATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTTCATAGCCGAATACACAAACATTATTATAATAAACACCCTCACCACTACAATCTTCCTAGGAACAACATATGACGCACTCTCCCCTGAACTCTACACAACATATTTTGTCACCAAGACCCTACTTCTGACCTCCCTGTTCTTATGAATTCGAACAGCATACCCCCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATACCCATTACAATCTCCAGCATTCCCCCTCAAACCTAAGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGTTTAAACCCCCTTATTTCTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAATTAATCCCCTGGCCCAACCCGTCATCTACTCTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCTCGCACTGATTTTTTACCTGAGTAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCAAAAAAATAAACCCTCGTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTCTAATAGCTATCCTCTTCAACAATATACTCTCCGGACAATGAACCATAACCAATACTACCAATCAATACTCATCATTAATAATCATAATGGCTATAGCAATAAAACTAGGAATAGCCCCCTTTCACTTCTGAGTCCCAGAGGTTACCCAAGGCACCCCTCTGACATCCGGCCTGCTTCTTCTCACATGACAAAAACTAGCCCCCATCTCAATCATATACCAAATCTCTCCCTCACTAAACGTAAGCCTTCTCCTCACTCTCTCAATCTTATCCATCATAGCAGGCAGTTGAGGTGGATTAAACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACATAGGATGAATAATAGCAGTTCTACCGTACAACCCTAACATAACCATTCTTAATTTAACTATTTATATTATCCTAACTACTACCGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTATCTCGCACCTGAAACAAGCTAACATGACTAACACCCTTAATTCCATCCACCCTCCTCTCCCTAGGAGGCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATTATCGAAGAATTCACAAAAAACAATAGCCTCATCATCCCCACCATCATAGCCACCATCACCCTCCTTAACCTCTACTTCTACCTACGCCTAATCTACTCCACCTCAATCACACTACTCCCCATATCTAACAACGTAAAAATAAAATGACAGTTTGAACATACAAAACCCACCCCATTCCTCCCCACACTCATCGCCCTTACCACGCTACTCCTACCTATCTCCCCTTTTATACTAATAATCTTATAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTGTAACAGCTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTACTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAACTGGCTTCAATCTACTTCTCCCGCCGCCGGGAAAAAAGGCGGGAGAAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAATCACCTCGGAGCTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCACTCAGCCATTTTACCTCACCCCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCTCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCGCTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAATATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTAGACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAAGATATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCTAAATCCTATATATCTTAATGGCACATGCAGCGCAAGTAGGTCTACAAGACGCTACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGATCACGCCCTCATAATCATTTTCCTTATCTGCTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTAACTAATACTAACATCTCAGACGCTCAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCATCCTAGTCCTCATCGCCCTCCCATCCCTACGCATCCTTTACATAACAGACGAGGTCAACGATCCCTCCCTTACCATCAAATCAATTGGCCACCAATGGTACTGAACCTACGAGTACACCGACTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTAGAACCAGGCGACCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATTCGTATAATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAATTCCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATACTACGGTCAATGCTCTGAAATCTGTGGAGCAAACCACAGTTTCATGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGAAATGCCCCAACTAAATACTACCGTATGGCCCACCATAATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAAACTACCACTTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGAACCAAAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCCTAGGCCTACCCGCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACTAATTACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATAGCCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACCCTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTACTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTAAGCCTCTACCTGCACGACAACACATAATGACCCACCAATCACATGCCTATCATATAGTAAAACCCAGCCCATGACCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAGCCATGTGATTTCACTTCCACTCCATAACGCTCCTCATACTAGGCCTACTAACCAACACACTAACCATATACCAATGATGGCGCGATGTAACACGAGAAAGCACATACCAAGGCCACCACACACCACCTGTCCAAAAAGGCCTTCGATACGGGATAATCCTATTTATTACCTCAGAAGTTTTTTTCTTCGCAGGATTTTTCTGAGCCTTTTACCACTCCAGCCTAGCCCCTACCCCCCAACTAGGAGGGCACTGGCCCCCAACAGGCATCACCCCGCTAAATCCCCTAGAAGTCCCACTCCTAAACACATCCGTATTACTCGCATCAGGAGTATCAATCACCTGAGCTCACCATAGTCTAATAGAAAACAACCGAAACCAAATAATTCAAGCACTGCTTATTACAATTTTACTGGGTCTCTATTTTACCCTCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGGACTTCACGTCATTATTGGCTCAACTTTCCTCACTATCTGCTTCATCCGCCAACTAATATTTCACTTTACATCCAAACATCACTTTGGCTTCGAAGCCGCCGCCTGATACTGGCATTTTGTAGATGTGGTTTGACTATTTCTGTATGTCTCCATCTATTGATGAGGGTCTTACTCTTTTAGTATAAATAGTACCGTTAACTTCCAATTAACTAGTTTTGACAACATTCAAAAAAGAGTAATAAACTTCGCCTTAATTTTAATAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACTACCACAACTCAACGGCTACATAGAAAAATCCACCCCTTACGAGTGCGGCTTCGACCCTATATCCCCCGCCCGCGTCCCTTTCTCCATAAAATTCTTCTTAGTAGCTATTACCTTCTTATTATTTGATCTAGAAATTGCCCTCCTTTTACCCCTACCATGAGCCCTACAAACAACTAACCTGCCACTAATAGTTATGTCATCCCTCTTATTAATCATCATCCTAGCCCTAAGTCTGGCCTATGAGTGACTACAAAAAGGATTAGACTGAGCCGAATTGGTATATAGTTTAAACAAAACGAATGATTTCGACTCATTAAATTATGATAATCATATTTACCAAATGCCCCTCATTTACATAAATATTATACTAGCATTTACCATCTCACTTCTAGGAATACTAGTATATCGCTCACACCTCATATCCTCCCTACTATGCCTAGAAGGAATAATACTATCGCTGTTCATTATAGCTACTCTCATAACCCTCAACACCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTTTTTGCCGCCTGCGAAGCAGCGGTAGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGTACATAACCTAAACCTACTCCAATGCTAAAACTAATCGTCCCAACAATTATATTACTACCACTGACATGACTCTCCAAAAAACACATAATTTGAATCAACACAACCACCCACAGCCTAATTATTAGCATCATCCCCCTACTATTTTTTAACCAAATCAACAACAACCTATTTAGCTGCTCCCCAACCTTTTCCTCCGACCCCCTAACAACCCCCCTCCTAATACTAACTACCTGACTCCTACCCCTCACAATCATGGCAAGCCAACGCCACTTATCCAGTGAACCACTATCACGAAAAAAACTCTACCTCTCTATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCACACTTATCCCCACCTTGGCTATCATCACCCGATGAGGCAACCAGCCAGAACGCCTGAACGCAGGCACATACTTCCTATTCTACACCCTAGTAGGCTCCCTTCCCCTACTCATCGCACTAATTTACACTCACAACACCCTAGGCTCACTAAACATTCTACTACTCACTCTCACTGCCCAAGAACTATCAAACTCCTGAGCCAACAACTTAATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATACGCCTCACACTCATTCTCAACCCCCTGACAAAACACATAGCCTACCCCTTCCTTGTACTATCCCTATGAGGCATAATTATAACAAGCTCCATCTGCCTACGACAAACAGACCTAAAATCGCTCATTGCATACTCTTCAATCAGCCACATAGCCCTCGTAGTAACAGCCATTCTCATCCAAACCCCCTGAAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACATTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTTCAAACTAGACTACTTCTCCATAATATTCATCCCTGTAGCATTGTTCGTTACATGGTCCATCATAGAATTCTCACTGTGATATATAAACTCAGACCCAAACATTAATCAGTTCTTCAAATATCTACTCATTTTCCTAATTACCATACTAATCTTAGTTACCGCTAACAACCTATTCCAACTGTTCATCGGCTGAGAGGGCGTAGGAATTATATCCTTCTTGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCCATTCAAGCAATCCTATACAACCGTATCGGCGATATCGGTTTCATCCTCGCCTTAGCATGATTTATCCTACACTCCAACTCATGAGACCCACAACAAATAGCCCTTCTAAACGCTAATCCAAGCCTCACCCCACTACTAGGCCTCCTCCTAGCAGCAGCAGGCAAATCAGCCCAATTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAGTCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAGTAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTATTCCAAAGACCACATCATCGAAACCGCAAACATATCATACACAAACGCCTGAGCCCTATCTATTACTCTCATCGCTACCTCCCTGACAAGCGCCTATAGCACTCGAATAATTCTTCTCACCCTAACAGGTCAACCTCGCTTCCCTACCCTTACTAACATTAACGAAAATAACCCCACCCTACTAAACCCCATTAAACGCCTGGCAGCCGGAAGCCTATTCGCAGGATTTCTCATTACTAACAACATTTCCCCCGCATCCCCCTTCCAAACAACAATCCCCCTCTACCTAAAACTCACAGCCCTCGCTGTCACTTTCCTAGGACTTCTAACAGCCCTAGACCTCAACTACCTAACCAACAAACTTAAAATAAAATCCCCACTATGCACATTTTATTTCTCCAACATACTCGGATTCTACCCTAGCATCACACACCGCACAATCCCCTATCTAGGCCTTCTTACGAGCCAAAACCTGCCCCTACTCCTCCTAGACCTAACCTGACTAGAAAAGCTATTACCTAAAACAATTTCACAGCACCAAATCTCCACCTCCATCATCACCTCAACCCAAAAAGGCATAATTAAACTTTACTTCCTCTCTTTCTTCTTCCCACTCATCCTAACCCTACTCCTAATCACATAACCTATTCCCCCGAGCAATCTCAATTACAATATATACACCAACAAACAATGTTCAACCAGTAACTACTACTAATCAACGCCCATAATCATACAAAGCCCCCGCACCAATAGGATCCTCCCGAATCAACCCTGACCCCTCTCCTTCATAAATTATTCAGCTTCCTACACTATTAAAGTTTACCACAACCACCACCCCATCATACTCTTTCACCCACAGCACCAATCCTACCTCCATCGCTAACCCCACTAAAACACTCACCAAGACCTCAACCCCTGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCGCTGTAGTATATCCAAAGACAACCATCATTCCCCCTAAATAAATTAAAAAAACTATTAAACCCATATAACCTCCCCCAAAATTCAGAATAATAACACACCCGACCACACCGCTAACAATCAATACTAAACCCCCATAAATAGGAGAAGGCTTAGAAGAAAACCCCACAAACCCCATTACTAAACCCACACTCAACAGAAACAAAGCATACATCATTATTCTCGCACGGACTACAACCACGACCAATGATATGAAAAACCATCGTTGTATTTCAACTACAAGAACACCAATGACCCCAATACGCAAAATTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAATGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGGATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGCAACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGTAATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTGAGGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTTGCCCTTCATTATTGCAGCCCTAGCAGCACTCCACCTCCTATTCTTGCACGAAACGGGATCAAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACACAATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATTCTCACCAGACCTCCTAGGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCCTCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCCTAACAAACTAGGAGGCGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCCCATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTATTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTACCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCTAATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGGCCTGTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGAGAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACAGTACCATAAATACTTGACCACCTGTAGTACATAAAAACCCAATCCACATCAAAACCCTCCCCCCATGCTTACAAGCAAGTACAGCAATCAACCTTCAACTGTCACACATCAACTGCAACTCCAAAGCCACCCCTCACCCACTAGGATATCAACAAACCTACCCACCCTTAACAGTACATAGCACATAAAGCCATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCCTCAGATAGGGGTCCCTTGACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCTACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG'
rsrs_seq = consts.RCRS

#def pato_lines(output_list, pato_list)
#    # pato dictionary has values that are lists of list.
#    for i in pato_list:
#        output_list.append

def check_mutation_effect(event, sitevar_table_row):
    # check if the variant maps in coding regions
    if sitevar_table_row[0] in mt_coding_loci:
        if event.mutation_type() not in ["Insertion", "Deletion"]:
            sitevar_table_row[3] = "syn"
        else:
            sitevar_table_row[3] = "frameshift"
    return sitevar_table_row


def pprint2datatype(event):
    """This function takes as input mutation events as they are pprinted (es. Deletion(514) is pprinted as 514-515d,
    given a deletion span of 2nt) and re-convert them to datatypes instances.
    Ambiguities must be taken into account: their pprint is 310C(Y), the ambiguity should be discarded
    """
    if '(' in event: # ambiguity
        event = event.split('(')[0]
    if event[-1] == "d": #514-515d
        e = event[:-1].split('-')
        if len(e) == 1:
            mut = Deletion("%d-%dd" % (int(e[0]), int(e[0])))
        else: # else, it should be only 2, no control so far and we have finger crossed
            mut = Deletion("%d-%dd" % (int(e[0]), int(e[1])))
    elif '.' in event: #310.CC
        e = event.split('.')
        mut = Insertion("%d.%s" % (int(e[0]), e[1]))
    else: #15746C
        pos = int(event[:-1])
        var = event[-1]
        # SNP_MixIn class does not take parameters, so attributes are added after its initialization
        mut = SNP_MixIn()
        mut.start = pos
        mut.change = var
## Difference between datatypes is not considered anymore. 15th May 2014.
#        ref = rsrs_seq[pos-1]
#        if (ref in consts.PUR and var in consts.PUR) or (ref in consts.PYR and var in consts.PYR):
#            mut = Transition(pos, ref=ref)
#            print mut.pprint()
#        elif (ref in consts.PUR and var in consts.PYR) or (ref in consts.PYR and var in consts.PUR):
#            mut = Transversion("%d%c" % (pos, var), ref=ref)
#        else:
#            mut = Unknown(pos)
#            mut.change = var
    return mut

def feature2datatype(dict, snp, ref_seq=rsrs_seq):
    """output should be a list [Transition(100) : [HF, CI_low, CI_up], Transition(101) : [HF, CI_low, CI_up]]
    """
    for j in zip(snp[3], snp[6], snp[7], snp[8]):
        if snp[-1] == 'del': #OK?
            # [515, ['ACA'], 525, ['A'], [10], [35.0], [0.019], [0.01176], 'del']
            mut = Deletion("%d-%dd" % (snp[0]+1, snp[0]+1+len(snp[1][0])-2))
        elif snp[-1] == 'mism': #OK
            # [73, ['G'], 'mism'] ma in realta' vorrei
            # [73, ['A', 'G'], 'mism'] dove 'A' sarebbe il nt in RSRS, 'G' la mutazione
            ref = ref_seq[snp[0]-1]
            var = j[0]
            #var = snp[1][0]
            # print snp[0]-1, ref, var
            # print (ref in consts.PUR and var in consts.PUR) or (ref in consts.PYR and var in consts.PYR)
            if (ref in consts.PUR and var in consts.PUR) or (ref in consts.PYR and var in consts.PYR):
                mut = Transition(snp[0], ref=ref)
            elif (ref in consts.PUR and var in consts.PYR) or (ref in consts.PYR and var in consts.PUR):
                mut = Transversion("%d%c" % (snp[0], var), ref=ref)
            else:
                mut = Unknown(snp[0])
                mut.change = var
        elif snp[-1] == 'ins':
            # [309, ['CCT'], 'ins']
            mut = Insertion("%d.%s" % (snp[0], snp[1]))
        try:
            dict[mut.pprint()] = [j[1], j[2], j[3]]
        except:
            print mut.pprint(), "Error", j
            sys.exit()

#def feature2datatype_old(dict, snp, ref_seq=rsrs_seq):
#    """output should be a list [Transition(100) : [HF, CI], Transition(101) : [HF, CI]]
#        """
#    for j in zip(snp[3], snp[6], snp[7]):
#        if snp[-1] == 'del': #OK?
#            # [515, ['ACA'], 525, ['A'], [10], [35.0], [0.019], [0.01176], 'del']
#            mut = Deletion("%d-%dd" % (snp[0]+1, snp[0]+1+len(snp[1][0])-2))
#        elif snp[-1] == 'mism': #OK
#            # [73, ['G'], 'mism'] ma in realta' vorrei
#            # [73, ['A', 'G'], 'mism'] dove 'A' sarebbe il nt in RSRS, 'G' la mutazione
#            ref = ref_seq[snp[0]-1]
#            var = j[0]
#            #var = snp[1][0]
#            # print snp[0]-1, ref, var
#            # print (ref in consts.PUR and var in consts.PUR) or (ref in consts.PYR and var in consts.PYR)
#            if (ref in consts.PUR and var in consts.PUR) or (ref in consts.PYR and var in consts.PYR):
#                mut = Transition(snp[0])
#            elif (ref in consts.PUR and var in consts.PYR) or (ref in consts.PYR and var in consts.PUR):
#                mut = Transversion("%d%c" % (snp[0], var))
#            else:
#                mut = Unknown(snp[0])
#                mut.change = var
#        elif snp[-1] == 'ins':
#            # [309, ['CCT'], 'ins']
#            mut = Insertion("%d.%s" % (snp[0], snp[1]))
#        dict[mut.pprint()] = [j[1], j[2]]


#definisco come riempire i dizionari
def fillDict(inhandle, dict, separator, listindex):
    for line in inhandle:
        x = line.strip().split(separator)
        a=x[0]
        b=[]
        b.extend(x[listindex:])
        dict[a]=b
    return dict

## special dictionary which takes datatypes.SNP as keys
#def fillPosDict(inhandle, dict, separator, listindex):
#    for line in inhandle:
#        x = line.strip().split(separator)
#        a=pprint2datatype(x[0])
#        b=[]
#        b.extend(x[listindex:])
#        dict[a]=b
#    return dict

# special dictionary which takes datatypes.SNP as keys
def fillPosDict(inhandle, dict, separator, listindex):
    for line in inhandle:
        x = line.strip().split(separator)
        a=pprint2datatype(x[0])
        b=[]
        b.extend(x[listindex:])
        dict[a]=b
    return dict

# special dictionary for patho table, since one position could be annotated for more than one locus
def fillPathoDict2(inhandle, dict, separator, listindex):
    for line in inhandle:
        x = line.strip().split(separator)
        a=pprint2datatype(x[0])
        b = x[listindex:]
        for ind,i in enumerate(b): # convert fields to float whenever possible
            try:
                i = float(i)
                b[ind] = i
            except:
                pass
        try:
            b[2] = int(b[2])
        except:
            #print e, len(b), x[0], b
            pass
        dict.setdefault(a, []).append(b)
    return dict

# special dictionary for sitevar table, keys must be integers
def fillSiteVarDict(inhandle, dict, separator, listindex):
    for line in inhandle:
        x = line.strip().split(separator)
        a=x[0]
        b=[]
        b.extend(x[listindex:])
        dict[int(a)]=b
    return dict

#definisco come ripulire le stringhe
def replace_all(text, dic):
    for r, t in dic.iteritems():
        text = text.replace(r, t)
    return text

def data_parsing(file_file, site_file, bestres_file, haptab_file):
    #apro i file di input
    #diff = open(diff_file, 'r')
    file = open(file_file, 'r')
    site = open(site_file, 'r')
    bestres = open(bestres_file, 'r')
    haptab = open(haptab_file, 'r')

    print "Parsing pathogenicity table..."
    #dizionario del tabellone pato
    # scarto la prima riga (intestazione)
    file = file.readlines()[1:]
    d={}
    d = fillPathoDict2(file, d, '\t', 2)
    
    print "Parsing variability data..."
    #dizionario sitevar
    # scarto la prima riga (intestazione)
    site = site.readlines()[1:]
    g={}
    g = fillSiteVarDict(site, g, '\t', 1)
    for pos in g.keys():
        # convert variability values (nt and aa, if available) in floats
        try:
            g[pos][1] = float(g[pos][1])
        except:
            pass
        try:
            g[pos][4] = float(g[pos][4])
        except:
            pass

    print "Parsing info about haplogroup-defining sites..."
    haplo = {}
    htree = tree.HaplogroupTree(pickle_data=open(data_file +'/data/phylotree_r16.pickle', 'rb').read())
    for haplogroup in htree._aplo_dict.keys():
        haplo[haplogroup] = []
        # change Transition and Transversion datatypes to SNP_MixIn,
        # because Transition and Transversion variants parsed from merged_diff are that type
        for event in htree.get_filtered_positions(haplogroup):
                if event.mutation_type() in ['Transition', 'Transversion']:
                    event_new = SNP_MixIn()
                    event_new.start = event.start
                    event_new.change = event.change
                    haplo[haplogroup].append(event_new)
                else:
                    haplo[haplogroup].append(event)
    # la funzione di cui sopra sostituisce la procedura qui sotto
    #haplo = [line.strip().split('\t') for line in haptab] #lista di liste aplogruppo/variante
    # DS .. convert strings of variants (eg '146T') in datatypes.SNP (eg Transition(146))
    #for i in haplo[1:]:
    #    i[1] = pprint2datatype(i[1])
    #print haplo[:4]
    #hapconto = hapcon.read() #file letto interamente che usero' per contare quante volte e' presente una variante
    haplo_sites = []
    for haplogroup in haplo.keys():
        haplo_sites.extend(haplo[haplogroup])
    hapconto = Counter(haplo_sites)
    #hapcon = open(haptab_file, 'r')
    #hapcont = [line.strip().split('\t') for line in hapcon]
    #hapconto = []
    #for l in hapcont:
    #    for e in l:
    #        hapconto.append(e)
    #print hapconto.keys()[:10]
    #print "7146A is found %d times in hapconto" % (hapconto.count('7146A'))

    print "Parsing info about haplogroup assignments..."
    #dizionario best results
    # scarto la prima riga (intestazione)
    bestres = bestres.readlines()[1:]
    best={}
    best = fillDict(bestres, best, ',', 1)
    #for s in best.keys():
    #    best[s] = [best[s][0].split(';')[0]]
    #print best
    return d, g, haplo, hapconto, best

#def main_functional_analysis(diff_file, file_file, site_file, bestres_file, haptab_file, PATH, FILENAME):
def get_HF(VCF_dict, sampleID, variant_site):
    # variant site is pprint
    HF, CI = ('', '')
    try:
        HF, CI_low, CI_up = VCF_dict[sampleID][variant_site]
        CI = ';'.join(str(item) for item in [CI_low, CI_up])
    except:
        #print "sample %s, variant %s not found." % (sampleID, variant_site)
        pass
    return HF, CI
#     main_functional_analysis(merg, d, g, best, haplo, hapconto, PATH, FILENAME, mutations_HF_dict)
def main_functional_analysis(merg, d, g, best, haplo, hapconto, PATH, FILENAME, mutations_HF_dict, samp):
    ###
    ###
    #diff = open(sys.argv[1],'r') #file diff
    #file = open(sys.argv[2],'r') #tabellone patogenicita'
    #site = open(sys.argv[3],'r') #tabella sitevar
    #bestres = open(sys.argv[4],'r') #file best results
    #haptab = open(sys.argv[5],'r') #tabella aplotipi
    #hapcon = open(sys.argv[5],'r') #tabella aplotipi da capo


    #prendo l'aplogruppo relativo nel best results
    #hapg=best[samp[0]][0]
    #hapgs = best[samp[0]][0].split(';')
    # DS
    #try:
    hapgs = best[samp][0].split(';')
    print "Best haplogroup predictions for sample", samp, ":", hapgs
    for hapgro in hapgs:
        print "Functional annotation for haplogroup", hapgro
        #apro il file di output
        #o = open("output."+samp[0]+hapgro+'.txt', 'w')
        #o = open('.'.join([sys.argv[1].split('.')[0], 'extended.csv']), 'w')
        # DS
        o = open(os.path.join(PATH, '.'.join([samp, hapgro, 'annotation.csv'])), 'w')
        #scrivo l'intestazione
        #aggiungere campi:
        # 
        header = ['Sample',"Variant Allele",'HF','CI_lower;CI_upper','RSRS','MHCS','rCRS','Haplogroup','Other Haplogroups','Locus','Nt Variability','Codon Position','Aa Change','Aa Variability','tRNA Annotation','Disease Score','RNA predictions','MutPred pred','MutPred prob','PolyPhen-2 HumDiv pred','PolyPhen-2 HumDiv prob','PolyPhen-2 HumVar pred','PolyPhen-2 HumVar prob','PANTHER pred','PANTHER prob','PhD-SNP pred','PhD-SNP prob','SNPs&GO pred','SNPs&GO prob','Mitomap Associated Disease(s)','Mitomap Homoplasmy','Mitomap Heteroplasmy','Somatic Mutations','SM Homoplasmy','SM Heteroplasmy','ClinVar','OMIM link','dbSNP ID','Mamit-tRNA link','PhastCons20Way','PhyloP20Way','AC/AN 1000 Genomes','1000 Genomes Homoplasmy','1000 Genomes Heteroplasmy'] 
        o.write('\t'.join(header)+'\n')


    # lista delle varianti che definiscono l'aplogruppo best
        varbest=[]
        for i in merg: # for each variant
            if mutations_HF_dict != {}:
                #print i.pprint()
                HF, CI = get_HF(mutations_HF_dict, samp, i.pprint())
                #print "%s, HF: %s, CI: %s" % (i.pprint(), HF, CI)
            else:
                HF, CI = ('', '')
            if i in haplo[hapgro]: # if it defines best haplogroup
                if hapconto[i] > 1: # if it defines also other haplogroups
                    try:
                        special_data = d[i] # data associated to variant in patho_table; if variant in patho table; rename it "special_data"
                        # da da da...
                    except KeyError:
                        special_data = [g[i.start]]
                        for row in special_data:
                            row = check_mutation_effect(i, row)
                    for row in special_data:
                        varbest_element = []
                        for element in [[samp, i.pprint()], [HF, CI], merg[i], [hapgro], ['+'], row]:
                            varbest_element.extend(element)
                elif hapconto[i] == 1:
                    try:
                        special_data = d[i]
                    except KeyError:
                        special_data = [g[i.start]]
                        for row in special_data:
                            row = check_mutation_effect(i, row)
                    for row in special_data:
                        varbest_element = []
                        for element in [[samp, i.pprint()], [HF, CI], merg[i], [hapgro], [''], row]:
                            varbest_element.extend(element)
            else: # if the variant does not define the best haplogroup
                if hapconto[i] > 0: # if it defines at least one haplogroup
                    try:
                        special_data = d[i] # data associated to variant in patho_table; if variant in patho table; rename it "special_data"
                        # da da da...
                    except KeyError:
                        special_data = [g[i.start]]
                        for row in special_data:
                            row = check_mutation_effect(i, row)
                    for row in special_data:
                        varbest_element = []
                        for element in [[samp, i.pprint()], [HF, CI], merg[i], [''], ['+'], row]:
                            varbest_element.extend(element)
                else: # if the variant does not define any haplogroup
                    try:
                        special_data = d[i]
                    except KeyError:
                        special_data = [g[i.start]]
                        for row in special_data:
                            row = check_mutation_effect(i, row)
                    for row in special_data:
                        varbest_element = []
                        for element in [[samp, i.pprint()], [HF, CI], merg[i], [''], [''], row]:
                            varbest_element.extend(element)
            #print "SPECIAL DATA:", special_data
            varbest_element.extend(repeat('', len(header)-len(varbest_element))) # padding list (row) content
            varbest.append(varbest_element)
            #print i
            #print "change is", i.print_snp()
            # if hapconto.count(i) > 1: #controllo per vedere se la variante definisce altri aplogruppi
            # DS: since i is a datatypes instance, I check its pprint in hapconto (which is a brutal string)

        #for x,i in enumerate(ord):
        #    print x,i
        #sys.exit()
        #mut = sorted(ord, key=itemgetter(14), reverse=True) # mutpred probability
        #mut = sorted(allvar, key=itemgetter(16), reverse=True) # mutpred probability
        #mut = sorted(varbest, key=itemgetter(16), reverse=True) # mutpred probability
        dis = sorted(varbest, key=itemgetter(15), reverse=True) # % disease (ex-mean pathogenicity score) 
        vrb = sorted(dis, key=itemgetter(10)) # nt variability
        hpg = sorted(vrb, key=itemgetter(7)) # genomes haplogroup defining
        mrg = sorted(hpg, key=itemgetter(5), reverse=True) # MHCS

        for i in mrg:
            i = [str(el) for el in i]
            o.write('\t'.join(i)+'\n')

        o.close()
#    except Exception, e:
#        print dir(Exception), e,
#        print "Error for sample:", samp
#        pass
if __name__ == '__main__':
    ### things to put in haplo_driver
    # USAGE
    # variants_functional_annotation.py ../mt_classification_best_results.csv
    path = os.getcwd()
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:p:s:t:v:")
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit()    
    bestres_file = 'mt_classification_best_results.csv'
    patho_file = data_file+'/data/patho_table.txt'
    #        file_file = open(sys.argv[1],'r') #tabellone patogenicita'
    site_file = data_file+'/data/sitevar_modified.txt'
    #        site_file = open(sys.argv[2],'r') #tabella sitevar
    #        bestres_file = os.path.join(path, folder, 'mt_classification_best_results')
    #        bestres_file = open(os.path.join(path, 'test', 'mt_classification_best_results'))
    #        bestres_file = open(sys.argv[3],'r') #file best results
    haptab_file = data_file+'/data/haplogroups.txt'
    #        haptab_file = open(sys.argv[3],'r') #tabella aplotipi
    #hapcon = open(sys.argv[5],'r') #tabella aplotipi da capo
    VCF_dict_tmp = 'VCF_dict_tmp'
    for o,a in opts:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-b": bestres_file = a
        elif o == "-p": patho_file = a
        elif o == "-s": site_file = a
        elif o == "-t": haptab_file = a
        elif o == "-v": VCF_dict_tmp = a
        else: assert False, "Unhandled option."
    # parsing VCF_dict_tmp
    mutations_HF_dict = {}
    try:
        # the resulting mutations_HF_dict will have this structure:
        # {'NA12878': {'3108d': [0.6386, 0.02548], '10688G': [1.0, 0.0], '2259T': [0.998, 0.00196], ...
        # where for each variant site, the first value in square bracket is HF and the second one is CI
        VCF_dict = ast.literal_eval(open(VCF_dict_tmp, 'r').read())
        mutations_HF_dict = {}
        for sampleID in VCF_dict.keys():
            mutation_events = VCF_dict[sampleID]
            mutations_HF_dict[sampleID] = {}
            for mut_feature_list in mutation_events:
                feature2datatype(mutations_HF_dict[sampleID], mut_feature_list) # fill 
    except:
        print "Heteroplasmy data file ('VCF_dict_tmp') not found. HF will not be reported in the output."
        pass
    #def data_parsing(file_file, site_file, bestres_file, haptab_file, VCF_dict = {})
    for sampleID in mutations_HF_dict.keys():
        for mut in mutations_HF_dict[sampleID].keys():
            pass
            #print sampleID, mut, mutations_HF_dict[sampleID][mut]
    #sys.exit()
    d, g, haplo, hapconto, best = data_parsing(patho_file, site_file, bestres_file, haptab_file)
    #folder = sys.argv[1]
    #    for infile in glob.glob(os.path.join(path, '*_dir', '*_merged_diff.csv')):

    for infile in glob.glob(os.path.join(path, 'OUT*', '*_merged_diff.csv')):
        (PATH, FILENAME) = os.path.split(infile)
        print infile
        diff_file = infile
        samp = FILENAME.replace('_merged_diff.csv', '')
        print "Parsing variant data for sample %s..." % samp
        #dizionario merged diff
        # scarto la prima riga (intestazione)
        try:
            diff = open(diff_file, 'r').readlines()[1:]
            merg={}
            merg = fillPosDict(diff, merg, ',', 1) # DS: merg becomes a dictionary { datatypes.SNP_MixIn/Insertion/Deletion : [yes, yes, yes] }
        except:
            print "Some problem occurred with file %s." % FILENAME
            continue
        main_functional_analysis(merg, d, g, best, haplo, hapconto, PATH, FILENAME, mutations_HF_dict, samp)
