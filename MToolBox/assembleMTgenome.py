#!/usr/bin/env python

"""
Written by Ernesto Picardi - e.picardi@biologia.uniba.it
Edited by Claudia Calabrese - claudia.calabrese23@gmail.com
          Domenico Simone - dome.simone@gmail.com
          Gonzalo S Nido - gonzalo.s.nido@gmail.com
"""

import getopt, sys, os, re, ast, math
import mtVariantCaller

def usage():
    print """Assembling MT-DNA from SAM/BAM/Pileup files
Version 2.0 - Written by Ernesto Picardi - 2011-2012
        Edited by Domenico Simone, Claudia Calabrese and Gonzalo Nido - 2013-2016

USAGE: assembleMTgenome.py -i [FILE.bam|FILE.sam|FILE.pileup] [OPTIONS]

Options:
    -i   FILE      Input File (.pileup .sam or .bam)
    -r   PATH      Path to fasta reference genomes [/usr/local/share/genomes]
    -f   FILE      Reference MT-DNA in fasta [chrRSRC.fasta]
    -a   FILE      Human Genome in fasta [hg19RSRC.fasta]
    -q   INT       Min. per base quality score [25]
    -c   FLOAT     Min. confidence level [0.80]
    -d   INT       Min. coverage depth [5]
    -g   INT       Min. gap lentgh for contig [10]
    -o   STRING    Output base name [mtDNAassembly]
    -s   FILE      samtools executable [/usr/local/bin/samtools]
    -t   INT>5     Minimum distance from read end(s) for indels to be detected. Values < 5 will be ignored [5]
    -v   INT       Minimum number of observations of an indel [5]
    -z   FLOAT     Heteroplasmy threshold for variants to be reported in consensus FASTA [0.8]
    -F             Generate fasta output [no]
    -C             Generate coverage file [no]
    -M             Generate variant list file [no]
    -P             Print out basic statistics on stdout [no]
    
    """

try:
    opts, args = getopt.getopt(sys.argv[1:], "hf:i:q:c:d:o:g:a:r:s:FCMPv:z:t:")
except getopt.GetoptError, err:
    print str(err) 
    usage()
    sys.exit()

fasta_dir='/usr/local/share/genomes'
mtdna_fasta='chrRSRS.fa'
inputfile=''
hgenome_fasta='hg19RSRS.fa'
mqual=25
clev=0.80
cov=5
glen=10
basename='mtDNAassembly'
sexe='samtools'
crf=0
crc=0
crm=0
pout=0
hf=float(0.8)
tail=5
for o,a in opts:
    if o == "-h":
        usage()
        sys.exit()
    elif o == "-a": hgenome_fasta = a
    elif o == "-c": clev = float(a)
    elif o == "-d": cov = int(a)
    elif o == "-f": mtdna_fasta = a
    elif o == "-g": glen = int(a)
    elif o == "-i": inputfile = a
    elif o == "-q": mqual = int(a)
    elif o == "-o": basename = a
    elif o == "-r": fasta_dir = a
    elif o == "-s": sexe = a
    elif o == "-t":
        if int(a)<5:
            tail = 5
        else:
            tail = int(a)
    elif o == "-z": hf = float(a)
    elif o == "-v": indel_obs = float(a)
    elif o == "-F": crf = 1
    elif o == "-C": crc = 1
    elif o == "-M": crm = 1
    elif o == "-P": pout = 1
    else:
        assert False, "Unhandled option."

# DS
mtdnafile=os.path.join(fasta_dir,mtdna_fasta)
hgenome=os.path.join(fasta_dir,hgenome_fasta)
print mtdnafile
print hgenome

sample_name = os.getcwd().split(os.sep)[-1].split('_')[1]
print "assembleMTgenome for sample", sample_name

if not os.path.exists(mtdnafile):
    usage()
    sys.exit('File %s does not exist.' %(mtdnafile))

if not os.path.exists(inputfile):
    usage()
    sys.exit('File %s does not exist.' %(inputfile))    

ext=inputfile.split('.')[-1]
basext=inputfile.replace('.'+ext,'')

if ext not in ['sam','bam','pileup']:
    usage()
    sys.exit('Input file name must contain: .sam, .bam or .pileup.')

samfile=basext+'.sam'
bamfile=basext+'.bam'
pileupfile=basext+'.pileup'

if ext in ['sam','bam'] and not os.path.exists(hgenome):
    sys.exit('Human genome file does not exist.')

if ext in ['sam','bam'] and not os.path.exists(hgenome+'.fai'):
    sys.exit('Human genome indices do not exist. Run samtools faidx fist.')

r=re.compile("#+")
r1=re.compile("\^.{1}")
rr=re.compile("[\+\-]{1}[0-9]+")

def normS(s,ref):
    c=re.finditer(rr,s)
    sl=list(s)
    cc=[(x.start(),x.end()) for x in c]
    for i in cc:
        n=int(''.join(sl[i[0]+1:i[1]]))
        sl[i[0]:i[1]+n]=['#' for xx in range(len(sl[i[0]:i[1]+n]))]
    ns=''.join(sl)
    ns=ns.replace('#','')
    ss=''
    for i in ns:
        if i in '.,ACGTNacgtn<>*': ss+=i
    return (ss.replace('.',ref)).replace(',',ref)
    
def nuc(seq):
    d={'A':0,'C':0,'G':0,'T':0,'N':0}
    for i in seq:
        if d.has_key(i): d[i]+=1
        else: d['N']+=1
    return d

dn={'A':'T','T':'A','C':'G','G':'C'}

def comp(s):
    ss=''
    for i in s:
        if dn.has_key(i): ss+=dn[i]
        else: ss+='N'
    return ss

def entropy(d):
    h=0
    f={}
    for i in d.keys():
        try: v=float(d[i])/sum(d.values())
        except: v=0.0
        if v > 0:
            f[i]=v
    for k in f.keys():
        h+=f[k]*math.log10(f[k])
    return h*(-1)

def ff(v,l):
    for i in l:
        x=0
        for j in i:
            if j in v: x+=1
        if x==len(v): return i
    return 0

dIUPAC={'AG':'R','CT':'Y','GC':'S','AT':'W','GT':'K','AC':'M','CGT':'B','AGT':'D','ACT':'H','ACG':'V'}

def getIUPAC(f):
    vv=''.join([i[1] for i in f if i[0]>0])
    k=ff(vv,dIUPAC.keys())
    if k!=0: return dIUPAC[k]
    else: return '#'

def freq(d):
    f=[]
    for i in 'ACGT':
        try: v=float(d[i])/sum(d.values())
        except: v=0.0
        f.append((v,i))
    f.sort()
    f.reverse()
    if f[0][0] >= clev:
        return f[0][1]
    elif f[0][0] >= clev/2 and f[1][0] >= clev/2:
        return getIUPAC(f[0:2])
    elif f[0][0] >= clev/3 and f[1][0] >= clev/3 and f[2][0] >= clev/3:
        return getIUPAC(f[0:3])
    elif f[0][0] >= clev/4 and f[1][0] >= clev/4 and f[2][0] >= clev/4 and f[3][0] >= clev/4:
        return 'X'
    else:
        return '#'

if mtdnafile==None:
    usage()
    sys.exit('Please insert a valid mtDNA file in fasta format.')

if inputfile==None:
    usage()
    sys.exit('Please insert a valid pileup file.')

if ext=='sam':
    print 'Converting SAM to BAM...'
    cmd='%s view -bt %s.fai %s > %s' %(sexe,hgenome,samfile,bamfile)
    print '`'+cmd+'`'
    os.system(cmd)
    ext='bam'

if ext=='bam':
    print 'Sorting and indexing BAM...'
    cmd1='%s sort %s.bam %s-sorted' %(sexe,basext,basext)
    cmd2='%s index %s-sorted.bam' %(sexe,basext)
    print '`'+cmd1+'`'
    os.system(cmd1)
    print '`'+cmd2+'`'
    os.system(cmd2)
    print 'Creating pileup...'
    cmd3='%s mpileup -B -d 500000 -f %s %s-sorted.bam > %s.pileup' %(sexe,hgenome,basext,basext)
    print '`'+cmd3+'`'
    os.system(cmd3)
    
mtdna={}
x=1

print 'Reading mtDNA sequence...'

try:
    f=open(mtdnafile)
except IOError:
    sys.stderr.write('mtDNA file could not be open')
    sys.exit()

# Reads the reference mtDNA to the mtdna dict, which has the positions
# as keys. Positions start at 1, not at 0
for i in f:
    if i.strip()=='': continue
    if i.startswith('>'): continue
    for j in i.strip():
        mtdna[x]=(j.upper(),['#',(0,0,0,0),0,0.0])
        x+=1

f.close()

if len(mtdna) == 0:
    sys.stderr.write('Ouch... mtDNA was not properly read!')
    sys.exit()

print 'Reading pileup file...'

f=open(pileupfile)

for i in f:
    if i.strip()=='': continue
    l=(i.strip()).split('\t')
    if l[0]!=mtdna_fasta.split('.')[0]: continue
    pos=int(l[1])
    if pos > len(mtdna):
        sys.stderr.write('Positions of reads hanging off the end of the chromosome, consider running samtools CleanSam to soft-clip\n')
        sys.stderr.write('Position read: %s\n' % pos)
        sys.stderr.write('Length  mtdna: %i\n' % len(mtdna))
    if len(l) == 6:
        ref=l[2]
        seq=normS(re.sub(r1,"",l[4]),l[2])
        qual=l[5]
        s=''
        q=0
        for j in range(len(seq)):
            if seq[j] not in '<>*' and ord(qual[j])-33 >= mqual:
                s+=seq[j].upper()
                q+=(ord(qual[j])-33)
        try: mq=float(q)/len(s)
        except: mq=0.0
        dnuc=nuc(s)
        mfreq=freq(dnuc)
        lnuc=(dnuc['A'],dnuc['C'],dnuc['G'],dnuc['T'])
        cnuc='#'
        if len(s) >= cov:
            cnuc=mfreq
        mtdna[pos][1][0]=cnuc   # Alternative base
        mtdna[pos][1][1]=lnuc   # Base counts
        mtdna[pos][1][2]=len(s) # Coverage
        mtdna[pos][1][3]=mq     # Quality

f.close()

# The `mtdna` dict() structure has the position (coordinate) as the key (starts with 1 not 
# with 0). The value is a tuple() of 2 members:
#    [0]  Reference base
#    [1]  list() of 4 members:
#         [0]  Consensus IUPAC base (N stands for "all nucleotides" and '#' stands for unknown)
#         [1]  list() of base counts (A, C, T, G) in the reads
#         [2]  Coverage of the base
#         [3]  Average QS

print 'Assembling...'
tablefile=basename+'-table.txt'
statfile=basename+'-statistics.txt'
coveragefile=basename+'-coverage.txt'
contigfile=basename+'-contigs.fasta'
mutlistfile=basename+'-mutlist.txt'
aseq=''
f=open(tablefile,'w')
f.write('Position\tRefNuc\tConsNuc\tCov\tMeanQ\tBaseCount(A,C,G,T)\n')
assb=0
totb=0
cop=0
maxCval=1
for i in range(len(mtdna)):
    line=[str(i+1),mtdna[i+1][0],mtdna[i+1][1][0],str(mtdna[i+1][1][2]),"%.2f" %(mtdna[i+1][1][3]),str(mtdna[i+1][1][1])]
    f.write('\t'.join(line)+'\n')
    if mtdna[i+1][1][0] !='#':
        aseq+=mtdna[i+1][0]
        cop+=mtdna[i+1][1][2]
        if mtdna[i+1][1][2] > maxCval:
            maxCval=mtdna[i+1][1][2]
        assb+=1
    else:
        aseq+='#'
    totb+=1
    
f.close()

try:
    passb=(float(assb)/totb)*100
except:
    passb=0.0

try:
    covmt=(float(cop)/assb)
except:
    covmt=0.0

fseq=aseq.replace('#','N')
rseq=comp(fseq)

bcomp=nuc(fseq)
bcomp2=nuc(rseq)

try:
    pa=float(bcomp['A'])/sum(bcomp.values())
except:
    pa=0.0

try:
    pc=float(bcomp['C'])/sum(bcomp.values())
except:
    pc=0.0

try:
    pg=float(bcomp['G'])/sum(bcomp.values())
except:
    pg=0.0

try:
    pt=float(bcomp['T'])/sum(bcomp.values())
except:
    pt=0.0

try:
    pgc=float(bcomp['G']+bcomp['C'])/sum(bcomp.values())
except:
    pgc=0.0

try:
    gcskl=float(bcomp['C']-bcomp['G'])/(bcomp['C']+bcomp['G'])
except:
    gcskl=0.0

try:
    gcskh=float(bcomp2['C']-bcomp2['G'])/(bcomp2['C']+bcomp2['G'])
except:
    gcskh=0.0

gaps=[]
for i in re.finditer(r,aseq):
    cc=(i.start()+1,i.end())
    if (cc[1]-cc[0])+1 >= glen: gaps.append(cc)

contigs=[]
if len(gaps)!=0:
    for i in range(len(gaps)-1):
        cc=(gaps[i][1]+1,gaps[i+1][0]-1)
        contigs.append((cc,fseq[cc[0]-1:cc[1]]))
    if gaps[0][0]!=1:
        cc=(1,gaps[0][0]-1)
        contigs.insert(0,(cc,fseq[cc[0]-1:cc[1]]))
    if gaps[-1][1]!=len(aseq):
        cc=(gaps[-1][1]+1,len(aseq))
        contigs.append((cc,fseq[cc[0]-1:cc[1]]))
    contigs.sort()
else:
    cc=(1,len(aseq))
    contigs=[(cc,fseq[cc[0]-1:cc[1]+1])]

# print out option
if pout:
    print '============================='
    print 'Basic statistics:'
    print 'Assembled bases: %.2f' %(passb)+'%'
    print 'Mean coverage depth: %.2f' %(covmt)
    print 'Number of Contigs: %i' %(len(contigs))
    print 'Base composition [A,C,G,T]: %.2f,%.2f,%.2f,%.2f' %(pa,pc,pg,pt)
    print '============================='
    print ""
#

sam_file = open(basext+'.sam', 'r')
sam = sam_file.readlines()
sam_file.close()

mt_table_file = open(tablefile, 'r')
mt_table = mt_table_file.readlines()
mt_table_file.close()

# Calling of indels and mismatches. In the case of indels, mt_table (file that was generated in a
# previous step) is not used. However, for calling mismatches I think it is.
print " -Calling mtvcf_main_analysis..."
mut_events = mtVariantCaller.mtvcf_main_analysis(mt_table, sam, sample_name, cov, indel_obs, tail)
print " -mtvcf_main_analysis DONE"
if os.path.exists('..'+os.sep+'VCF_dict_tmp'):
    VCF_dict = ast.literal_eval(open('..'+os.sep+'VCF_dict_tmp', 'r').read()) # global VCF dict
    print "Mutation events will be appended to existing global VCF dict ../VCF_dict_tmp"
else:
    VCF_dict = {} # global VCF dict
    print "Creating new global VCF dict ../VCF_dict_tmp"

if mut_events:
    print "Updating the VCF dict..."
    VCF_dict.update(mut_events)

mut_events_cellar = open('..'+os.sep+'VCF_dict_tmp', 'w')
mut_events_cellar.write(str(VCF_dict))
mut_events_cellar.close()


if crf:
    position=1
    f=open(contigfile,'w')
    print "Generating fasta output..."
    for i in contigs:
        string_seq = i[1]
        nuc_index = i[0][0]
        dict_seq = {}
        for nuc in string_seq:
            dict_seq[nuc_index] = nuc
            nuc_index += 1
        # This only gathers consensus bases for the mut_events
        consensus_single = mtVariantCaller.get_consensus_single(mut_events[mut_events.keys()[0]],hf=hf)
        for p_info in consensus_single:
            if p_info[0] in dict_seq.keys():
                if p_info[-1] == 'mism':
                    dict_seq[p_info[0]] = p_info[1][0] # check THIS
                elif p_info[-1] == 'ins':
                    dict_seq[p_info[0]+'.1'] = p_info[1][0][1:]
                elif p_info[-1] == 'del':
                    for deleted_pos in p_info[1]:
                        if deleted_pos < len(dict_seq):
                            del(dict_seq[deleted_pos])
        # sort positions in dict_seq and join to have the sequence
        contig_seq = ''
        for j in sorted(dict_seq.keys()):
            contig_seq += dict_seq[j]
        new_i = ((i[0][0], i[0][1]), contig_seq)
        f.write('>Contig.%i|%i-%i\n' %(position,new_i[0][0],new_i[0][1]))
        for j in range(0,len(new_i[1]),60):
            f.write(new_i[1][j:j+60]+'\n')
        position+=1
    f.close()

if crc:
    f=open(coveragefile,'w')
    print "Writing coverage file..."
    dass={}
    for i in contigs:
        dass[i[0]]=[0,0,0,0,0]
    dann={(1,578):['D-Loop1',0,0],(16025,16571):['D-Loop2',0,0],(650,1603):['RNR1',0,0],(1673,3230):['RNR2',0,0],(3308,4263):['ND1',0,0],(4470,5512):['ND2',0,0],(5904,7446):['COX1',0,0],(7586,8270):['COX2',0,0],(8527,9208):['ATP6',0,0],(8366,8573):['ATP8',0,0],(9207,9991):['COX3',0,0],(10059,10405):['ND3',0,0],(10470,10767):['ND4L',0,0],(10760,12138):['ND4',0,0],(12337,14149):['ND5',0,0],(14149,14674):['ND6',0,0],(14747,15888):['CYTB',0,0]}
    for i in range(len(mtdna)):
        for j in dass:
            if j[0]<=i+1<=j[1]:
                if mtdna[i+1][1][2]>0:
                    dass[j][1]+=1
                    dass[j][2]+=mtdna[i+1][1][2]
                if mtdna[i+1][1][0] !='#':
                    dass[j][3]+=1
                    dass[j][4]+=mtdna[i+1][1][2]
                    dass[j][0]+=1
        for j in dann:
            if j[0]<=i+1<=j[1]:
                if mtdna[i+1][1][0] !='#':
                    dann[j][1]+=1
                    dann[j][2]+=mtdna[i+1][1][2]            
    x=1
    f.write('Coverage of Assembled mtDNA. Coverage %.4f - Per base depth %.3f.\n' %(passb,covmt))
    f.write('Name\tType\tStart\tEnd\tLength\tCoverage\tMeanDepth\n')
    for i in contigs:
        vv=dass[i[0]]
        #cv1=float(vv[1])/vv[0]
        #cvd1=float(vv[2])/vv[0]
        #cv2=float(vv[3])/vv[0]
        cvd2=float(vv[4])/vv[0]
        contiglen=(i[0][1]-i[0][0])+1
        fcov=(float(vv[0])/contiglen)*100
        f.write('Contig.%i\tContig\t%i\t%i\t%i\t%.3f\t%.3f\n' %(x,i[0][0],i[0][1],contiglen,fcov,cvd2))
        x+=1
    x=1
    for i in gaps:
        gaplen=(i[1]-i[0])+1
        f.write('Gap.%i\tGap\t%i\t%i\t%i\t0.000\t0.000\n' %(x,i[0],i[1],gaplen))
        x+=1
    for i in dann:
        vv=dann[i]
        flen=(i[1]-i[0])+1
        try: fcov=(float(vv[1])/flen)*100
        except: fcov=0.0
        try: cvd=float(vv[2])/vv[1]
        except: cvd=0.0
        f.write('%s\tAnnotation\t%i\t%i\t%i\t%.3f\t%.3f\n' %(vv[0],i[0],i[1],flen,fcov,cvd))
    f.close()

if crm:
    f=open(mutlistfile,'w')
    print "Writing file with list of mutations..."
    for i in mut_events[mut_events.keys()[0]]:
        if i[-1] == 'mism':
            for j in range(0,len(i[3])):
                f.write('mism\t%i\t%s\t%s\t%i\t%i\t%.4f\t%.4f\t%.4f\n' %(i[0],i[1],i[3][j],i[4][j],i[2],i[6][j],i[7][j],i[8][j])) 
        elif i[-1] == 'del':
            for j in range(0,len(i[1])):
                f.write('del\t%i\t%s\t%s\t%i\t%i\t%.4f\t%.4f\t%.4f\n' %(i[0],i[1][j],i[3][j],i[4][j],i[2],i[6][j],i[7][j],i[8][j]))
        elif i[-1] == 'ins':
            for j in range(0,len(i[3])):
                f.write('ins\t%i\t%s\t%s\t%i\t%i\t%.4f\t%.4f\t%.4f\n' %(i[0],i[1],i[3][j],i[4][j],i[2],i[6][j],i[7][j],i[8][j])) 
    f.close()

print "Finished with this ID ----\n"

