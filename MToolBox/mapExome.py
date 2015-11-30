#!/usr/bin/python

"""
	Written by Ernesto Picardi - e.picardi@biologia.uniba.it
	"""

import getopt, sys, os

def usage():
	print """Map FASTQ onto mtDNA
		Options:
		-a		Input Fastq
		-b 		Input Fastq for pair-end (optional)
		-c		Input Fastq for single-end
		-s		Use gmap for read alignment [default is to use gsnap]
		-g		GSNAP/GMAP executable [/usr/local/bin/gsnap]
		-D		GSNAP/GMAP database location [/usr/local/share]
		-M		GSNAP/GMAP database for mtDNA [chrRSRS]
		-H		GSNAP/GMAP database for complete human genome [hg19RSRS]
		-t		GSNAP/GMAP threads [8]
		-o		Out folder
		"""

try:
	opts, args = getopt.getopt(sys.argv[1:], "hsa:b:c:g:D:M:H:t:o:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit()

if len(opts)==0:
	usage()
	sys.exit()
fastq1=None
fastq2=None
fastq3=None
gsnapexe='/usr/local/bin/gsnap'
gsnapdb='/usr/local/share/gmapdb'
use_gsnap=True
mtdb='chrRSRS'
humandb='hg19RSRS'
mqual=30
thread=8
folder=os.path.join(os.getcwd(),'OUTfolder2')
for o,a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-a": fastq1 = a.replace("*", "")
	elif o == "-b": fastq2 = a.replace("*", "")
	elif o == "-c": fastq3 = a.replace("*", "")
	elif o == "-s": use_gsnap = False
	elif o == "-g": gsnapexe = a
	elif o == "-D": gsnapdb = a
	elif o == "-M": mtdb = a
	elif o == "-H": humandb = a
	#elif o == "-q": mqual = int(a)
	elif o == "-t": thread = int(a)
	elif o == "-o": folder = a
	else:
		assert False, "unhandled option"

print "  Input files: %s, %s, %s" % (fastq1, fastq2, fastq3)

def rev(seq):
	d={'A':'T','T':'A','C':'G','G':'C','N':'N'}
	s=''.join([d[x] for x in seq])
	return s[::-1]

if not os.path.exists(folder): os.mkdir(folder)

"""
	map1cmd='%s -D %s --gunzip -d %s -A sam --nofails --pairmax-dna=500 --query-unk-mismatch=1 -n 1 -Q -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,thread,fastq1,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	print 'Mapping onto mtDNA...'
	os.system(map1cmd)
	"""
RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample'

if use_gsnap == False:
	gsnapexe = os.path.join(os.path.dirname(gsnapexe), 'gmap')
	# If input files are compressed, they will be uncompressed and the uncompressed versions deleted after alignment
	L = [fastq1, fastq2, fastq3]
	fD = {}
	for p in L:
		if p:
			if p.endswith(".gz"):
				os.system("gzip -dc %s > %s" % (p, p[:-3]))
				fD[p] = p[:-3]
			elif p.endswith(".bz") or p.endswith(".bz2"):
				os.system("bzip2 -dc %s > %s" % (p, '.'.join(p.split(".")[:-1])))
				fD[p] = '.'.join(p.split(".")[:-1])
			else:
				fD[p] = p
	if fastq2!=None and fastq3!=None:
		map1cmd='%s -D %s -d %s -f samse --nofails %s -n 1 -O -t %i %s %s %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fD[fastq1],fD[fastq2],fD[fastq3],os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	elif fastq2!=None and fastq3== None:
		map1cmd='%s -D %s -d %s -f samse --nofails %s -n 1 -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fD[fastq1],fD[fastq2],os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	else:
		map1cmd='%s -D %s -d %s -f samse --nofails %s -n 1 -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fD[fastq1],os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	print '  Mapping onto mtDNA...'
	#print map1cmd
	os.system(map1cmd)
	for k in fD.keys(): # remove uncompressed files after alignment
		if k.endswith(".gz") or k.endswith(".bz") or k.endswith(".bz2"):
			os.remove(fD[k])
else:
	if fastq2!=None and fastq3!=None:
		map1cmd='%s -D %s --gunzip --bunzip2 -d %s -A sam --nofails --pairmax-dna=500 --query-unk-mismatch=1 %s -n 1 -Q -O -t %i %s %s %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fastq1,fastq2,fastq3,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	elif fastq2!=None and fastq3== None:
		#map1cmd='%s -D %s --gunzip --bunzip2 -d %s -A sam --nofails --pairmax-dna=500 --query-unk-mismatch=1 %s -n 1 -Q -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fastq1,fastq2,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
		map1cmd='%s -D %s -d %s -A sam --nofails --pairmax-dna=500 --query-unk-mismatch=1 %s -n 1 -Q -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fastq1,fastq2,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	else:
		map1cmd='%s -D %s --gunzip --bunzip2 -d %s -A sam --nofails --pairmax-dna=500 --query-unk-mismatch=1 %s -n 1 -Q -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,RG_tag,thread,fastq1,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
	print '  Mapping onto mtDNA...'
	#print map1cmd
	os.system(map1cmd)

#

print '  Extracting FASTQ from SAM...'
mtoutsam=os.path.join(folder,'outmt.sam')
dics={}
f=open(mtoutsam)
for i in f:
	# original version
	# if i.strip()=='': continue
	if i.strip()=='' or i.startswith('@'): continue
	l=(i.strip()).split('\t')
	if l[2]=='*': continue
	if dics.has_key(l[0]): dics[l[0]].append(l)
	else: dics[l[0]]=[l]
f.close()
single,pair1,pair2=[],[],[]

for i in dics:
	ll=dics[i]
	if len(ll)==1:
		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
		single.append(entry)
	else:
		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
		pair1.append(entry)
		strand,seq,qual=int(ll[1][1]) & 16,ll[1][9],ll[1][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[1][0],seq,'+',qual])+'\n'
		pair2.append(entry)

sig,pai=0,0
if len(single)!=0:
	mtoutfastq=os.path.join(folder,'outmt.fastq')
	out=open(mtoutfastq,'w')
	out.writelines(single)
	out.close()
	sig=1
if len(pair1)!=0:
	mtoutfastq1=os.path.join(folder,'outmt1.fastq')
	out=open(mtoutfastq1,'w')
	out.writelines(pair1)
	out.close()
	mtoutfastq2=os.path.join(folder,'outmt2.fastq')
	out=open(mtoutfastq2,'w')
	out.writelines(pair2)
	out.close()
	pai=1


if sig:
	print '  Mapping onto complete human genome... single reads'
	if use_gsnap == False:
		map2cmd='%s -D %s -d %s -f samse --nofails -n 1 -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,humandb,thread,mtoutfastq,os.path.join(folder,'outhumanS.sam'),os.path.join(folder,'loghumanS.txt'))
	else:
		map2cmd='%s -D %s -d %s -A sam --nofails --query-unk-mismatch=1 -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,humandb,thread,mtoutfastq,os.path.join(folder,'outhumanS.sam'),os.path.join(folder,'loghumanS.txt'))
	#print map2cmd
	os.system(map2cmd)
if pai:
	print '  Mapping onto complete human genome... pair reads'
	if use_gsnap == False:
		map3cmd='%s -D %s -d %s -f samse --nofails -n 1 -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,humandb,thread,mtoutfastq1,mtoutfastq2,os.path.join(folder,'outhumanP.sam'),os.path.join(folder,'logmt.txt'))
	else:
		map3cmd='%s -D %s -d %s -A sam --nofails --query-unk-mismatch=1 -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,humandb,thread,mtoutfastq1,mtoutfastq2,os.path.join(folder,'outhumanP.sam'),os.path.join(folder,'loghumanP.txt'))
	#print map3cmd
	os.system(map3cmd)

print '  Reading Results...'
if sig:
	hgoutsam=os.path.join(folder,'outhumanS.sam')
	dicsingle={}
	f=open(hgoutsam)
	for i in f:
		if i.strip()=='': continue
		l=(i.strip()).split('\t')
		if l[2]=='*': continue
		if dicsingle.has_key(l[0]):
			dicsingle[l[0]].append(l)
		else:
			dicsingle[l[0]]=[l]
	f.close()
if pai:
	hgoutsam2=os.path.join(folder,'outhumanP.sam')
	dicpair={}
	f=open(hgoutsam2)
	for i in f:
		if i.strip()=='': continue
		l=(i.strip()).split('\t')
		if l[2]=='*': continue
		if dicpair.has_key(l[0]):
			dicpair[l[0]].append(l)
		else:
			dicpair[l[0]]=[l]
	f.close()

print '  Filtering reads...'
good=[]
for i in dics:
	ll=dics[i]
	if len(ll)==1:
		if dicsingle.has_key(i):
			r=dicsingle[i]
			#print ll
			#print r
			if len(r)==1:
				if r[0][2]==ll[0][2] and ll[0][3]==r[0][3]: good.append('\t'.join(ll[0])+'\n')
		else: good.append('\t'.join(ll[0])+'\n')
	else:
		if dicpair.has_key(i):
			r=dicpair[i]
			if len(r) == 2:
				if r[0][2]==ll[0][2] and ll[0][3]==r[0][3] and r[1][2]==ll[1][2] and ll[1][3]==r[1][3]:
					good.append('\t'.join(ll[0])+'\n')
					good.append('\t'.join(ll[1])+'\n')
		else:
			good.append('\t'.join(ll[0])+'\n')
			good.append('\t'.join(ll[1])+'\n')

"""
	hgoutsam=os.path.join(folder,'outhuman.sam')
	goodhits={}
	f=open(hgoutsam)
	for i in f:
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	nn=l[0].split()[0]
	if l[2]=='chrM': mt=1
	else: mt=2
	if goodhits.has_key(nn): goodhits[nn].append(mt)
	else: goodhits[nn]=[mt]
	f.close()
	"""
finalsam=os.path.join(folder,'OUT.sam')
out=open(finalsam,'w')
# add SAM Header and Read Group for GATK Indels realignment
out.write("@SQ	SN:%s	LN:16569\n" % mtdb)
out.write("@RG	ID:sample	PL:sample	PU:sample	LB:sample	SM:sample\n")
out.writelines(good)
out.close()

"""
	f=open(mtoutsam)
	xx=1
	for i in f:
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	nn=l[0].split()[0]+'.'+str(xx)
	if not goodhits.has_key(nn): out.write(i)
	else:
	if sum(goodhits[nn])==1: out.write(i)
	xx+=1
	f.close()
	out.close()
	"""
print '  Outfile saved on %s.' %(finalsam)
print 'Success.'

