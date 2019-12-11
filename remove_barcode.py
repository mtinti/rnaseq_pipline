
#This Python script requires Biopython 1.51 or later
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
import gzip
from Bio.Seq import Seq
import pandas as pd
import sys
import numpy as np
import tqdm
import shutil
import os

r_barcode = 'TCGCGAGGC'
f_barcode ='GCCTCGCGA'
def get_len(in_file):
	def blocks(files, size=65536):
		while True:
			b = files.read(size)
			if not b: break
			yield b
	
	with gzip.open(in_file,"rU") as f:
		res = sum(bl.count("\n") for bl in blocks(f))
		res = float(res)/3   
		return res

def find_barcode(b_seq = f_barcode, in_seq=''):
	f_in_seq = Seq(in_seq)
	r_in_seq = f_in_seq.reverse_complement()
	if b_seq in f_in_seq:
         return 1,'f'
	elif b_seq in r_in_seq:
         return 1,'r'
	else:
         return 0,'n'

def remove_barcode(seq, quol, barcode_orient):
    if barcode_orient == 'f':
        start = seq.find(f_barcode)+len(f_barcode)
        seq = seq[start:]
        quol = quol[start:]
    if barcode_orient == 'r':
        start = seq.find(r_barcode)
        seq = seq[:start]
        quol = quol[:start]
    return   seq,  quol  
         
         
def parse(file_name, inFileLength):
#Setup variables (could parse command line args instead)
	file_f = file_name+"f1.fastq.gz"
	file_r = file_name+"f2.fastq.gz"
	new_f_barcode = open(file_name+"1_barcode.fastq",'w')
	new_r_barcode = open(file_name+"2_barcode.fastq",'w') 
	new_f_removed = open(file_name+"1_removed.fastq",'w')
	new_r_removed = open(file_name+"2_removed.fastq",'w') 
	#len_file = get_len(file_f)
	#print len_file, 'records'       
	ids = []
	counts_f = []
	counts_r = []
	code_f = []
	code_r = []  
	score_f = []
	score_r = []
	len_f = []
	len_r = []               
	count = 0
	count_both = 0
	count_deleted = 0
	count_kept = 0
	f_iter = FastqGeneralIterator(gzip.open(file_f,"rt"))
	r_iter = FastqGeneralIterator(gzip.open(file_r,"rt"))
	#r1_iter = FastqGeneralIterator(gzip.open(file_r1,"rt"))
	#r2_iter = FastqGeneralIterator(gzip.open(file_r2,"rt"))
	for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in tqdm.tqdm(zip(f_iter, r_iter), total=int(inFileLength), miniters=500000):
		assert f_id.split('/')[0] == r_id.split('/')[0]
		count += 1
		#print (count)
		cf, codef =   find_barcode(in_seq=f_seq) 
		cr, coder =   find_barcode(in_seq=r_seq)
		#counts_f.append(cf)
		#code_f.append(codef)
  
		#counts_r.append(cr)
		#code_r.append(coder)  
		#ids.append(f_id.split('/')[0])
		#f_q_average = np.mean([ord(n) for n in f_q])
		#r_q_average = np.mean([ord(n) for n in r_q])
		#score_f.append(f_q_average)
		#score_r.append(r_q_average)   

		if cf == 0 and cr == 0:
			new_f_removed.write("@%s\n%s\n+\n%s\n" % (f_id, f_seq, f_q))
			new_r_removed.write("@%s\n%s\n+\n%s\n" % (r_id, r_seq, r_q))
			count_deleted+=1
		else:
			new_f_barcode.write("@%s\n%s\n+\n%s\n" % (f_id, f_seq, f_q))
			new_r_barcode.write("@%s\n%s\n+\n%s\n" % (r_id, r_seq, r_q))
			count_kept+=1
		#elif cf == 0 and cr == 1:
			#new_f_barcode.write("@%s\n%s\n+\n%s\n" % (f_id, f_seq, f_q))
			#r_seq, r_q = remove_barcode(r_seq, r_q, coder)
			#new_r_barcode.write("@%s\n%s\n+\n%s\n" % (r_id, r_seq, r_q))
   
		#elif cf == 1 and cr == 0:
			#f_seq, f_q = remove_barcode(f_seq, f_q, codef)
			#new_f_barcode.write("@%s\n%s\n+\n%s\n" % (f_id, f_seq, f_q))   
			#new_r_barcode.write("@%s\n%s\n+\n%s\n" % (r_id, r_seq, r_q))

		#else:
			#f_seq, f_q = remove_barcode(f_seq, f_q, codef)
			#r_seq, r_q = remove_barcode(r_seq, r_q, coder)      
			#new_f_barcode.write("@%s\n%s\n+\n%s\n" % (f_id, f_seq, f_q))
			#new_r_barcode.write("@%s\n%s\n+\n%s\n" % (r_id, r_seq, r_q))
			#count_both+=1
   
		#len_f.append(len(f_seq))
		#len_r.append(len(r_seq))
		
		
	new_f_barcode.close()
	new_r_barcode.close()
	new_f_removed.close()
	new_r_removed.close()                            
		
	def gzip_fastq(infile):
		with open(infile, 'rb') as f_in:
			with gzip.open(infile+'.gz', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
				os.remove(infile)

	gzip_fastq(new_f_barcode.name)
	gzip_fastq(new_r_barcode.name)
	gzip_fastq(new_f_removed.name)
	gzip_fastq(new_r_removed.name)

	print ('all:',count, 'kept:', count_kept, 'deleted:', count_deleted)
	#res=pd.DataFrame()
	#res['ids']=ids
	#res['counts_f']=counts_f
	#res['counts_r']=counts_r
	#res.set_index('ids',inplace=True)
	#res['sum']=res.sum(axis=1)
	#res['dif']=res['counts_f']-res['counts_r']

	#report = open(file_name+'barcode_report.txt','w')
	#report.write('all,f,r,b\n')
	#report.write( str(count)+','+str(np.sum(counts_f))+','+str(np.sum(counts_r))+','+str(np.sum(count_both)) )
	#report.write('\n\n')
	#report.write(str(res['dif'].value_counts()))
	#report.close()
	#res['code_f']=code_f
	#res['code_r']=code_r 
	#res['score_f']=score_f
	#res['score_r']=score_r
	#res['len_f']=len_f
	#res['len_r']=len_r

	#res.to_csv(file_name+'df.csv.gz',compression ='gzip')
 	#print res.describe()
	#print res.head()
	#del res

if __name__ == '__main__':
    
	
    parse(sys.argv[1], sys.argv[2])
    #parse('14256_1#1_')

