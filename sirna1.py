import sys
import os
import re
from Bio.SeqUtils import GC
from Bio.Seq import Seq

seq = sys.argv[2]
seq = seq.replace('\n','').replace('\t','')

if len(seq) > 200:
    seq = seq[100:]
    seq = seq[0:len(seq)-99]
else:
    seq = seq

seq_len = len(seq)
index = [m.start() for m in re.finditer(r'AA',seq)]
index_list = []
for i in index:
	if i + 20 > seq_len:
		continue
	else:
		index_list.append(i)

pre_sirna = []
for index in index_list:
	pre_sirna.append(seq[index:index + 23])

pre_gc_sirna = []
index_gc_sirna = []
for i,seq1 in enumerate(pre_sirna):
	if (30 < round(GC(seq1),2) < 58) and ('AAAA' not in seq1[2:21]) and ('TTTT' not in seq1[2:21]):
		pre_gc_sirna.append(seq1)
		index_gc_sirna.append(index_list[i])

num = len(pre_gc_sirna)
num_mat = [100] * num

SS = []
for seq2 in pre_gc_sirna:
	SS.append(seq2[2:])

for i,seq2 in enumerate(SS):
	if seq2.endswith('TT'):
		num_mat[i] += 40
	
	if seq2[0:3].count('A') != seq2[19:].count('A') or seq2[0:3].count('T') != seq2[19:].count('T'):
		num_mat[i] += 1
	
	if (seq2[19:].count('A') + seq2[19:].count('T')) - (seq2[0:3].count('A') + seq2[0:3].count('T')) >= 1:
		num_mat[i] += 1
	
	if seq2[0] == 'G' or seq2[0] == 'C':
		num_mat[i] += 20
	
	if seq2[5] == 'A':
		num_mat[i] += 10
	
	if seq2[18] == 'A' or seq2[18] == 'T':
		num_mat[i] += 20
	
	if seq2[18] == 'G':
		num_mat[i] -= 5
			
	if seq2[0] == 'T':
		num_mat[i] -= 5	

score_seq = []
for i,seq3 in enumerate(SS):
	score_seq.append([num_mat[i],str(index_gc_sirna[i] + 101), 
					pre_gc_sirna[i], str(Seq(seq3).transcribe()), 
					str(Seq(pre_gc_sirna[i][0:21]).reverse_complement().transcribe())])

sort = sorted(score_seq,reverse = True)

type_dict = {'1':'rat', '2':'pig', '3':'mouse','4':'human', '5':'Golden_Hamster', '6':'Ecoli', '7':'duck',
            '8':'chicken', '9':'Angiostrongylus_cantonensis', '10':'cattle', '11':'sheep'}

if len(sort) == 0:
    print('Not Found!')
else:
    f1 = open('/var/www/html/siRNA/sense_seq_1.fa','w')
    f1.write('>1' + '\n' + sort[0][2][2:])
    f1.close()
    str1 = '/home/db/genomes/' + type_dict[sys.argv[1]] +'/' + type_dict[sys.argv[1]]
    
    os.system("sudo /home/daniu/project/software/bowtie-1.2/bowtie  -f -v 3 " + str1 + " /var/www/html/siRNA/sense_seq_1.fa /var/www/html/siRNA/sense_seq_1.log")

    with open('/var/www/html/siRNA/sense_seq_1.log') as f2:
        for line in f2:
            line_list = line.strip('\n').split('\t')
            length = 0
            if line_list[-1].split(',') == ['']:
                length = len(sort[0][2][2:])
            else:
                length = len(sort[0][2][2:]) - len(line_list[-1].split(','))
    
    total_length = len(sort[0][2][2:])
    match = str(round(float(length) / total_length,2) *100) + "%"
    
    print('blast(%)' + '@' + 'position' + '@' + 'Forward' + '@' + 'Reverse')
    print(match + '@' + str(int(sort[0][1])+2) + '-' + str(int(sort[0][1])+23) + '@' + sort[0][3] + '@' + sort[0][4])
    
    #sys.stdout = savestd
    #f1 = open('/var/www/html/siRNA/sense_seq_1.fa','w')
    f3 = open('/var/www/html/siRNA/sense_seq_2.fa','w')
    f3.write('>1' + '\n' + sort[1][2][2:])
    f3.close()
    str1 = '/home/db/genomes/' + type_dict[sys.argv[1]] +'/' + type_dict[sys.argv[1]]
    
    os.system("sudo /home/daniu/project/software/bowtie-1.2/bowtie  -f -v 3 " + str1 + " /var/www/html/siRNA/sense_seq_2.fa /var/www/html/siRNA/sense_seq_2.log")

    with open('/var/www/html/siRNA/sense_seq_2.log') as f2:
        for line in f2:
            line_list = line.strip('\n').split('\t')
            length = 0
            if line_list[-1].split(',') == ['']:
                length = len(sort[1][2][2:])
            else:
                length = len(sort[1][2][2:]) - len(line_list[-1].split(','))
    
    total_length = len(sort[1][2][2:])
    match = str(round(float(length) / total_length,2) *100) + "%"
    
    print('blast(%)' + '@' + 'position' + '@' + 'Forward' + '@' + 'Reverse')
    print(match + '@' + str(int(sort[1][1])+2) + '-' + str(int(sort[1][1])+23) + '@' + sort[1][3] + '@' + sort[1][4])

    f4 = open('/var/www/html/siRNA/sense_seq_3.fa','w')
    f4.write('>1' + '\n' + sort[2][2][2:])
    f4.close()
    str1 = '/home/db/genomes/' + type_dict[sys.argv[1]] +'/' + type_dict[sys.argv[1]]
    
    os.system("sudo /home/daniu/project/software/bowtie-1.2/bowtie  -f -v 3 " + str1 + " /var/www/html/siRNA/sense_seq_3.fa /var/www/html/siRNA/sense_seq_3.log")

    with open('/var/www/html/siRNA/sense_seq_3.log') as f2:
        for line in f2:
            line_list = line.strip('\n').split('\t')
            length = 0
            if line_list[-1].split(',') == ['']:
                length = len(sort[2][2][2:])
            else:
                length = len(sort[2][2][2:]) - len(line_list[-1].split(','))
    
    total_length = len(sort[2][2][2:])
    match = str(round(float(length) / total_length,2) *100) + "%"
    
    print('blast(%)' + '@' + 'position' + '@' + 'Forward' + '@' + 'Reverse')
    print(match + '@' + str(int(sort[2][1])+2) + '-' + str(int(sort[2][1])+23) + '@' + sort[2][3] + '@' + sort[2][4])
          
          	      		
    



