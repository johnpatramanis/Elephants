TABLE1=open('Galaxy1-[imported__mammoth_SNPs].gd_snp','r')
TABLE2=open('Galaxy10-[imported__mammoth_coding_SNPs].gd_sap','r')

PROT_FILE=open('Ivory proteins.txt','r')


print(TABLE1,TABLE2,PROT_FILE)


OUT1=open('AMINO_PER_IND.txt','w')

PROTS_OF_INT=[]
for LINE in PROT_FILE:
    LINE=LINE.strip()
    PROTS_OF_INT.append(LINE)


print(PROTS_OF_INT)




POS_OF_INT={}

TABLE2.readline()

for LINE in TABLE2:
    line=LINE.strip().split()
    if (line[3] in PROTS_OF_INT) and (line[4]!=line[6]):
        POS_OF_INT['_'.join([line[0],line[1]])]=[line[3],line[4],line[6]]
    
print(POS_OF_INT)




TABLE1.readline()
TABLE1.readline()

for line in TABLE1:
    line=line.strip().split()
    LOC='_'.join([line[0],line[1]])
    # print(POS_OF_INT,LOC)
    if LOC in POS_OF_INT:
        POS_OF_INT[LOC].append(['2',line[8],line[12],line[16],line[28],line[32]])
        print(POS_OF_INT[LOC])
        

SAMPS=['AFRICAN','ASIAN','ASIAN','ASIAN','MAMMOTH','MAMMOTH']

OUT1.write('SCAFFOLD\tPOSITION\tPROT\tREF_AMIN\tALT_AMIN\t')
OUT1.write('\t'.join(SAMPS))
OUT1.write('\n')
for J in POS_OF_INT.keys():
    MYLIST=[]
    for K in POS_OF_INT[J][3]:
        if K=='2' or K=='-1':
            MYLIST.append('REF')
        if K=='0' or K=='1':
            MYLIST.append('ALT')
        
    OUT1.write('{}\t{}\t{}\t{}\t{}\t'.format(J.split('_')[1],J.split('_')[2],POS_OF_INT[J][0],POS_OF_INT[J][1],POS_OF_INT[J][2]))
    OUT1.write('\t'.join(MYLIST))
    OUT1.write('\n')