

PROT_LIST=['COL1A1', 'COL1A2', 'AHSG', 'CHAD', 'SERPINF1', 'BGN', 'THBS1', 'LUM', 'SPP1','DCN', 'MMP2', 'VTN', 'OMD', 'POSTN', 'PCOLCE', 'APOA4', 'SERPINC1', 'APOE', 'COL5A1', 'COL5A2', 'COL17A1', 'DSPP', 'F2', 'COML5', 'SPOCK1', 'ALB']

FILE=open('CODING_SNPS.gd_sap','r')
VARIANTS={}

for line in FILE:
    print(line)
    line=line.strip().split()
    if line[3] not in VARIANTS:
        VARIANTS[line[3]]=0
    if line[3] in VARIANTS:
        VARIANTS[line[3]]+=1


for J in PROT_LIST:
    if J in VARIANTS:
        print(J,VARIANTS[J])
    if J not in VARIANTS:
        print(J,0)
        
OUT=open('VARIANT_PROT_COUNTS','w')

for J in PROT_LIST:
    if J in VARIANTS:
        OUT.write('{}\t{}\n'.format(J,VARIANTS[J]))
    if J not in VARIANTS:
        OUT.write('{}\t{}\n'.format(J,0))