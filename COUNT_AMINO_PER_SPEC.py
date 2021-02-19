FILE=open('AMINO_PER_IND.txt','r')


FILE.readline()
UNIQUES={}
ALL_ACIDS=[]
for LINE in FILE:
    line=LINE.strip().split()
    print(line)
    ALL_ACIDS.append([line[2],line[5:11]])
    
print(ALL_ACIDS)

for J in ALL_ACIDS:
    if J[0] in UNIQUES:
        pass
    if J[0] not in UNIQUES:
        UNIQUES[J[0]]=[0,0,0]  #AFR,ASIAN,MAMM




for J in ALL_ACIDS:
    PROT=J[0]
    GENOS=J[1]
    # print(GENOS)
    #MAMMOTH
    if ( ( GENOS[5]!=GENOS[0] and GENOS[5]!=GENOS[1] and GENOS[5]!=GENOS[2] and GENOS[5]!=GENOS[3] ) or   ( GENOS[4]!=GENOS[0] and GENOS[4]!=GENOS[1] and GENOS[4]!=GENOS[2] and GENOS[4]!=GENOS[3] ) ):
        UNIQUES[PROT][2]+=1
    #ASIAN
    if ( ( GENOS[1]!=GENOS[0] and GENOS[1]!=GENOS[4] and GENOS[1]!=GENOS[5] ) or   ( GENOS[2]!=GENOS[0] and GENOS[2]!=GENOS[4] and GENOS[2]!=GENOS[5] ) or ( GENOS[3]!=GENOS[0] and GENOS[3]!=GENOS[4] and GENOS[3]!=GENOS[5] ) ):
        UNIQUES[PROT][1]+=1
    #AFRICAN
    if ( GENOS[0]!=GENOS[1] and GENOS[0]!=GENOS[2] and GENOS[0]!=GENOS[3] and GENOS[0]!=GENOS[4] and GENOS[0]!=GENOS[5] ):
        UNIQUES[PROT][0]+=1



UNIQUE_TABLE=open('UNIQUE_PER_SPECIES.txt','w')

UNIQUE_TABLE.write('PROTEIN\tUNIQUE_AFRICAN\tUNIQUE_ASIAN\tUNIQUE_MAMMOTH\n')
for J in UNIQUES.keys():
    UNIQUE_TABLE.write('{}\t{}\t{}\t{}\n'.format(J,UNIQUES[J][0],UNIQUES[J][1],UNIQUES[J][2]))
    
########################################################################################################################################################################################################


#SAPS vs Asian
VSASIAN={}
for J in ALL_ACIDS:
    if J[0] in VSASIAN:
        pass
    if J[0] not in VSASIAN:
        VSASIAN[J[0]]=[0,0]  #AFR,MAMM

print(VSASIAN)
    
for J in ALL_ACIDS:
    PROT=J[0]
    GENOS=J[1]
    # print(GENOS)
    if (((GENOS[1]=='ALT') or (GENOS[2]=='ALT') or (GENOS[3]=='ALT') ) & ( GENOS[4]=='ALT' or GENOS[5]=='ALT')):
        VSASIAN[PROT][0]+=1
        
    
    if ( ( (GENOS[1]=='ALT') or (GENOS[2]=='ALT') or (GENOS[3]=='ALT') ) & ( (GENOS[4]=='REF') & (GENOS[5]=='REF') ) ):
        VSASIAN[PROT][1]+=1
        VSASIAN[PROT][0]+=1
    
    if ( ((GENOS[1]=='REF') & (GENOS[2]=='REF') & (GENOS[3]=='REF')) & ( (GENOS[4]=='ALT') or (GENOS[5]=='ALT')) ):
        VSASIAN[PROT][1]+=1
    
    if ( ((GENOS[1]=='REF') & (GENOS[2]=='REF') & (GENOS[3]=='REF')) & ( (GENOS[4]=='REF') & (GENOS[5]=='REF') ) ):
        pass
    


A_M_DIFF=open('ASIAN_MAMMOTH_DIFFS.txt','w')
A_M_DIFF.write('DIFF_ASIAN_AFRICAN\tDIFF_ASIAN_MAMMOTH\n')

for J in VSASIAN.keys():
    A_M_DIFF.write('{}\t{}\t{}\n'.format(J,VSASIAN[J][0],VSASIAN[J][1]))
    










################################################################################################################################################################
#SAPS vs Asian






################################################################################################################################################################
#SAps vs Mammoth


