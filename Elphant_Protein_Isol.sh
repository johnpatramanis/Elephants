#!/bin/bash

#########################
#########################
HOME_DIR='/home/rjt939/Elephants/';

GENE='COL1A2';

SCAFOLD="scaffold_5";
START='55153479';
STOP='55189012';

REFERENCE_GENOME='/home/rjt939/Elephants/REFERENCE_GENOME/Loxodonta_africana.loxAfr3.dna.toplevel.fa';

IE_TABLE='';

cd $HOME_DIR

samtools faidx $REFERENCE_GENOME;

samtools faidx $REFERENCE_GENOME $SCAFOLD > 'REFERENCE_'$SCAFOLD'.fa'; ## Creat sub_file with just scaffol of iterest

cat Galaxy1-[imported__mammoth_SNPs].gd_snp | grep -P  "^""$SCAFOLD""\t">$SCAFOLD".snp" # ATLERNATIVE?

R -e '

args = commandArgs(trailingOnly=TRUE)

library(ShortRead)

directory=args[1]
setwd(directory)


fa<-readFasta(args[2])
sam=args[3]
gene=args[3]
START=as.numeric(args[4])
STOP=as.numeric(args[5])
SCAF=args[6]




VAR_TABLE=read.table(paste0(SCAF,".snp")) #load sub-table
POS=as.numeric(VAR_TABLE[,2])

REF_SEQ=as.character(DNAString(as.character(sread(fa))))


if (STOP>=START){
    VARIANTS=which(POS>=START & POS<=STOP)
    }else{
    VARIANTS=which(POS>=STOP & POS<=START)
   }
VARIANTS=VAR_TABLE[VARIANTS,]


NAMES=c("Elephas_maximus1","Elephas_maximus2","Elephas_maximus3","Mammuthus_primigenius1","Mammuthus_primigenius2")
SEQUENCES=rep(REF_SEQ,length(NAMES))


for (V in 1:dim(VARIANTS)[1]){

    VAR=VARIANTS[V,]
	REF=as.character(VAR[3][1,])
	ALT=as.character(VAR[4][1,])
	LOC=as.numeric(VAR[2])
	GENOTYPES=c(VAR[9],VAR[13],VAR[17],VAR[29],VAR[33])  #make sure they are in correct order, correspond to the the names/species
    
	# print(c(REF,ALT,LOC,V))
    print(V)



    for (i in 1:length(SEQUENCES)){
		
		GENOTYPE=GENOTYPES[[i]]

		if (GENOTYPE==0){ #homozygous for alternative?
		    print("ALTERNATIVE")
            substr(SEQUENCES[i], LOC, LOC) <- ALT #replace with variant for that sample
            #substr(SEQUENCES[i], LOC-10, LOC+10)       ##check if positions changed
            
		}else{
            print("REFERENCE")
		}

 
	}
}


#### NOW ISOLATE CHUNK OF INTEREST AND PRINT (prbly another loop)
seq<-as.character(DNAString(as.character(sread(fa)))[START:STOP]) ## Isolate location of gene in chromosome, start-stop


newseq<-ShortRead(sread=DNAStringSet(seq), id=BStringSet(paste0(gene,"_Loxodonta_africana"))) # prepare fasta sequence

writeFasta(newseq, paste0(gene,"_Loxodonta_africana",".fa")) # write it out
writeFasta(newseq, paste0(gene,"_Elephas_maximus",".fa")) # write it out



################################################################################################################
################################################################################################################
#SPLICING


' --args $HOME_DIR 'REFERENCE_'$SCAFOLD'.fa' $GENE $START $STOP $SCAFOLD



##### Change code to do add sample name!
## Sample name should be second!




#Clean any leftover files from previous runs!
rm *spliced.* 
rm *translated.fa

R -e '

################################################################################################################
################################################################################################################
#SPLICING


args = commandArgs(trailingOnly=TRUE)

library(ShortRead)

directory=args[1]
gene=args[3]
setwd(directory)


fas<-dir(pattern=paste0("^",gene,".*\\.fa$"))  #Grab all fasta files with the gene name
info<-read.table("./starts.txt", h=T, as.is=T) # name of each gene/ where each gene starts / which strand




for(x in 1:length(fas)){

    samp<-strsplit(gsub(".fa", "", fas[x]), "_")[[1]][3]
    pop<-strsplit(gsub(".fa", "", fas[x]), "_")[[1]][2]

    fa<-readFasta(fas[x]) #readFasta loads the fasta file into R
    fa<-DNAString(as.character(sread(fa))) #Turn iNAString
    tab<-read.table(paste0("./EIT/", gene, ".txt"), as.is=T, sep="\t") # Exon / Intron file

    if(info[info[,1]==gene,3]=="+"){ # If the gene of the fasta file is on the (+) strand
        starts<-info[info[,1]==gene,2] #grab first start position
        ends<-NULL
        for(i in 1:length(tab[,1])){   #use intron position to create chunks of introns
            ends<-c(ends, ((starts[length(starts)]+tab[i,2])-1)) # get end position from last start position + length of intron/exon -1
            starts<-c(starts, (starts[length(starts)]+tab[i,2])) # next start position is last end position (+1)
        }
        starts<-starts[1:length(tab[,1])] # remove the last start position
        
        starts<-starts[-grep("Intron", tab[,1])] # remove intron starts from list
        ends<-ends[-grep("Intron", tab[,1])]    #remove intron ends from list
        
    }else{ # if gene of the fasta file is on (-) strand, same as above but the reverse logic (move from right to left)
        ends<-info[info[,1]==gene,2] # what previously would be start is here the end
        starts<-NULL #the same logic as above
        for(i in 1:length(tab[,1])){ #the same logic as in the above loop, but we are moving to the left, so starts are bigger numbers than their ends
            starts<-c(starts, ends[length(ends)]-tab[i,2]+1)
            ends<-c(ends, ends[length(ends)]-tab[i,2])
        }
        
        ends<-ends[1:length(tab[,1])]   #again remove last part
        #remove introns
        starts<-starts[-grep("Intron", tab[,1])]
        ends<-ends[-grep("Intron", tab[,1])]
        starts<-rev(starts) #flip them into canonical order, from smaller number to larger
        ends <-rev(ends)
    }

    seqs<-NULL

    #now we use those starts / ends pairs to isolate the exons only from the sequence

    for(i in 1:length(starts)){
        seqs<-paste(sep="", seqs, as.character(fa[starts[i]:ends[i]]))
    }
    
    seqs<-gsub(" ", "", seqs)

    ids<-paste0(gene, "_spliced")

    newseq<-ShortRead(sread=DNAStringSet(seqs), id=BStringSet(ids))

    writeFasta(newseq, paste0(gene,"_",pop,"_",samp,"_spliced.fa"))


    }


' --args $HOME_DIR 'REFERENCE_'$SCAFOLD'.fa' $GENE $START $STOP


#######################################################################################
#######################################################################################
##BLASTING

for i in *_spliced.fa; do /home/rjt939/BLAST/ncbi-blast-2.6.0+/bin/makeblastdb -dbtype nucl -in $i; done #creates a BLAST database for each file!




rm *.blast

ls *spliced.fa |cut -f 2,3 -d "_" |sort |uniq>SAMPLES #save the names/populations of the files/samples #1KG_samples



cat SAMPLES |while read sams;
    do /home/rjt939/BLAST/blast-2.2.26/bin/blastall  -p tblastn -i ./REFERENCE_PROTEINS/$GENE".fa"  -d $GENE"_"$sams"_spliced.fa" -o $GENE"_"$sams"_spliced.blast" -F F -E 32767 -G 32767 -n T -m 0 -M PAM70;
        
    done;






R -e '

################################################################################################################
################################################################################################################
#SPLICING


args = commandArgs(trailingOnly=TRUE)

library(ShortRead)

directory=args[1]
setwd(directory)
gene=args[3]


f<-dir(pattern=paste0("^",gene,".*\\_spliced.blast$")) #load all blast files for this gene




for(i in 1:length(f)){ #for each blast file

    samp<-strsplit(gsub(".fa", "", f[i]), "_")[[1]][3]
    pop<-strsplit(gsub(".fa", "", f[i]), "_")[[1]][2]
    fout<-paste0(gene,"_",pop,"_",samp, "_translated.fa")  #name of output fasta (protein)


    blout<-f[i]   #grab each file
    zz<-pipe(paste0("grep -e \"Identities\" -e \"Query\" -e \"Sbjct\" -e \"Length of query\" ", blout, " | grep -v \"Query=\""))  #create pipe connection object for file
    a<-scan(zz, what="", sep="\n") #scan pipe ?
    close(zz) #
    zz<-pipe(paste0("grep \"letters)\" ", blout)) #new pipe
    len<-scan(zz, what="", sep="\n")  #sequence
    close(zz)
    len<-as.numeric(strsplit(strsplit(len, "\\(")[[1]][2], " ")[[1]][1]) #length of sequence
    separator<-grep("Identities", a) #

    if(length(separator)>1){
        a<-a[(separator[1]+1):(separator[2]-1)]
        b<-a[grep("Query",a)]
        tot<-as.numeric(strsplit(b[length(b)], " ")[[1]][length(strsplit(b[length(b)], " ")[[1]])])
        b<-as.numeric(strsplit(b[1], " ")[[1]][2])
        a<-a[grep("Sbjct",a)]
        a<-gsub("Sbjct: (\\d+)( *+)", "", a, perl=TRUE)
        a<-gsub(" (\\d+)", "", a, perl=TRUE)
        if(b!=1){
            a<-c(paste(collapse="", rep("X", b-1)), a)
        }
        if(tot<len){
            a<-c(a, paste(collapse="", rep("X", len-tot)))
        }
        a<-paste(collapse="", a)
        newseq<-AAStringSet(a)
        names(newseq)<-gene
        writeXStringSet(newseq, fout)
    
    }else{#if separator=0
        if(length(a)==1|length(a)==0){ #if no results at all! # seems to be called only for samples that lack AMELY, good!
            a=a[1] 
            print(paste0("No Results for 1 BLAST ",gene))
            # b<-a[grep("Query",a,ignore.case=TRUE)]
            # tot<-as.numeric(strsplit(b[length(b)], " ")[[1]][length(strsplit(b[length(b)], " ")[[1]])])
            # a=rep("X", tot)
            # a<-paste(collapse="", a)            
            # newseq<-AAStringSet(a)
            # names(newseq)<-gene
            # writeXStringSet(newseq, fout)
            
            
            }else{
            a<-a[(separator[1]+1):(length(a))]
            b<-a[grep("Query",a,ignore.case=TRUE)]
            tot<-as.numeric(strsplit(b[length(b)], " ")[[1]][length(strsplit(b[length(b)], " ")[[1]])])
            b<-as.numeric(strsplit(b[1], " ")[[1]][2])
            a<-a[grep("Sbjct",a)]
            a<-gsub("Sbjct: (\\d+)( *+)", "", a, perl=TRUE)
            a<-gsub(" (\\d+)", "", a, perl=TRUE)
            if(b!=1){
                a<-c(paste(collapse="", rep("X", b-1)), a)
            }
            if(tot<len){
                a<-c(a, paste(collapse="", rep("X", len-tot)))
            }
            a<-paste(collapse="", a)
            newseq<-AAStringSet(a)
            names(newseq)<-gene
            writeXStringSet(newseq, fout)
            }
    }
}





' --args $HOME_DIR 'REFERENCE_'$SCAFOLD'.fa' $GENE $START $STOP

