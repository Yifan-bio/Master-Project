#2. Extract the TPM values from the Salmon output
multipleFieldSelection.py -i ./salmon/Input/*/quant.sf -k 1 -f 4 -o ./salmon/Input/iso_tpm.txt

#3. Before running SUPPA, we need to calculate the AS events on the hg19 annotation
	#3.1: Generate the events: 
	suppa.py generateEvents -i ./salmon/Input/gencode.v40.primary_assembly.annotation.gtf -o ./salmon/Input/gencodev40_events -e SE SS MX RI FL -f ioe
	#3.2: Put all the ioe events in the same file:
	awk '
	    FNR==1 && NR!=1 { while (/^seqname/) getline; }
	    1 {print}
	' ./*ioe > ./hg19_ensembl_events_all_events.ioe

#4. Run SUPPA for getting the psi values of the events:
suppa.py psiPerEvent -i ~/annotation/hg19_ensembl_events_all_events.ioe -e ~/tra2/Salmon/quantification/iso_tpm.txt -o ~/tra2/SUPPA/events

#5. Run SUPPA for obtaining the Differential splicing analysis
	#5.1: Split the PSI and TPM files between the 2 conditions
	~/scripts/split_file.R ~/tra2/Salmon/quantification/iso_tpm.txt SRR1513329,SRR1513330,SRR1513331 SRR1513332,SRR1513333,SRR1513334 ~/tra2/Salmon/quantification/CTRL.tpm ~/tra2/Salmon/quantification/KD.tpm 
	~/scripts/split_file.R ~/tra2/SUPPA/events.psi SRR1513329,SRR1513330,SRR1513331 SRR1513332,SRR1513333,SRR1513334 ~/tra2/SUPPA/CTRL.psi ~/tra2/SUPPA/KD.psi
	#5.2: Run SUPPA
	python ~/suppa.py diffSplice -m empirical -i ~/annotation/hg19_ensembl_events_all_events.ioe -e ~/tra2/Salmon/quantification/CTRL.tpm ~/tra2/Salmon/quantification/KD.tpm -p ~/tra2/SUPPA/CTRL.psi ~/tra2/SUPPA/KD.psi -o ~/tra2/SUPPA/CTRL_KD
