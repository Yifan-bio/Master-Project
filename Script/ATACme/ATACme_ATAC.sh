
trim() {
    trim_galore --paired --retain_unpaired --output_dir $trim_dir $R1 $R2
}

trim_ATACmix() {
    cutadapt --length $ATAC_size -o $trim_dir/${prefix}_R1.fq.gz -p $trim_dir/${prefix}_R2.fq.gz $R1 $R2
}

bsbolt_align() {
    bsbolt -F1 $R1 -F2 $R2 -DB $Index -O $WDIR/${prefix}_bsbolt.bam -t 4 -OT 2
}

bsbolt_filter() {
    samtools fixmate -O bam $WDIR/${prefix}_bsbolt.bam $WDIR/${prefix}_bsbolt_fix.bam
    samtools view -@ 4 -h $WDIR/${prefix}_bsbolt_fix.bam | grep -v chrM | samtools sort -@ 4 -O bam -o $WDIR/${prefix}_bsbolt_rmChrM.bam
    picard MarkDuplicates I=$WDIR/${prefix}_bsbolt_rmChrM.bam O=$WDIR/${prefix}_bsbolt_rmChrM_rmDup.bam M=$WDIR/${prefix}_bsbolt_rmChrM_rmDup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    samtools view -h -b -q 30 -@ 4 -F 1804 -f 2 -o ${prefix}.final.bam ${prefix}.rmChrM_dedup.bam
    samtools index ${prefix}.final.bam
}

