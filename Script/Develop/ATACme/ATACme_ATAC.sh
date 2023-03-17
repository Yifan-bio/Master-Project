
trim() {
    trim_galore --paired --retain_unpaired --output_dir $trim_dir $R1 $R2
}

bsbolt_align() {
    bsbolt -F1 $R1 -F2 $R2 -DB $Index -O $WDIR/${prefix}_bsbolt.bam -t 4 -OT 2
}

bsbolt_filter() {
    samtools fixmate -O bam $WDIR/${prefix}_bsbolt.bam $WDIR/${prefix}_bsbolt_fix.bam
    samtools view -@ 4 -h $WDIR/${prefix}_bsbolt_fix.bam | grep -v chrM | samtools sort -@ 4 -O bam -o $WDIR/${prefix}_bsbolt_rmChrM.bam
    samtools markdup -r -@ 4 $WDIR/${prefix}_bsbolt_rmChrM.bam $WDIR/${prefix}.rmChrM_dedup.bam
    samtools view -h -b -q 20 -@ 4 -F 3852 -f 2 -o ${prefix}.final.bam ${prefix}.rmChrM_dedup.bam
    samtools index ${prefix}.final.bam
}

Meth_call() {
    MethylDackel
}

Peak_call() {
    
}
