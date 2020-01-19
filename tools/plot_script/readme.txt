heatmap generator
same chr:
perl Stat_share.pl bam mapq bin chr start end flank prefix withchr cut
diff chr:
perl Stat_share_dif.pl bam mapq bin chr1 start1 end1 chr2 start2 end2 flank prefix withchr cut

mapq: filt reads with lower mapq;
bin: window to count the barcode;
flank: flank size by bp;
withchr: if input parameter (chr/chr1/chr2) is same as bam ,set 0;else will add "chr" ahead;
cut:remove the most ratio of share barcode nums(0-1),0 means not remove;


