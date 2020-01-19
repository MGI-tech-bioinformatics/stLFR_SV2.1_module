Usage: perl make_Haplotype.pl bam mapq bin chr1 bp1 chr2 bp2 flank hapcut_out phase_out ratio vcf

bam: stLFR bam after rmdup;
mapq: mapping quality cut-off;
bin: bin size for extract share barcode,bp;
chr1 bp1: first break point position;
chr2 bp2: second break point position;
flank: flank length for both sides;
hapcut_out: path for hapcut2 raw result;
phase_out: path for barcode phased result(the same as SV2.0 pipeline input files);
ratio: ratio cut_off to judge the SV occor on which Haplotype
(for example 0.75 means: if over 75% barcode support Haplotype 1,the SV should on Haplotype 1, otherwise, if 25% ~ 75% barcode support Haplotype 1, the Haplotype can't be judged);
vcf:vcf file for extract homozygous mutation

out files:
Haplotype_withoutsv_1.final.list
Haplotype_withoutsv_2.final.list
Haplotype_withsv_1.final.list
Haplotype_withsv_2.final.list



