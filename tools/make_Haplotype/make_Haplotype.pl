#! /usr /bin/perl -w
use strict;
die "Usage: perl $0 bam mapq bin chr1 bp1 chr2 bp2 flank hapcut_out phase_out ratio vcf\n" if @ARGV <12;

my ($bam,$mapq,$bin,$chr1,$bp1,$chr2,$bp2,$flank,$phase_out,$sv_out,$cut_off,$vcf)=@ARGV;
my (%hash1,%hash2);
my($s1,$e1,$s2,$e2);
$s1=$bp1-$flank;
if($s1<0){
	$s1=0;
}
$e1=$bp1+$flank;
$s2=$bp2-$flank;
if($s2<0){
	$s2=0;
}
$e2=$bp2+$flank;
open IN,"samtools view -F 1024 -q $mapq $bam $chr1:$s1-$e1|"or die $!;
while(<IN>){
	chomp;
	my @array=split;
	my $bar=(split/#/,$array[0])[-1];
	next if ($bar eq "0_0_0");
	my $num=int($array[3]/$bin)*$bin;
	if(exists $hash1{$num}){
		push @{$hash1{$num}},$bar;
	}else{
		$hash1{$num}=[];
		push @{$hash1{$num}},$bar;
	}
}
close IN;
print "load $chr1 finish\n";
open IN2,"samtools view -F 1024 -q $mapq $bam $chr2:$s2-$e2|"or die $!;
while(<IN2>){
	chomp;
	my @array=split;
	my $bar=(split/#/,$array[0])[-1];
	next if ($bar eq "0_0_0");
	my $num=int($array[3]/$bin)*$bin;
	if(exists $hash2{$num}){
		push @{$hash2{$num}},$bar;
	}else{
		$hash2{$num}=[];
		push @{$hash2{$num}},$bar;
	}
}
close IN2;

print "load $chr2 finish\n";

for my $k ( sort {$a <=>$b} keys %hash1){
	my %uniq;
	@uniq{@{$hash1{$k}}} = ( );
	my @uniq_array=sort keys %uniq;
	$hash1{$k}=\@uniq_array;
}


for my $j ( sort {$a <=>$b} keys %hash2){
	my %uniq;
	@uniq{@{$hash2{$j}}} = ( );
	my @uniq_array=sort keys %uniq;
	$hash2{$j}=\@uniq_array;
}
my %shareb;
for my $kk ( sort {$a <=>$b} keys %hash1){
	for my $jj ( sort {$a <=>$b} keys %hash2){
		for my $bar( @{$hash2{$jj}}){
			if ((grep m/$bar/,@{$hash1{$kk}})>0){
				$shareb{$bar}+=1;
			}
		}
	}
}
print "extract share barcodes finish\n";


my ($mc,$pc);
$mc=$pc=0;
open BP1,"<$sv_out/$chr1.barcode.phase" or die $!;
while(<BP1>){
	chomp;
	my @l=split;
	next unless exists $shareb{$l[0]};
	if($l[-1] eq "PAT"){
		$pc+=1;
	}elsif($l[-1] eq "MAT"){
		$mc+=1;
	}
}
close BP1;

my($type1,$type2,$type3,$type4);
if($pc/($pc+$mc)>=$cut_off){
	$type1=1;
	$type3=2;
}elsif($pc/($pc+$mc)<=(1-$cut_off)){
	$type1=2;
	$type3=1;
}else{
	die "error:$chr1 can not judge phase!!\n";
}

print "judge $chr1 phase finish\n";


$mc=$pc=0;
open BP2,"<$sv_out/$chr2.barcode.phase" or die $!;
while(<BP2>){
	chomp;
	my @l=split;
	next unless exists $shareb{$l[0]};
	if($l[-1] eq "PAT"){
		$pc+=1;
	}elsif($l[-1] eq "MAT"){
		$mc+=1;
	}
}
close BP2;
if($pc/($pc+$mc)>=$cut_off){
	$type2=1;
	$type4=2;
}elsif($pc/($pc+$mc)<=(1-$cut_off)){
	$type2=2;
	$type4=1;
}else{
	die "error:$chr2 can not judge phase!!\n";
}

print "judge $chr2 phase finish\n";
$s1=$e1=$s2=$e2=0;
open REG1,"<$sv_out/$chr1.region" or die $!;
while(<REG1>){
	chomp;
	my @l=split;
	if($bp1>=$l[0] && $bp1<=$l[1]){
		$s1=$l[0];
		$e1=$l[1];
		last;
	}
}
if($s1!=0 && $e1!=0){
	print "search $chr1 block finish: start:$s1\tend:$e1\n";
}else{
	die "error:can not search $chr1 block !!!\n";
}

open REG2,"<$sv_out/$chr2.region" or die $!;
while(<REG2>){
	chomp;
	my @l=split;
	if($bp2>=$l[0] && $bp2<=$l[1]){
		$s2=$l[0];
		$e2=$l[1];
		last;
	}
}
if($s2!=0 && $e2!=0){
	print "search $chr2 block finish: start:$s2\tend:$e2\n";
}else{
	die "error:can not search $chr2 block !!!\n";
}

my($sl1,$sl2,$sl10,$sl20);
my $line1=`grep -w -n $s1 $phase_out/hapblock_L0_$chr1`;
$sl1=(split/:/,$line1)[0];
my $line2=`grep -n BLOCK $phase_out/hapblock_L0_$chr1|awk -F ":" '{print \$1}'`;
my @n1=split/\n/,$line2;
for my $k(0..scalar(@n1)-2){
	if($sl1>$n1[$k] && $sl1<$n1[$k+1]){
		$sl10=$n1[$k]+1;
		last;
	}
}
print "search $chr1 block start line finish:$sl1/$sl10\n";

my $line3=`grep -w -n $s2 $phase_out/hapblock_L0_$chr2`;
$sl2=(split/:/,$line3)[0];
my $line4=`grep -n BLOCK $phase_out/hapblock_L0_$chr2|awk -F ":" '{print \$1}'`;
my @n2=split/\n/,$line4;
for my $k(0..scalar(@n2)-2){
	if($sl2>$n2[$k] && $sl2<$n2[$k+1]){
		$sl20=$n2[$k]+1;
		last;
	}
}
print "search $chr2 block start line finish:$sl2/$sl20\n";

open OUT1,">Haplotype_withsv_1.tmp" or die $!;
open OUT2,">Haplotype_withsv_2.tmp" or die $!;
open OUT3,">Haplotype_withoutsv_1.tmp" or die $!;
open OUT4,">Haplotype_withoutsv_2.tmp" or die $!;
open BK1,"<$phase_out/hapblock_L0_$chr1" or die $!;
open BK2,"<$phase_out/hapblock_L0_$chr2" or die $!;
my ($index1,$index2);
$index1=$index2=0;
while(<BK1>){
	$index1+=1;
	next unless $index1>=$sl10;
	chomp;
	my @l=split;
	next if ($l[1] eq "-" && $l[2] eq "-");
	last if ($l[4]>$bp1);
	print OUT1 "$l[3]\t$l[4]\t$l[$type1]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
	print OUT3 "$l[3]\t$l[4]\t$l[$type3]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
}

while(<BK2>){
	$index2+=1;
	next unless $index2>=$sl20;
	last if /^\*\*/;
	chomp;
	my @l=split;
	next if ($l[1] eq "-" && $l[2] eq "-");
	if($l[4]<=$bp2){
		print OUT2 "$l[3]\t$l[4]\t$l[$type2]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
	}else{
		print OUT1 "$l[3]\t$l[4]\t$l[$type2]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
	}
	print OUT4 "$l[3]\t$l[4]\t$l[$type4]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
}
close BK2;

while(<BK1>){
	$index1+=1;
	last if /^\*\*/;
	chomp;
	my @l=split;
	next if ($l[1] eq "-" && $l[2] eq "-");
	next unless ($l[4]>$bp1);
	print OUT2 "$l[3]\t$l[4]\t$l[$type1]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
	print OUT3 "$l[3]\t$l[4]\t$l[$type1]\t$l[5]\t$l[6]\t$l[-1]\t$l[7]\n";
}
close BK1;

print "extarct het mutations finish\n";

if($vcf=~/\.gz$/){
	open VCF,"gzip -dc $vcf|" or die $!;
}else{
	open VCF,"<$vcf" or die $!;
}
while(<VCF>){
	next if /^#/;
	next unless $_=~/AF=1;/;
	chomp;
	my @l=split;
	next unless($l[0] eq $chr1|| $l[0] eq $chr2);
	if($l[0] eq $chr1){
		next unless ($l[1]>=$s1 && $l[1] <=$e1);
		print OUT3 "$l[0]\t$l[1]\t1\t$l[3]\t$l[4]\tHOM\t$l[-1]\n";
		if($l[1]>$bp1){
			print OUT2 "$l[0]\t$l[1]\t1\t$l[3]\t$l[4]\tHOM\t$l[-1]\n";
		}else{
			print OUT1 "$l[0]\t$l[1]\t1\t$l[3]\t$l[4]\tHOM\t$l[-1]\n";
		}
	}
	
	if($l[0] eq $chr2){
		next unless ($l[1]>=$s2 && $l[1] <=$e2);
		print OUT4 "$l[0]\t$l[1]\t1\t$l[3]\t$l[4]\tHOM\t$l[-1]\n";
		if($l[1]>$bp2){
			print OUT1 "$l[0]\t$l[1]\t1\t$l[3]\t$l[4]\tHOM\t$l[-1]\n";
		}else{
			print OUT2 "$l[0]\t$l[1]\t1\t$l[3]\t$l[4]\tHOM\t$l[-1]\n";
		}
	}
}
close VCF;
print "extract hom mutations finish\n";
close OUT1;
close OUT2;
close OUT3;
close OUT4;


open F1,">Haplotype_withsv_1.final.tmp" or die $!;
print F1 "chr\tposition\tmutation\tbase\treference_allele\tvariant_allele\tmismatch_quality\tvarInfo\n";
`sort -nk1 -nk2 Haplotype_withsv_1.tmp >>Haplotype_withsv_1.final.tmp`;
close F1;


open F2,">Haplotype_withsv_2.final.tmp" or die $!;
print F2 "chr\tposition\tmutation\tbase\treference_allele\tvariant_allele\tmismatch_quality\tvarInfo\n";
`sort -k 1nr -k 2n Haplotype_withsv_2.tmp >>Haplotype_withsv_2.final.tmp`;
close F2;


open F3,">Haplotype_withoutsv_1.final.tmp" or die $!;
print F3 "chr\tposition\tmutation\tbase\treference_allele\tvariant_allele\tmismatch_quality\tvarInfo\n";
`sort -nk1 -nk2 Haplotype_withoutsv_1.tmp >>Haplotype_withoutsv_1.final.tmp`;
close F3;

open F4,">Haplotype_withoutsv_2.final.tmp" or die $!;
print F4 "chr\tposition\tmutation\tbase\treference_allele\tvariant_allele\tmismatch_quality\tvarInfo\n";
`sort -nk1 -nk2 Haplotype_withoutsv_2.tmp >>Haplotype_withoutsv_2.final.tmp`;
close F4;
for my $k(1..2){
	open SV,"<Haplotype_withsv_$k.final.tmp" or die  $!;
	open OUT,">Haplotype_withsv_$k.final.list" or die $!;
	my $head=<SV>;
	print OUT $head;
	while(<SV>){
		chomp;
		my @l=split;
		if($l[2]==0){
			print OUT "$l[0]\t$l[1]\t$l[2]\t$l[3]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\n";
		}else{
			print OUT "$l[0]\t$l[1]\t$l[2]\t$l[4]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\n";
		}
	}
	close OUT;
	close SV;
}

for my $k(1..2){
	open NSV,"<Haplotype_withoutsv_$k.final.tmp" or die  $!;
	open OUT,">Haplotype_withoutsv_$k.final.list" or die $!;
	my $head=<NSV>;
	print OUT $head;
	while(<NSV>){
		chomp;
		my @l=split;
		if($l[2]==0){
			print OUT "$l[0]\t$l[1]\t$l[2]\t$l[3]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\n";
		}else{
			print OUT "$l[0]\t$l[1]\t$l[2]\t$l[4]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\n";
		}
	}
	close OUT;
	close NSV;
}
`rm *.tmp`;
print "all finish\n";


