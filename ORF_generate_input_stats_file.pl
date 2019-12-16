#! /usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use Array::Utils qw(:all);
use Data::Dumper;

my $inputfile1 = shift @ARGV;

unless (defined $inputfile1){
    print "\nPlease provide the path of the ORF header file: \n";
    print "\n\nUsage: \n";
    print "          perl $0 ~/ZZG/final.contigs.fa.fna.out.header \n\n";
        
    exit;
}

my $path = abs_path($inputfile1);
my ($file_name, $file_dir) = fileparse($path);

my $outputfile1 = $file_dir. '/'.  'ORF_lengths.txt';
my $outputfile2 = $file_dir. '/'.  'ORF_partials.txt';
my $outputfile3 = $file_dir. '/'.  'ORF_start_types.txt';
my $outputfile4 = $file_dir. '/'.  'ORF_gc_conts.txt';
my (@lengths, @partials, @start_types) = ();

#open (OFH1, OFH2, OFH3, OFH4), <($outputfile1, $outputfile2, $outputfile3, $outputfile4) or die $!;
open (OFH1, ">$outputfile1") or die $!;
open (OFH2, ">$outputfile2") or die $!;
open (OFH3, ">$outputfile3") or die $!;
open (OFH4, ">$outputfile4") or die $!;

print "Generating files for:\n", $outputfile1, "\n", $outputfile2, "\n", $outputfile3, "\n", $outputfile4, "\n\n";

open (FH1, "<$inputfile1") or die $!;
while(<FH1>){
    chomp;
    my @array = split(/#/);
    my $start = $array[1]; my $end = $array[2];
    my $gene_length = abs($end - $start) + 1;
    push @lengths, $gene_length;
    #print OFH1 $gene_length, "\t", $array[0], "\n";
    print OFH1 $gene_length, "\n";
    my @array2 = split(/=|;/, $array[4]);
    my $partial_mark = $array2[3];
    my $start_type = $array2[5];
    my $gc_cont = $array2[11];
    print OFH2 $partial_mark, "\n";
    print OFH3 $start_type, "\n";
    print OFH4 $gc_cont, "\n";
    #push @partials, $partial_mark;
    #push @start_types, $start_type;
}
close(FH1);
close(OFH1);close(OFH2);close(OFH3);close(OFH4);


__DATA__
>k99_48_1 # 3 # 215 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.545
>k99_48_2 # 233 # 1105 # 1 # ID=1_2;partial=01;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10bp;gc_cont=0.505
>k99_63_1 # 1 # 885 # 1 # ID=2_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.304
>k99_63_2 # 1048 # 1620 # 1 # ID=2_2;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.216
>k99_69_1 # 1 # 1635 # 1 # ID=3_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.601
>k99_74_1 # 2 # 253 # 1 # ID=4_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.563
>k99_74_2 # 250 # 723 # 1 # ID=4_2;partial=00;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10bp;gc_cont=0.616
>k99_74_3 # 750 # 1172 # 1 # ID=4_3;partial=00;start_type=TTG;rbs_motif=GGA/GAG/AGG;rbs_spacer=11-12bp;gc_cont=0.532
>k99_74_4 # 1307 # 1561 # 1 # ID=4_4;partial=01;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.557
>k99_88_1 # 3 # 1211 # 1 # ID=5_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.602