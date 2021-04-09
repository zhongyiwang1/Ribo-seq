#!/usr/bin/perl
#use warnings;
use strict;

#######################################################################################################################
# v3.0 by Zhongyi Wang (zhongyi.wang.ds@gmail.com)
# Jan 13, 2015

#  useage: perl ### Processing3_Obtain_Read_Distributions_for_5UTR_CDS_3UTR.pl ### Transcript ID ### distribution_of_ribosome_footprints OR distribution_of_mRNA_reads ### output path
#######################################################################################################################

my ( $transcript_file, $distribution_dir, $output_path, $fl, $name, $flname,
    $oflname, $ensembl, $sum1, $sum2, $a, $b, $c, $i, $j, $k, $ip, $jp, $kp );
my ( @trans, @array1, @array2, @seq, @seqs, %hash );

print "#########################################\n\n";
$transcript_file = $ARGV[0];
if ( $transcript_file eq "" ) {
    print "No Transcript File specified!!!\n\n";
}
$distribution_dir = $ARGV[1];
if ( $distribution_dir eq "" ) {
    print "No Profile Directory specified!!!\n\n";
}
$output_path = $ARGV[2];
if ( $output_path eq "" ) {
    print "No Output Path specified!!!\n\n";
}
print "#########################################\n\n";

open IN, "$transcript_file";
while (<IN>) {
    chomp($_);
    @trans = split /\s+/, $_;
    $hash{ $trans[1] } = 1;
}
close IN;

if ($distribution_dir =~ /ribosome\_footprints/) {
    open OUT, ">>$output_path/distribution_of_ribosome_footprints_for_5UTR_CDS_3UTR.txt";
}
else {
    open OUT, ">>$output_path/distribution_of_mRNA_reads_for_5UTR_CDS_3UTR.txt";
}

print OUT "Ensembl_ID\t" . "5'UTR_reads\t" . "CDS_reads\t" . "3'UTR_reads\t" . "Total_reads\n";

$sum1 = $a = $b = $c = 0;
opendir DH, $distribution_dir;
readdir DH;
while ( $fl = readdir DH ) {
    chomp($fl);
    foreach $name ( keys %hash ) {
        chomp($name);
        if ( $fl =~ /$name/ ) {
            $flname = $distribution_dir . "/" . $fl;
            print $flname. "\n";
            if ( -f $flname ) {
                open FH, $flname;
                @seq = <FH>;
                push @seqs, @seq;
            }
            if ( -f $flname ) {
                open FH, $flname;
                $flname =~ /(distribution_of_ribosome_footprints|distribution_of_mRNA_reads)\/(.+?).txt/;
                $ensembl = $2;
                while (<FH>) {
                    chomp;
                    @array1 = split /\s+/, $_;
                    if ( $array1[7] eq "5'UTR" ) {
                        $a += $array1[5];
                    }
                    if (   $array1[7] eq "start_codon"
                        or $array1[7] eq "CDS"
                        or $array1[7] eq "stop_codon" )
                    {
                        $b += $array1[5];
                    }
                    if ( $array1[7] eq "3'UTR" ) {
                        $c += $array1[5];
                    }
                }

                $sum1 = $a + $b + $c;
                print OUT "$ensembl\t" . "$a\t" . "$b\t" . "$c\t" . "$sum1"
                  . "\n";
                $ensembl = $sum1 = $a = $b = $c = 0;
            }
        }
    }
}
closedir DH;
close FH;

@seqs = uniq(@seqs);

$sum2 = $i = $j = $k = 0;
foreach (@seqs) {
    chomp;
    @array2 = split /\t/, $_;
    if (   $array2[7] eq "5'UTR" ) { $i += $array2[5]; }
    if (   $array2[7] eq "start_codon"
        or $array2[7] eq "CDS"
        or $array2[7] eq "stop_codon" )
    {
        $j += $array2[5];
    }
    if (   $array2[7] eq "3'UTR" ) { $k += $array2[5]; }
}
$sum2 = $i + $j + $k;
$ip   = $i * 100 / $sum2;
$jp   = $j * 100 / $sum2;
$kp   = $k * 100 / $sum2;
$ip   = sprintf( "%.2f", $ip );
$jp   = sprintf( "%.2f", $jp );
$kp   = sprintf( "%.2f", $kp );

open SUM, ">$output_path/how_many_reads_are_mapped_to_transcripts_interested.txt";

print SUM "TotalTrans Reads: $sum2\n5UTR Reads: $i ($ip\%)\nCDS Reads: $j ($jp\%)\n3UTR Reads: $k ($kp\%)";

close SUM;

print "done! (^_^)\n";

###############################
#   SUBROUTINES AND MODULES   #
###############################
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}
