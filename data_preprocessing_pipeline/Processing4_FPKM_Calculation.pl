#!/usr/bin/perl
#use warnings;
use strict;

#######################################################################################################################
# v3.0 by Zhongyi Wang (zhongyi.wang.ds@gmail.com)
# Jan 13, 2015

#  useage: perl ### Processing4_FPKM_Calculation.pl ### Lengths of 5UTR, CDS and 3UTR  ### Footprint distributions for 5'UTR, CDS and 3'UTR ### millions of alignments mapped to transcriptome ### output path
#######################################################################################################################

my ( $CDS_UTR_length, $read_distribution, $type, $millions_of_alignments, $output_path, $i, $j, $output_name, $k, $FPKM, $RPM );
my ( @array1, @array2, @array3, @arr, @uniq );
my ( %hash ) ;

print "#########################################\n\n";
$CDS_UTR_length = $ARGV[0];
if ( $CDS_UTR_length eq "" ) {
    print "No Transcript Length File specified!!!\n\n";
}
$read_distribution = $ARGV[1];
if ( $read_distribution eq "" ) {
    print "No Transcript File Specified!!!\n\n";
}
$type = $ARGV[2];
if ( $type eq "" ) {
    print "No Library Type Specified!!!\n\n";
}
$millions_of_alignments = $ARGV[3];
if ( $millions_of_alignments eq "" ) {
    print "No Total Reads specified!!!\n\n";
}
$output_path = $ARGV[4];
if ( $output_path eq "" ) {
    print "No Output Path specified!!!\n\n";
}

if ( $read_distribution =~ /footprints/ ) {
    $output_name = "ribosome_footprints";
}
elsif ( $read_distribution =~ /reads/ ) {
    $output_name = "mRNA_reads";
}
print "#########################################\n\n";

open IN1, "$read_distribution";
while (<IN1>) {
    chomp;
    @array1 = split /\s+/;
    # Ensembl_ID : 5'UTR_reads : CDS_reads : 3'UTR_reads : Total_reads
    $hash{ $array1[0] } = "$array1[1]\t$array1[2]\t$array1[3]\t$array1[4]";
}
close IN1;

open IN2, "$CDS_UTR_length";
while (<IN2>) {
    chomp;
    @array2 = split /\s+/;
    if ($_ !~ /Ensembl_ID/) {
        if ( $hash{ $array2[0] } ) {
            # Ensembl_ID : 5'UTR_reads : CDS_reads : 3'UTR_reads : Total_reads : 5'UTR_length : CDS_length : 3'UTR_length : Total_length
            push( @arr, "$array2[0]\t$hash{$array2[0]}\t$array2[1]\t$array2[2]\t$array2[3]\t$array2[4]" );
        }
        else {
            push( @arr, "$array2[0]\t0\t0\t0\t0\t$array2[1]\t$array2[2]\t$array2[3]\t$array2[4]" );
        }
    }
}
@uniq = uniq(@arr);
close IN2;

open OUT, ">$output_path/FPKM_and_density.txt";
print OUT "Ensembl_ID\t5UTR_reads\tCDS_reads\t3UTR_reads\tTotal_reads\t5UTR_length\tCDS_length\t3UTR_length\tTotal_length\t5UTR_ratio\tCDS_ratio\t3UTR_ratio\tRPM\tFPKM\n";

foreach ( @uniq ) {
    chomp;
    @array3 = split /\s+/;
    # Ensembl_ID : 5'UTR_reads : CDS_reads : 3'UTR_reads : Total_reads : 5'UTR_length : CDS_length : 3'UTR_length : Total_length
    
    if ( $array3[5] > 0 ) {
        $i = $array3[1] / $array3[5]; # 5'UTR
    } else { $i = "NA"; }
    
    if ( $array3[6] > 0 ) {
        $j = $array3[2] / $array3[6]; # CDS
    } else { $j = "NA"; }
    
    if ( $array3[7] > 0 ) {
        $k = $array3[3] / $array3[7]; # 3'UTR
    } else { $k = "NA"; }
    
    # if ($type eq "RP") { ## In case if I want to use all transcript reads, instead of cds reads, for FPKM calculation for TR.
    
    $RPM  = $array3[2] / $millions_of_alignments;  # alignments per million of CDS reads (RPM)
    
    if ($j ne "NA") {
        $FPKM = ( $j * 1000 ) / $millions_of_alignments; # alignments per kilobase per million of CDS reads (FPKM)
    }
    else { $FPKM = "NA"; }
    
    if ( $i ne "NA" ) { $i = sprintf( "%.4f", $i ); }
    else              { $i = "NA"; }
    if ( $j ne "NA" ) { $j = sprintf( "%.4f", $j ); }
    else              { $j = "NA"; }
    if ( $k ne "NA" ) { $k = sprintf( "%.4f", $k ); }
    else              { $k = "NA"; }
    if ( $RPM ne "NA" ) { $RPM  = sprintf( "%.4f", $RPM ); }
    else                { $RPM = "NA"; }
    if ( $FPKM ne "NA" ) { $FPKM = sprintf( "%.4f", $FPKM ); }
    else                 { $FPKM = "NA"; }
    
    print OUT "$array3[0]\t$array3[1]\t$array3[2]\t$array3[3]\t$array3[4]\t$array3[5]\t$array3[6]\t$array3[7]\t$array3[8]\t$i\t$j\t$k\t$RPM\t$FPKM\n";
    
    @array3 = ();
}
close OUT;


###############################
#   SUBROUTINES AND MODULES   #
###############################
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

