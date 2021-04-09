#!/usr/bin/perl
#use warnings;
use strict;

#######################################################################################################################
# v3.0 by Zhongyi Wang (zhongyi.wang.ds@gmail.com)
# Jan 19, 2015

#  useage: perl ### Processing1_Refine_Index_Sam_File.pl ### SAM file ### Lower limit of read size ### Upper limit of read size ### Output_path

# ribosome footprint lengths are between 26 and 34
#          mRNA read lengths are between 20 and 29
#######################################################################################################################
my ( $sam_file, $threshold_low, $threshold_high, $output_path, $length, $key );
my ( @sam, %hash );
print "#########################################\n\n";
$sam_file = $ARGV[0];
if ( $sam_file eq "" ) {
    print "invalid SAM file\n";
}
$threshold_low = $ARGV[1];
if ( $threshold_low eq "" ) {
    print "invalid number\n";
}
$threshold_high = $ARGV[2];
if ( $threshold_high eq "" ) {
    print "invalid number\n";
}
$output_path = $ARGV[3];
if ( $output_path eq "" ) {
    print "invalid path\n";
}
print "#########################################\n\n";

open SAM, "$sam_file" or die($!);
while (<SAM>) {
    chomp;
    @sam = split /\s{1,20}/, $_;
    
    if ( $sam[1] eq "0" ) { # for Evo-Devo sam files, $sam[1] eq "16" indicates forward strand
        
        if ( $sam[5] =~ /^(\d+)M$/ )
        {    # for reads mapping within a single exon
            $length = $1;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"+\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
        elsif ( $sam[5] =~ /^(\d+)M(\d+)N(\d+)M$/ )
        {    # for reads spanning two exons
            $length = $1 + $3;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"+\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
        elsif ( $sam[5] =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/ )
        {    # for reads spanning three exons
            $length = $1 + $3 + $5;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"+\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
        elsif ( $sam[5] =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/ )
        {    # for reads spanning four exons
            $length = $1 + $3 + $5 + $7;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"+\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
    }
    elsif ( $sam[1] eq "16" ) { # for Evo-Devo sam files, $sam[1] eq "0" indicates reverse strand
        
        if ( $sam[5] =~ /^(\d+)M$/ )
        {    # for reads mapping within a single exon
            $length = $1;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"-\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
        elsif ( $sam[5] =~ /^(\d+)M(\d+)N(\d+)M$/ )
        {    # for reads spanning two exons
            $length = $1 + $3;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"-\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
        elsif ( $sam[5] =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/ )
        {    # for reads spanning three exons
            $length = $1 + $3 + $5;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"-\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
        elsif ( $sam[5] =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/ )
        {    # for reads spanning four exons
            $length = $1 + $3 + $5 + $7;
            if ( $length >= $threshold_low and $length <= $threshold_high ) {
                $hash{"-\t$sam[2]\t$sam[3]\t$length"}++;
                $length = 0;
            }
        }
    }
}

open OUT, ">$output_path/unique_mappers_frame_analysis_index.sam" or die($!);
print OUT "strand\tchromosome\tgenome_coord\tread_length\treads\n";
foreach $key ( sort keys %hash ) {
    print OUT "$key\t$hash{$key}\n";
}
close OUT;

