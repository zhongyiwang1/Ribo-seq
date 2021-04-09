#!/usr/bin/perl
#use warnings;
use strict;

#######################################################################################################################
# v2.0 by Zhongyi Wang (zhongyi.wang.ds@gmail.com)
# 8 August, 2016

#  useage: perl ### Processing2_Obtain_Read_Distributions_for_Each_Gene.pl ### Transcript_ID_file ### GTF_exon_file ### GTF_refined_file ### SAM_file ### Species ### ### Data_type ### Output_path
#######################################################################################################################
my $start  = `date +%s`;
my $start1 = `date`;
chomp( $start, $start1 );
#############################
print "\n#######################################################################################\n";
print "Script started at: ($start1)\n";
print "#######################################################################################\n\n";

my (
@sam,         @gtf,    @array,      @array1,
@array2,      @array3, @gtf_array1, @gtf_array2,
@array_reads, @arr1,   @TRANS
);
my (
%ensembl_1, %ensembl_2, %ensembl_3,  %ensembl_4, %ensembl_5,
%ensembl_6, %ensembl_7, %ensembl_8,  %ensembl_9, %ensembl_10,
%ensembl_11,%ensembl_12,%ensembl_13, %ensembl_14,%ensembl_15,
%hash1,     %hash2,     %hash3,      %hash4,     %hash_trans,
%hash_num,  %hash_exon, %hash_start, %hash_stop
);
my (
$transcript_file,  $gtf_exon_file,   $gtf_refined_file,
$data_type,        $output_path,     $sam_file,
$dir_name,         $strand_sam,      $num_extended,
$chr_sam,          $ex,              $num,
$transcript,       $number,          $trans_length,
$total,            $pos_num,         $trans_id,
$sum_reads,        $i,               $j,
$m,                $n,               $item,
$pos,              $frame,           $length_plus,
$total_trans,      $each_trans,      $trans_hash,
$chunk_num,        $x,
$species,          $offset_26,
$offset_27,        $offset_28,       $offset_29,
$offset_30,        $offset_31,       $offset_32,
$offset_33,        $offset_34,       $start_nucle_5end,
$start_nucle_3end, $stop_nucle_5end, $stop_nucle_3end
);

######### Check variables #########

$transcript_file = $ARGV[0];
if ( $transcript_file eq "" ) {
    print "No Transcript File specified!!!\n\n";
}
$gtf_exon_file = $ARGV[1];
if ( $gtf_exon_file eq "" ) {
    print "No Exon GTF specified!!!\n\n";
}
$gtf_refined_file = $ARGV[2];
if ( $gtf_refined_file eq "" ) {
    print "No Refined GTF specified!!!\n\n";
}
$sam_file = $ARGV[3];
if ( $sam_file eq "" ) {
    print "No SAM specified!!!\n\n";
}
$species = $ARGV[4];
if ( $species eq "" ) {
    print "No Species specified!!!\n\n";
}
$data_type = $ARGV[5];
if ( $data_type eq "" ) {
    print "No Data type specified!!!\n\n";
}
$output_path = $ARGV[6];
if ( $output_path eq "" ) {
    print "No Output Path specified!!!\n\n";
}

####################################
#        SCRIPT   START            #
####################################
#First create a directory to store all of the files to generate.

if ( $data_type eq "RP" ) {
    $dir_name = "distribution_of_ribosome_footprints";
}
elsif ( $data_type eq "TR" ) {
    $dir_name = "distribution_of_mRNA_reads";
}
else { print "There is no data type specified!"; }

system("mkdir $output_path/$dir_name");
print "\nAll outputs will be stored in directory \"$dir_name\"\n";

$total_trans = `wc -l $transcript_file | awk {'print \$1 + 1'}`;
print "The total number of transcripts is $total_trans\n";

if ( $total_trans > 2000 ) {
    $each_trans = $total_trans / 15;
    $each_trans = sprintf( "%.0f", $each_trans );
    print "The number of transcripts for each chunk is $each_trans\n";
    $number = 1;
    open TRANSCRIPT, "$transcript_file";
    while (<TRANSCRIPT>) {
        chomp;
        @array = split /\s+/;
        if ( $number > 0 * $each_trans and $number <= 1 * $each_trans ) {
            $ensembl_1{ $array[1] } = 1;
        }
        elsif ( $number > 1 * $each_trans and $number <= 2 * $each_trans ) {
            $ensembl_2{ $array[1] } = 1;
        }
        elsif ( $number > 2 * $each_trans and $number <= 3 * $each_trans ) {
            $ensembl_3{ $array[1] } = 1;
        }
        elsif ( $number > 3 * $each_trans and $number <= 4 * $each_trans ) {
            $ensembl_4{ $array[1] } = 1;
        }
        elsif ( $number > 4 * $each_trans and $number <= 5 * $each_trans ) {
            $ensembl_5{ $array[1] } = 1;
        }
        elsif ( $number > 5 * $each_trans and $number <= 6 * $each_trans ) {
            $ensembl_6{ $array[1] } = 1;
        }
        elsif ( $number > 6 * $each_trans and $number <= 7 * $each_trans ) {
            $ensembl_7{ $array[1] } = 1;
        }
        elsif ( $number > 7 * $each_trans and $number <= 8 * $each_trans ) {
            $ensembl_8{ $array[1] } = 1;
        }
        elsif ( $number > 8 * $each_trans and $number <= 9 * $each_trans ) {
            $ensembl_9{ $array[1] } = 1;
        }
        elsif ( $number > 9 * $each_trans and $number <= 10 * $each_trans ) {
            $ensembl_10{ $array[1] } = 1;
        }
        elsif ( $number > 10 * $each_trans and $number <= 11 * $each_trans ) {
            $ensembl_11{ $array[1] } = 1;
        }
        elsif ( $number > 11 * $each_trans and $number <= 12 * $each_trans ) {
            $ensembl_12{ $array[1] } = 1;
        }
        elsif ( $number > 12 * $each_trans and $number <= 13 * $each_trans ) {
            $ensembl_13{ $array[1] } = 1;
        }
        elsif ( $number > 13 * $each_trans and $number <= 14 * $each_trans ) {
            $ensembl_14{ $array[1] } = 1;
        }
        elsif ( $number > 14 * $each_trans and $number <= 16 * $each_trans ) {
            $ensembl_15{ $array[1] } = 1;
        }
        else { next; }
        $number++;
    }
    @TRANS = (  \%ensembl_1,  \%ensembl_2,  \%ensembl_3,  \%ensembl_4,  \%ensembl_5,
    \%ensembl_6,  \%ensembl_7,  \%ensembl_8,  \%ensembl_9,  \%ensembl_10,
    \%ensembl_11, \%ensembl_12, \%ensembl_13, \%ensembl_14, \%ensembl_15  );
    close TRANSCRIPT;
}
else {
    open TRANSCRIPT, "$transcript_file";
    while (<TRANSCRIPT>) {
        chomp;
        @array = split /\s+/;
        $ensembl_1{ $array[1] } = 1;    #optional: 0 or 1
    }
    @TRANS = ( \%ensembl_1 );
    close TRANSCRIPT;
}

########## STEP 1 ##########
$chunk_num = 1;
foreach $trans_hash (@TRANS) {
    print "\n\n##############\n";
    print "Chunk $chunk_num\n";
    print "##############\n\n";
    
    my $step1starttime = `date`;
    chomp($step1starttime);
    print "\n#######################################################################################\n";
    print "STEP 1 (Process GTF File) started at: ($step1starttime)\n";
    print "#######################################################################################\n\n";
    
    $number = 1;
    foreach $transcript ( sort keys %$trans_hash ) {
        print "No.$number\t$transcript\t";
        
        open GTF, "$gtf_exon_file" or die($!);
        while (<GTF>) {
            chomp;
            @gtf = split /\s+/, $_;
            chomp($transcript);
            
            if ( $gtf[6] eq $transcript ) {
                if ( $gtf[1] eq "exon" ) {
                    for ( $ex = $gtf[2] ; $ex <= $gtf[3] ; $ex++ ) {
                        $hash1{$ex} = 1;
                    }
                }
                $strand_sam = $gtf[4];
                $chr_sam    = $gtf[0];
                $chr_sam =~ s/^chr//i; # remove "chr" from the chromosome name if exists.
            }
        }
        print "Strand\:$strand_sam\tChromosome\: $chr_sam\n";
        $trans_length = keys %hash1;
        $hash_trans{$transcript} = $trans_length;
        
        $num = 1;
        if ( $strand_sam eq "+" ) {
            foreach $pos_num ( sort { $a <=> $b } keys %hash1 ) {
                push @gtf_array1, "$transcript\t$strand_sam\t$chr_sam\t$pos_num\t$num";
                $num++;
            }
        }
        elsif ( $strand_sam eq "-" ) {
            foreach $pos_num ( sort { $b <=> $a } keys %hash1 ) {
                push @gtf_array1, "$transcript\t$strand_sam\t$chr_sam\t$pos_num\t$num";
                $num++;
            }
        }
        
        $number++;
        %hash1      = ();
        @gtf        = ();
        $strand_sam = $chr_sam = $trans_length = $ex = $pos_num = 0;
        
    }
    close GTF;
    
    foreach (@gtf_array1) {
        chomp;
        @gtf_array2 = split /\s{1,20}/, $_;
        
        # Ensembl_id--strand--chromosome--genome_coord--transcript_coord
        $hash2{ $gtf_array2[1] }{ $gtf_array2[2] }{ $gtf_array2[3] } = $gtf_array2[4];    # values are transcript coordinates
        $hash3{ $gtf_array2[1] }{ $gtf_array2[2] }{ $gtf_array2[3] } = $gtf_array2[0];    # values are transcript Ensembl IDs
    }
    
    ########## STEP 2 ##########
    my $step2starttime = `date`;
    chomp($step2starttime);
    print "\n#######################################################################################\n";
    print "STEP 2 (Process SAM File) started at: ($step2starttime)\n";
    print "#######################################################################################\n\n";
    
    if ( $data_type eq "RP" ) {    ## Ribo-Seq data
        
        $offset_26 = 15;
        $offset_27 = 15;
        $offset_28 = 15;
        $offset_29 = 15;
        $offset_30 = 15;
        $offset_31 = 15;
        $offset_32 = 15;
        $offset_33 = 15;
        $offset_34 = 15;
        
        open SAM, "$sam_file" or die($!);
        while (<SAM>) {
            chomp;
            @sam = split /\s{1,20}/, $_;
            
            # strand--chromosome--genome_coord--read_length--reads
            if ( $sam[0] eq "+" ) {
                if ( $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } ) {
                    if ( $sam[3] == 26 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_26;
                    }
                    elsif ( $sam[3] == 27 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_27;
                    }
                    elsif ( $sam[3] == 28 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_28;
                    }
                    elsif ( $sam[3] == 29 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_29;
                    }
                    elsif ( $sam[3] == 30 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_30;
                    }
                    elsif ( $sam[3] == 31 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_31;
                    }
                    elsif ( $sam[3] == 32 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_32;
                    }
                    elsif ( $sam[3] == 33 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_33;
                    }
                    elsif ( $sam[3] == 34 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $offset_34;
                    }
                    else { next; }
                    
                    push @{ $hash4{"$hash3{$sam[0]}{$sam[1]}{$sam[2]}\t$num_extended"} }, $sam[4];
                }
            }
            elsif ( $sam[0] eq "-" ) {
                if ( $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } ) {
                    if ( $sam[3] == 26 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_26;
                    }
                    elsif ( $sam[3] == 27 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_27;
                    }
                    elsif ( $sam[3] == 28 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_28;
                    }
                    elsif ( $sam[3] == 29 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_29;
                    }
                    elsif ( $sam[3] == 30 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_30;
                    }
                    elsif ( $sam[3] == 31 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_31;
                    }
                    elsif ( $sam[3] == 32 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_32;
                    }
                    elsif ( $sam[3] == 33 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_33;
                    }
                    elsif ( $sam[3] == 34 ) {
                        $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + 1 + $offset_34;
                    }
                    else { next; }
                    
                    push @{ $hash4{"$hash3{$sam[0]}{$sam[1]}{$sam[2]}\t$num_extended"} }, $sam[4];
                }
            }
        }
    }
    elsif ( $data_type eq "TR" ) {    ## RNA-Seq data
        
        open SAM, "$sam_file" or die($!);
        
        while (<SAM>) {
            
            chomp;
            
            @sam = split /\s{1,20}/, $_;
            
            if ( $sam[0] eq "+" ) {
                
                if ( $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } ) {
                    
                    # $length_plus = int( $sam[3] / 2 + 0.5 );    # middle point of reads
                    
                    $length_plus = 15;    # assumed Asite of reads
                    
                    $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } + $length_plus;
                    
                    push @{ $hash4{"$hash3{$sam[0]}{$sam[1]}{$sam[2]}\t$num_extended"} }, $sam[4];
                }
            }
            elsif ( $sam[0] eq "-" ) {
                
                if ( $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } ) {
                    
                    # $length_plus = int( $sam[3] / 2 + 0.5 );    # middle point of reads
                    
                    $length_plus = 15;    # assumed Asite of reads
                    
                    $num_extended = $hash2{ $sam[0] }{ $sam[1] }{ $sam[2] } - $sam[3] + $length_plus + 1;
                    
                    push @{ $hash4{"$hash3{$sam[0]}{$sam[1]}{$sam[2]}\t$num_extended"} }, $sam[4];
                }
            }
        }
    }
    else { last; }
    
    ########## STEP 3 ##########
    my $step3starttime = `date`;
    chomp($step3starttime);
    print "\n#######################################################################################\n";
    print "STEP 3 (Get Footprints For Each Transcript) started at: ($step3starttime)\n";
    print "#######################################################################################\n\n";
    
    $total     = 0;
    $number    = 1;
    $sum_reads = 0;
    foreach $trans_id ( sort keys %hash_trans ) {
        chomp($trans_id);
        print "No.$number\t$trans_id\t";
        for ( $num = 1 ; $num <= $hash_trans{$trans_id} ; $num++ ) {
            foreach ( @{ $hash4{"$trans_id\t$num"} } ) {
                chomp;
                $total = $total + $_;
            }
            push @array1, "$num\t$total\n";
            $total = 0;
        }
        
        foreach (@array1) {
            chomp;
            @array_reads = split /\s+/;
            $sum_reads += $array_reads[1];
            $hash_num{ $array_reads[0] } = $array_reads[1];
        }
        
        if ( $sum_reads > 0 ) {    # keep transcripts with over 0 alignments
            
            open GTF, "$gtf_refined_file" or die($!);    ## refined GTF file
            while (<GTF>) {
                chomp;
                @gtf = split /\s+/, $_;
                
                if ( $gtf[6] eq $trans_id ) {
                    if ( $gtf[1] eq "exon" ) {
                        for ( $i = $gtf[2] ; $i <= $gtf[3] ; $i++ ) {
                            $hash_exon{$i} = 1;
                        }
                    }
                    if ( $gtf[1] eq "mainORF" or $gtf[1] eq "CDS" ) { # mainORF for Evgeny's Amniote gtf annotation; CDS for Ensembl's Drosophila gtf annotation
                        if ( $gtf[4] eq "+" ) {
                            for ( $m = $gtf[2] ; $m <= $gtf[2] + 2 ; $m++ ) {
                                $hash_start{$m}   = 1;
                                $start_nucle_5end = $gtf[2];
                                $start_nucle_3end = $gtf[2] + 2;
                            }
                            for ( $n = $gtf[3] + 1 ; $n <= $gtf[3] + 3 ; $n++ )
                            {    # stop codon is not included in the ORF
                                $hash_stop{$n}   = 1;
                                $stop_nucle_5end = $gtf[3] + 1;
                                $stop_nucle_3end = $gtf[3] + 3;
                            }
                        }
                        elsif ( $gtf[4] eq "-" ) {
                            for ( $m = $gtf[3] ; $m >= $gtf[3] - 2 ; $m-- ) {
                                $hash_start{$m}   = 1;
                                $start_nucle_5end = $gtf[3];
                                $start_nucle_3end = $gtf[3] - 2;
                            }
                            for ( $n = $gtf[2] - 1 ; $n >= $gtf[2] - 3 ; $n-- )
                            {    # stop codon is not included in the ORF
                                $hash_stop{$n}   = 1;
                                $stop_nucle_5end = $gtf[2] - 1;
                                $stop_nucle_3end = $gtf[2] - 3;
                            }
                        }
                    }
                    $strand_sam = $gtf[4];
                    $chr_sam    = $gtf[0];
                    $chr_sam =~ s/^chr//i; # remove "chr" from the chromosome name if there is a.
                }
            }
            
            $num = 1;
            if ( $strand_sam eq "+" ) {
                foreach $item ( sort { $a <=> $b } keys %hash_exon ) {
                    if (    $item > $start_nucle_3end and $item < $stop_nucle_5end   )
                    {
                        push( @arr1, "$num\t$item\tCDS\n" );
                        # trans_coord--genomic_coord--feature
                    }    #Ensembl GTF CDS includes start codon, not stop codon
                    elsif ( $hash_start{$item} ) {
                        push( @arr1, "$num\t$item\tstart_codon\n" );
                    }
                    elsif ( $hash_stop{$item} ) {
                        push( @arr1, "$num\t$item\tstop_codon\n" );
                    }
                    elsif ( $item < $start_nucle_5end ) {
                        push( @arr1, "$num\t$item\t5\'UTR\n" );
                    }
                    elsif ( $item > $stop_nucle_3end ) {
                        push( @arr1, "$num\t$item\t3\'UTR\n" );
                    }
                    $num++;
                }
                
            }
            elsif ( $strand_sam eq "-" ) {
                foreach $item ( sort { $b <=> $a } keys %hash_exon ) {
                    if (    $item < $start_nucle_3end and $item > $stop_nucle_5end  )
                    {
                        push( @arr1, "$num\t$item\tCDS\n" );
                    }    #Ensembl GTF CDS includes start codon, not stop codon
                    elsif ( $hash_start{$item} ) {
                        push( @arr1, "$num\t$item\tstart_codon\n" );
                    }
                    elsif ( $hash_stop{$item} ) {
                        push( @arr1, "$num\t$item\tstop_codon\n" );
                    }
                    elsif ( $item > $start_nucle_5end ) {
                        push( @arr1, "$num\t$item\t5\'UTR\n" );
                    }
                    elsif ( $item < $stop_nucle_3end ) {
                        push( @arr1, "$num\t$item\t3\'UTR\n" );
                    }
                    $num++;
                }
            }
            
            $frame = 1;    # add the frame information
            open OUT, ">$output_path/$dir_name/$trans_id.txt";
            
            print OUT "Transcript\t" . "Strand\t" . "Chromosome\t" . "Geno_Coord\t" . "Trans_Coord\t" . "Count\t" . "Frame\t" . "Region\n";
            
            foreach (@arr1) {
                chomp;
                @array2 = split /\s+/;
                
                # Trans_Name--Strand--Chromosome--Geno_Coord--Trans_Coord--Count--Frame--Region(Feature)
                print OUT "$trans_id\t$strand_sam\t$chr_sam\t$array2[1]\t$array2[0]\t$hash_num{$array2[0]}\t$frame\t$array2[2]\n";
                
                $frame++;
                if ( $frame > 3 ) {
                    $frame = 1;
                }
            }
            
            print "Total reads\: $sum_reads\n";
            $number++;
        }
        else {
            print "Total reads\: $sum_reads\n";
            $number++;
        }
        
        @gtf = @gtf_array1 = @array1 = @array2 = @array3 = @arr1 = ();
        %hash_num = %hash_exon = %hash_start = %hash_stop = ();
        $strand_sam = $chr_sam = $total = $pos_num = $trans_id = $sum_reads =
        $i = $j = $m = $n = $item = $pos = 0;
    }
    %hash2 = %hash3 = %hash4 = %hash_trans = ();
    $chunk_num++;
}

my $finish  = `date +%s`;
my $finish1 = `date`;
chomp( $finish, $finish1 );
#############################
print "\n#######################################################################################\n";
print "Script finished at: ($finish1)\n";
print "#######################################################################################\n\n";

print "\n#######################################################################################";
runtime();
print "#######################################################################################\n\n";

###############################
#   SUBROUTINES AND MODULES   #
###############################
#This subroutine prints out the total runtime of the program. Just add runtime() at the end of the script.
sub runtime {    #Make sure that "$start = `date +%s`;" comes after sub def.
    my $hours;
    $hours = $finish - $start;
    printf(
    "\nTotal running time: %02d:%02d:%02d\n",
    int( $hours / 3600 ),
    int( ( $hours % 3600 ) / 60 ),
    int( $hours % 60 )
    );
}
