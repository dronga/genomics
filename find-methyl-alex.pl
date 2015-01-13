#!/usr/bin/perl -w 

#use strict;

my $fastafile = shift or die "FASTA file not specified\n";



open my $fastafh, "$fastafile" or die "could not open fasta file $fastafile\n";

my $markercount = 0;
while (my $line = <$fastafh>) {
#       chomp($line);             .
        $markercount++;
        # do line-by-line processing.
        if(! substr($line,0,1) eq '>') {
        die "defline didn't start with > (defline was $line)\n";
        }
        my %dl1 ;
        ($dl1{chr},$dl1{bp},$dl1{rsid},$dl1{allele},my @lrf)= split(/_/, $line);
        my $seq1 = <$fastafh>;
                
        chomp($seq1);
        $line = <$fastafh>;
#       chomp($line);
        my %dl2;
        ($dl2{chr},$dl2{bp},$dl2{rsid},$dl2{allele},@lrf)= split(/_/, $line);
#       if($dl1{bp} ne $dl2{bp} or $dl1{rsid} ne $dl2{rsid}) { die "second record doesn't match first: $dl1"}
        my $seq2 = <$fastafh>;
        chomp($seq2);
        my $cpgstatus = 0;
        if($seq1 =~ m/CG/) {$cpgstatus=1}
        if($seq2 =~ m/CG/) {$cpgstatus=2}
        if($seq1 =~ m/CG/ and $seq2 =~ m/CG/)  {$cpgstatus=3}



        #if($cpgstatus>0){print "$markercount $dl1{bp} $dl1{rsid} $seq1 $seq2 $cpgstatus\n";}
#       }
        if($cpgstatus>0){print "$markercount $seq1 $seq2 $cpgstatus\n";}
}
close($fastafh);
