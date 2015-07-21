## Convert rs IDs to chr:pos in plink .clumped files using a bim reference file                                                                                                                             
## written by Alex Drong 21/07/15                                                                                                                                                                           
## usage : perl clump_rs_to_chrpos.pl ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bim PCOS-ALL-23andMe_clumped_5e-8_noLD_1KGv3_AllVariants.clumped > PCOS-ALL-23andMe_clumped_5\
e-8_noLD_1KGv3_AllVariants.clumped.chrpos                                                                                                                                                                   

#!/usr/bin/perl                                                                                                                                                                                             

#use strict;                                                                                                                                                                                                


my $bimfile = shift or die "bim file not specified\n";
open $bimfile, "$bimfile" or die "could not open bim file $bimfile\n";

my %bim ;
while (my $line = <$bimfile>) {
 chomp($line);
 (my $chr, my $rs, my $strand, my $pos, my $a1, my $a2) = split(/\t/, $line);
 my $chrpos = $chr . ":" . $pos;
 $bim{$rs} = $chrpos;

}


my $clumpfile = shift or die "clump file not specified\n";
open $clumpfile, "$clumpfile" or die "could not open fasta file $clumpfile\n";

my $header = <$clumpfile>;
print $header;

while (my $line = <$clumpfile>) {
    chomp($line);
    my @aline = split(/ +/, $line);
    if($#aline==12){
        $aline[3] = $bim{$aline[3]};
        $aline[12] =~ s/\(1\)//g;
        @snps = split(",", $aline[12]);
        my @asnps = @bim{@snps};
        $aline[12] = join(",", @asnps);
        shift @aline;
        print join("\t", @aline), "\n";
    }
}
