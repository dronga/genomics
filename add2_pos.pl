# This adds the second position to the results file to make it a bed file 

$snpfile = shift or die "FASTA file not specified\n";
open $snpfile, "$snpfile" or die "could not open fasta file $snpfile\n";
while (my $line = <$snpfile>) {
  chomp($line);
  
  (my @tarray) = split(/[ \t]/, $line);
  splice @tarray, 2, 0, ($tarray[1] + 1);
#  print($tarray[0] . "\t"  . $tarray[1] . "\t" .  ($tarray[1] +1) . "\n"); 
#  print($tarray[0] . "\t"  . $tarray[1] . "\t " . ($tarray[1] + 1) . "\t" . join("\t", @tarray[2 .. $#array]) . "\n");
  print join("\t", @tarray) . "\n";
}
