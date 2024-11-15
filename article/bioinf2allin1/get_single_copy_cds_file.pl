#!/usr/bin/perl -w
use FindBin qw($Bin $Script);

BEGIN{
     unshift (@INC,"$Bin/Soft/");
}

my %hash;
my %seq;
use Bio::SeqIO;
my $nuc=$ARGV[0];
my $single=$ARGV[1];
$ina = Bio::SeqIO->new(-file => $nuc, -format => 'fasta');
while(my $obj = $ina->next_seq()){
      my $id = $obj->id;
	
      $seq{$id}=$obj->seq;
}
my @file=glob("$single/OG*.fa");
for(my $i=0;$i<=$#file;$i++){
   if($file[$i]=~/OG\d+\.fa$/){
   open(OUT,">$file[$i].mafft.fa.cds");
    my $ina = Bio::SeqIO->new(-file => $file[$i], -format => 'fasta');
    while(my $obj = $ina->next_seq()){
        my $id = $obj->id;
         print OUT ">$id\n$seq{$id}\n";
    }
    close OUT;
   }



}
