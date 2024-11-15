#!/usr/bin/perl -w
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
BEGIN{

unshift (@INC,"$Bin/");
}
use strict;
use Bio::SeqIO;
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fgff,$fgenome,$pcfg, $index,$filter_flag, $outdir);
$pcfg="$Bin/../parameter.cfg";
GetOptions(
				"help|?" =>\&USAGE,
				"gff:s"=>\$fgff,
				"genome:s"=>\$fgenome,
				"filter:s"=>\$filter_flag,
				"od:s"=>\$outdir,
				"index:s"=>\$index,
				"pcfg:s"=>\$pcfg,
				) or &USAGE;
&USAGE unless ($fgff and $fgenome );
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
`mkdir $outdir/Raw` unless (-d "$outdir/Raw");
my %gff;


# ------------------------------------------------------------------
# gff convert
# ------------------------------------------------------------------
#convert:
#1: only gene mRNA CDS
#2: format convert to "ID=xxx;Parent=xxx"
#3: ID: cut "-TA/-PA", cut "mRNA:/gene:"
#4: ID: "|" convert to "_"
my %hash_gene_CDS;
open (GFF, ">$outdir/Raw/$index.gff3.tmp") or die $!;
open (IN, $fgff) or die $!;
my $gene_id ;
my %have;
my $mRNA_id;
my $is_gene=0;
my $gff_num = 0;
my $non_gene = 0; ## mark gene is real gene or not
my $pos;
while (<IN>){
	chomp;
	next if (m/^\#/) ;
	next if($_!~/\w+/);
	my @unit = split /\t/, $_;
	my ($chr, undef, $type, $star, $end, undef, $ori, $phase, $info) = @unit;
	next if($type !~/gene/ && $type ne "mRNA" && $type ne "CDS" &&$type ne "exon");
	# filter TE gene, filter
	if ($type eq "transposable_element_gene"||($type eq  "pseudogene")||$info=~/gene-Mir/){
		$non_gene = 1;
		next;
	}elsif($type =~/gene/ ){
		$non_gene = 0;
	}
	next if ($non_gene == 1);
	if ($type eq "gene") {
		$is_gene = 1;
	}
	### check star end
	if ($star > $end) {
		print "Wrong Gff format: star > end!\n";
		exit(1);
	}
	if ($type eq "gene") {
		$pos="";
		if($info=~/Name=([^;]+)/){
			$gene_id=$1;
 			$gene_id=~s/\(/\_/g;
                       $gene_id=~s/\)/\_/g;

			my($idgene)=$info=~/ID=([^;]+)/;
			if($gene_id=~/$idgene/||$idgene=~/$gene_id/){
				$gene_id=$idgene;
			}
		}elsif ($info=~/ID=/) {
			($gene_id) = $info=~/ID=([^;]+)/;
		}elsif($info=~/Gene /){
			($gene_id)= $info=~/Gene (\S+)/;
		}else{
			print "Gff Format is Something else, please check!!!\n" ;
			exit(1) ;
		}
		$gene_id=~s/gene://;
		$gene_id=~s/:/_/g;
		$gene_id=~s/\|/_/g;
		$gene_id=~s/\s+/\_/g;
		 $gene_id=~s/\%//g;
		if(defined $have{$gene_id}){
			$have{$gene_id}++;
			$gene_id="$gene_id"."_$have{$gene_id}";
		}else{
			$have{$gene_id}++;
		}
		$pos="$chr,$star,$end";
		$gene_id=~s/\(/\[/g;
		$gene_id=~s/\)/\]/g;
		$gene_id=~s/"//g;
		#$gene_id="$index\_$gene_id";
		$unit[8]="ID=$gene_id;";
		print GFF join("\t", @unit),"\n";
		$mRNA_id="NO";
	}
	if ($type eq "mRNA") {
		$gff_num++;
		my $gid;
		if($info=~m/ID=/){
			($mRNA_id, $gid)=&extract_id($info);
			if (!defined $gid ) {
				print "Warn: gff mRNA line does not exists Parent=!\n";
			}
		}elsif($info=~/mRNA\s+/){
			($mRNA_id, $gid) = $info=~/mRNA (\S+?) ; Gene (\S+?) ;/;
		}else{
			print "Gff Format is Something else, please check!!!\n" ;
			exit(1) ;
		}
		## convert
		$mRNA_id=~s/mRNA:/\.1/;
		$mRNA_id=~s/_mRNA/\.1/;
		$mRNA_id=~s/\|/_/;
		$mRNA_id=~s/\:/_/;
		$mRNA_id=~s/\s+/\_/g;
		$mRNA_id=~s/\%//g;

		if (defined $gene_id&&$gene_id ne $mRNA_id) {
			$gff{$mRNA_id}=$gene_id;
		}else{  ## no gene line, add
			$gff{$mRNA_id}=$mRNA_id;
			my @gene_unit = @unit;
			$gene_unit[2] = "gene";
			$gene_unit[8]="ID=$mRNA_id";
			$gff{$mRNA_id}=$mRNA_id;
			print GFF join("\t", @gene_unit),"\n";
		}
	
		$unit[8]="ID=$mRNA_id;Parent=$gff{$mRNA_id};";
		print GFF join("\t", @unit),"\n";
	}

	if ($type eq "CDS") {
		$hash_gene_CDS{$pos}++;
		next if($_=~/^NC.*transl_table=11/);
		if($mRNA_id eq "NO"){
			$mRNA_id=$gene_id;
		}
		$unit[8]="Parent=$mRNA_id;";
		print GFF join("\t", @unit),"\n";
	}
      if ($type eq "exon") {
                next if(!defined $mRNA_id);
                next if($_=~/gbkey=ncRNA/);
		
		next if($_=~/^NC.*transl_table=11/);
                if($mRNA_id eq "NO"){
                        $mRNA_id=$gene_id;
                }
                $unit[8]="Parent=$mRNA_id;";
                print GFF join("\t", @unit),"\n";
        }


}
close (IN) ;
close (GFF) ;
open (GFF1, "$outdir/Raw/$index.gff3.tmp") or die $!;
open (GFF2, ">$outdir/Raw/$index.gff3") or die $!;
while (<GFF1>){
	chomp;
        next if (m/^\#/) ;
        next if($_!~/\w+/);
        my @unit = split /\t/, $_;
       my ($chr, undef, $type, $star, $end, undef, $ori, $phase, $info) = @unit;
	next if($type eq "gene"&&!defined $hash_gene_CDS{"$chr,$star,$end"});
	print GFF2 join("\t", @unit),"\n";
}
close GFF1;
close GFF2;


#########Get Longest CDS and Pep file
open(SEQ,">$outdir/Raw/$index.genome.fa");
my $ina = Bio::SeqIO->new(-file => $fgenome, -format => 'fasta');
while(my $obj = $ina->next_seq()){
        my $id = $obj->id;
        my $seq = $obj->seq;
        print SEQ ">$id\n$seq\n";
}
close SEQ;
`LD_LIBRARY_PATH='' &&gffread $outdir/Raw/$index.gff3 -g $outdir/Raw/$index.genome.fa  -x $outdir/Raw/$index.cds.fa -y $outdir/Raw/$index.protein.fa`;
print "perl $Bin/get_longest_pepV2.pl -index $index  -pep  $outdir/Raw/$index.protein.fa -cds $outdir/Raw/$index.cds.fa -gff $outdir/Raw/$index.gff3 -od $outdir/Longest\n";
`perl $Bin/get_longest_pepV2.pl -index $index  -pep  $outdir/Raw/$index.protein.fa -cds $outdir/Raw/$index.cds.fa -gff $outdir/Raw/$index.gff3 -od $outdir/Longest`;
my $fcds="$outdir/Longest/$index.longest.cds.fa";
my $fpep="$outdir/Longest/$index.longest.protein.fa";
my $fcds_name = basename($fcds);
my $fpep_name = basename($fpep);
# ------------------------------------------------------------------
# cds pep check
# ------------------------------------------------------------------
my ($cds_num) =split/\s+/, `grep ">" $fcds|wc -l`;
my ($pep_num) =split/\s+/, `grep ">" $fpep|wc -l`;
my ($gene_ori)=split/\s+/, `awk '{if(\$3 =="gene")print }' $fgff |wc -l`;
print "cds num = $cds_num\n";
print "pep num = $pep_num\n";
if(($gene_ori-$cds_num)>3000){
	($gene_ori)=split/\s+/, `awk '{if(\$3 =="gene")print }' $fgff |grep "protein_coding"|grep -v "partial=true"|wc -l`;
}

print "gene num =$gene_ori\n";
if ($cds_num != $pep_num || $cds_num<100||($gene_ori-$cds_num)>3000) {
	print " $cds_num != $pep_num || $cds_num<100||($gene_ori-$cds_num)>3000 -> $index gff and genome are not consistent or the gff3 is not standard format,please check the longest cds sequence number  is consistent with the gene number!\n";
	exit(-1);
}
`perl $Bin//accordingIDgetGFF.pl  -prot $outdir/Longest/$index.id -combine $outdir/Raw/$index.gff3  -alter -o $outdir/Longest/$index.longest.gff3`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub extract_id {#
	my ($info)=@_;
	my ($mRNA_id, $gid);
	my @unit = split /;/, $info;

	#mRNA_id: 1: Name=(), 2: ID=()
	#gid: Parent=(), or not defined
	for (my $i=0; $i<@unit; $i++) {
		if ($unit[$i]=~/ID=/) {
			($mRNA_id)=$unit[$i]=~/ID=(\S+)/;
		}
		if($unit[$i]=~/Name=/){
			next if (defined $mRNA_id);
			($mRNA_id)=$unit[$i]=~/Name=(\S+)/;
		}
		if($unit[$i]=~/Parent=/){
			($gid)=$unit[$i]=~/Parent=(\S+)/;
		}
	}
	return ($mRNA_id, $gid);
}

sub AbsolutePath
{		#��ȡָ��Ŀ¼���ļ��ľ���·��
        my ($type,$input) = @_;

        my $return;
		$/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<zenghp\> 
	
Function: Check and Convert gff cds pep file,	

Usage:
  Options:
  -gff <file>        gff file, forced
  -genome <file>     genome file, forced
  -index <string>    index file,forced 
  -pcfg  <file>  the parameter config file 
  -od  <dir>         Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
