#!/usr/bin/perl -w
# Copyright (c) SXB 2015/6/28
# Writer:         SXB <SXBSXB>
# Program Date:   2015/6/28.
# Modifier:       SXB <SXBSXB>
# Last Modified:  2015/6/28.

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Bio::SeqIO;


my $programe_dir=basename($0);
my $path=dirname($0);

my $ver    = "1.0";
my $Writer = "SXB <SXB\SXB>";
my $Data   = "2015/6/28";
my $BEGIN=time();

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($prot,$out,$format,$combine,$alter,$geneIDformat,$mRNAIDformat);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"prot:s"=>\$prot,
			"alter"=>\$alter,
			"geneIDformat!"=>\$geneIDformat,
			"mRNAIDformat!"=>\$mRNAIDformat,
			"format!"=>\$format,
			"combine:s"=>\$combine,
			) || &help;
&help unless ($prot && $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: ......
    Usage:
        -prot              prot file file you want to get
        -combine         combined file; 
        -alter            for alterspliced gff file
        -o               out file
        -h               Help document
	Usage End.
	exit;
}
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
################

my @genegff;
my %hash;
open(IN,"$prot");
while(<IN>){
	chomp;
	my @tmp=split/\s+/;
        $hash{$tmp[1]}=1;
	$hash{$tmp[0]}=1;
	}
close IN;



open (OUT, ">$out")|| die "cannot open $out:$!";
open (IN2, "<$combine")|| die "cannot open $combine:$!";
my @merge;
while(<IN2>){
	chomp;
	next if(/^\#/);
	next if(/^$/);
	my @tmp=split/\t+/;
	if($tmp[2]eq "gene"){
		my ($id)=$tmp[8]=~/ID=([^;]+)/;
		if(defined $hash{$id}){
			print OUT "$_\n";
		}
	}elsif($tmp[2] eq "mRNA"){
		my ($id)=$tmp[8]=~/ID=([^;]+)/;
                if(defined $hash{$id}){
                        print OUT "$_\n";
                }
	}else{
		my ($id)=$tmp[8]=~/Parent=([^;]+)/;
                if(defined $hash{$id}){
                        print OUT "$_\n";
                }
	}
		
}		


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime # &Runtime($BEGIN);
{
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe_dir elapsed time : [",&sub_time($t),"]\n";
}
sub sub_time
{
	my ($T)=@_;chomp $T;
	my $s=0;my $m=0;my $h=0;
	if ($T>=3600) {
		my $h=int ($T/3600);
		my $a=$T%3600;
		if ($a>=60) {
			my $m=int($a/60);
			$s=$a%60;
			$T=$h."h\-".$m."m\-".$s."s";
		}else{
			$T=$h."h-"."0m\-".$a."s";
		}
	}else{
		if ($T>=60) {
			my $m=int($T/60);
			$s=$T%60;
			$T=$m."m\-".$s."s";
		}else{
			$T=$T."s";
		}
	}
	return ($T);
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
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
        elsif($type eq 'f_dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/";
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
