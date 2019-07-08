package BLAT;

use INSERTAG;
use TOOLS;
use Bio::SeqIO;
use Bio::Seq;

sub runblat
{
	my ($query, $subject,$prefix) = @_;
	my $aligner=Blat->new(-PRG=>"blat",-OPTIONS=>" -tileSize=15 -stepSize=1 -out=blast8 -minScore=15 -noHead ");
	
	my $query_io=Bio::SeqIO->new(-file=>">$prefix.query.fa",-format=>'fasta');
	my $subject_io=Bio::SeqIO->new(-file=>">$prefix.subject.fa",-format=>'fasta');
	
	my $query_seq=Bio::Seq->new(-seq=>$query,-id=>"QUERY_$prefix");
	my $subject_seq=Bio::Seq->new(-seq=>$subject,-id=>"SUBJECT_$prefix");
	
	$query_io->write_seq($query_seq);
	$subject_io->write_seq($subject_seq);

	$aligner->run(	-TARGET=>"$prefix.subject.fa",
			-QUERY=>"$prefix.query.fa",
			-OUTPUT=>"$prefix.blat8");

	open(OUT,"$prefix.blat8");
	my $i=1;
	my %blast;
	while(my $hit=<OUT>){
		my @hit=INSERTAG::getBLAST8($hit);
		#(0)identity (1)aln_length (2)query_start (3)query_end (4)target_start (5)target_end (6)target_name
		$blast{$i}=\@hit;
		$i++;
	}
	
	system("rm $prefix.subject.fa $prefix.query.fa $prefix.blat8");
		
	return \%blast;
}


sub runblat_psl
{
	my ($query, $subject,$prefix) = @_;
	my $aligner=Blat->new(-PRG=>"blat",-OPTIONS=>" -tileSize=15 -stepSize=1 -minScore=15 -noHead ");
	my $query_io=Bio::SeqIO->new(-file=>">$prefix.query.fa",-format=>'fasta');
	my $subject_io=Bio::SeqIO->new(-file=>">$prefix.subject.fa",-format=>'fasta');
	my $query_seq=Bio::Seq->new(-seq=>$query,-id=>"QUERY_$prefix");
	my $subject_seq=Bio::Seq->new(-seq=>$subject,-id=>"SUBJECT_$prefix");
	
	$query_io->write_seq($query_seq);
	$subject_io->write_seq($subject_seq);

	$aligner->run(	-TARGET=>"$prefix.subject.fa",
			-QUERY=>"$prefix.query.fa",
			-OUTPUT=>"$prefix.psl");

	my $i=1;
	my %blast;
	open(PSL,"$prefix.psl");
	while(<PSL>){
		if($_=~/^(\d+)/){
			chomp $_;
			my $line=$_;
			$blast{$i}=\$line;
			$i++;
		}
	}
	
	system("rm $prefix.subject.fa $prefix.query.fa $prefix.psl");
	return \%blast;
		
}

sub runblat_contig
{
	my ($query,$subject,$prefix) = @_;
	my $psl=BLAT::runblat_psl($query,$subject,$prefix);
	my $pslhit=0;
	if($psl->{1}){
		my @hit=split(/\t/,${$psl->{1}});
		$pslhit=\@hit;
	}
	return $pslhit;
}


1; 

