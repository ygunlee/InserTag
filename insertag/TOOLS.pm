package TOOLS;

use Bio::SearchIO;
use strict;
use warnings;

sub new {
	my $class = shift;
	my %param = @_;
	my $self = {};
	$self->{PRG} = undef;
	$self->{SUB} = undef;
	$self->{OPTIONS} = undef;
	foreach my $k (keys %param) {
		my $k1 = $k;
		$k1 = substr($k, 1) if($k =~ m/^-/);
		$self->{$k1} = $param{$k};
	}

	bless $self, ref($class) || $class;
	return $self;
}

sub program {
	my ($self, $value) = @_;
	$self->{PRG} = $value if($value);
	return $self->{PRG};
}

sub options {
	my ($self, $value) = @_;
	$self->{OPTIONS} = $value if($value);
	return $self->{OPTIONS};
}

sub run {
	print STDERR "you need implement your own run method";
}

package Blat;
use strict;
use Carp;
use English;
use Bio::SearchIO;
use Bio::SeqIO;
#use MISC;
use base qw(TOOLS);

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	$self->{PRG} = "blat" if(!$self->{PRG});
	$self->{OPTIONS} = "-tileSize=9 -stepSize=1 -out=blast" 
		if(!$self->{OPTIONS});
	return $self;
}

# return best hit for each query reads
sub run {
	my $self = shift;
	my %param = @_;
	croak "Missing TARGET parameter for $self->{PRG}" if(!$param{-TARGET});
	croak "Missing QUERY parameter for $self->{PRG}" if(!$param{-QUERY});
	croak "Missing OUTPUT parameter for $self->{PRG}" if(!$param{-OUTPUT});
	my $log = "1> /dev/null";
	system( join(" ", ($self->{PRG}, $param{-TARGET}, $param{-QUERY}, $param{-OUTPUT}, $self->{OPTIONS}, $log)));
	#return $output if($param{-FILE});
	#return _find_best_hit($output);
	#return find_best_blast($output) if(-e $output);
}

sub parseBlat {
	my ($file,$min_percent)=@_;
	my $in=Bio::SearchIO->new(-file=>$file,-format=>'blast');
	while(my $by_query=$in->next_result()){
	  while(my $by_db=$by_query->next_hit()){
	    while(my $by_hit=$by_db->next_hsp()){
		my $query_id = $by_query -> query_name();
		my $query_length = $by_query -> query_length();
		my $db_name = $by_db -> name();
		my $db_score = $by_db -> raw_score();
		my $db_sig = $by_db -> significance();
		my $hit_evalue = $by_hit -> evalue();
		my $hit_query_string = $by_hit -> query_string();
		my $hit_db_string = $by_hit -> hit_string();
		my $hit_homology_string = $by_hit -> homology_string();
		my $hit_homology_length = $by_hit -> length('total');
		my $hit_query_length = $by_hit -> length('query');
		my $hit_db_length = $by_hit -> length('hit');
		my $hit_db_strand = $by_hit -> strand('hit');
		my $hit_query_strand = $by_hit -> strand('query');
		my $hit_query_start = $by_hit -> start('query');
		my $hit_query_end = $by_hit -> end('query');
		my $hit_db_start = $by_hit -> start('hit');
		my $hit_db_end = $by_hit -> end('hit');
		my $hit_number = $by_db -> num_hsps();
		my @info=($query_id,$query_length,$db_name,$db_score,$db_sig,$hit_evalue,$hit_query_string,$hit_db_string,$hit_homology_string,$hit_homology_length,$hit_query_length,$hit_db_length,$hit_db_strand,$hit_query_strand,$hit_query_start,$hit_query_end,$hit_db_start,$hit_db_end,$hit_number);
		my $info=join("\t",@info);
		my $ratio=$hit_query_length/$query_length;
		if($db_sig == $hit_evalue && $ratio <= $min_percent){
		print $info,"\t",$ratio,"\n";
		}
	}}}}

sub maskBlat {
	my ($file,$seq,$min_overlap)=@_;
	my $in=Bio::SearchIO->new(-file=>$file,-format=>'blast');
	my $new=$seq;
	while(my $by_query=$in->next_result()){
	  while(my $by_db=$by_query->next_hit()){
	    while(my $by_hit=$by_db->next_hsp()){
		if($by_hit->length('query') >= $min_overlap){
		#print $by_hit->start('query') - 1 ,"\t",$by_hit->length('query'),"\t",length($new),"\n";
		substr($new,$by_hit->start('query') - 1,$by_hit->length('query'),"N"x$by_hit->length('query'));		
		}}}}
	return $new;
}

sub maskBlat2 {
	my ($file,$seq,$min_overlap)=@_;
	my $in=Bio::SearchIO->new(-file=>$file,-format=>'blast');
	my $new=$seq;
	my @start;
	my @start2;
	while(my $by_query=$in->next_result()){
	  while(my $by_db=$by_query->next_hit()){
	    while(my $by_hit=$by_db->next_hsp()){
		if($by_hit->length('query') >= $min_overlap){
		#print $by_hit->start('query') - 1 ,"\t",$by_hit->length('query'),"\t",length($new),"\n";
		substr($new,$by_hit->start('query') - 1,$by_hit->length('query'),"N"x$by_hit->length('query'));		
		push(@start,$by_hit->start('hit'));
		push(@start,$by_hit->end('hit'));
		push(@start2,$by_hit->start('query'));
		push(@start2,$by_hit->end('query'));
		}}}}
	#my @re=($new,\@start,\@start2);
	return ($new,\@start,\@start2);
}

sub multimaskBlat {
	my ($file1,$file2,$min_overlap)=@_;
	my $blatIn=Bio::SearchIO->new(-file=>$file1,-format=>'blast');
	my $fastaIn=Bio::SeqIO->new(-file=>$file2,-format=>'fasta');
	my %seq;
	while(my $seqIn=$fastaIn->next_seq()){
		$seq{$seqIn->id}=$seqIn->seq();
	}
	while(my $by_query=$blatIn->next_result()){
		while(my $by_db=$by_query->next_hit()){
			while(my $by_hit=$by_db->next_hsp()){
				if($by_hit->length('query') >= $min_overlap){
					print $by_query->query_accession(),"\n";
					print $seq{$by_query->query_accession()},"\n";
				}
			}
		}
	}
}
					
sub bestBlat { 
	my ($file,$min)=@_;
	my $in=Bio::SearchIO->new(-file=>$file,-foramt=>'blast');
	my @info;
	while(my $by_query=$in->next_result()){
	  while(my $by_db=$by_query->next_hit()){
	    while(my $by_hit=$by_db->next_hsp()){
		my $db_sig=$by_db->significance();
		my $hit_evalue=$by_hit->evalue();
		my $query_length=$by_query->query_length();
		my $hit_query_length=$by_hit->length('query');
		my $ratio = $hit_query_length/$query_length;	
		if($db_sig == $hit_evalue && $ratio >= $min){
			@info=($db_sig,$by_hit->start('hit'),
				$by_hit->end('hit'),
				$by_hit->start('query'),
				$by_hit->end('query'),
				$ratio*100);
		}
	}}}
	return @info;
}
	
1;


