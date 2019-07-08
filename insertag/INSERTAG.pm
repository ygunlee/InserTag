package INSERTAG;

##getOabed 
##localAssembly 
##getAFG
##parseAFG
##readVelvet

use warnings;
use strict;

use List::Util qw( min max sum);
use Bio::SeqIO;
use Bio::SearchIO;

sub getAFG  # copy from AmosLib.pm by Mihai pop of TIGR
{
    my $file = shift;
    my $level = 0;
    my $block = "";
    while (<$file>){
      if (/^\s*\{/){
          $level++;
      }
      if (/^\s*\}/){
          $level--;
      }
      $block .= $_;
      if ($level == 0){
          last;
      }
    }
    if ($level != 0){
      die ("end of file reached before end of block\n");
    }
    if ($block ne ""){
      return $block;
    } else {
      return undef;
    }
}

sub getContig
{
 my $record=shift;
 my @lines=split(/\n/,$record);
 my @contig=split(/\t/,$lines[2]);
 return @contig;
} 

sub getEach
{
 my $record = shift;
 my @lines=split(/\n/,$record);
 my $last= $#lines - 2;
 my %store;
 my %offset;
 my $even=0;
 my $odd=0;
 my @qual;
 my @box;
 my %class=(0=>0,1=>0,2=>0);
 my $decoy=0;
 my %start;
 my %end;
 my %split;
 for my $i (3 .. $last){
	my @read = split(/\t/,$lines[$i]);
	my @str = split(/,/,$read[15]);
	my $direction=0;
	if(($str[1] - $str[0]) > 0){$direction=1;} #0->100
	else{$direction=-1}; #100->0
	if($read[6] % 2 == 0 && $read[8]=~/hs37d5/){$decoy++;}
	
	$store{$read[6]}=$direction;
	$offset{$read[6]}=[$read[14],$str[1]-$str[0],$read[5]]; #offset,complement,cigar
	push(@qual,$read[4]);
	push(@box,$read[1]);
	push(@box,$read[2]);
	$start{$read[6]}=$read[1];
	$end{$read[6]}=$read[2];

	if($read[7]==0){$even++;}
	else{$odd++;}
	$class{$read[16]}++;
	if($read[16] == 3 && $read[7] == 0){ ##USAGE OF SPLIT READ
		$split{$read[6]}=1;
		$start{$read[6]}=$read[9];
		$end{$read[6]}=$read[10];
		${$offset{$read[6]}}[2]=$read[11];
	}else{$split{$read[6]}=0;}	
 }
 my $min= min @box;
 my $max= max @box;
 my $mean= (sum @qual)/(scalar @qual);

 return (\%store,$even,$odd,$min,$max,$mean,\%class,\%offset,$decoy,\%start,\%end,\%split);
}


sub getReads
{
 my $record = shift;
 my @lines=split(/\n/,$record);
 my $last= $#lines - 2;

 my @contig = split (/\t/,$lines[2]);
 my $seq = $contig[4];
 my $total = $contig[2];
 my $oh = $contig[5];
 my $ac = $contig[6];
 my $size = $contig[3];

 my @box;
 my @qual;
 my $chr;
 my $str;
 my %count;
 for my $i (3 .. $last){
	my @read = split(/\t/,$lines[$i]);
	push(@box,$read[1]);
	push(@box,$read[2]);
	push(@qual,$read[4]);
	$chr=$read[0];
	$str=$read[3];
	my $remain=$read[6]%2;
	$count{($read[6]-$remain)/2}++;
 }
 my $min= min @box;
 my $max= max @box;
 my $mean= (sum @qual)/($last-2);
 my $count= scalar keys (%count);

 return ($chr,$min,$max,$str,$mean,$total,$oh,$ac,$size,$seq,$count);
}

sub maskContig
{
 my ($seq,$qstart,$qend,$i) = @_; 
 my $start;
 my $size;

 if($qstart>$qend){
  $start=$qend;
  $size=$qstart-$qend+1;
 }else{
  $start=$qstart;
  $size=$qend-$qstart+1;
 }

 substr($seq,$start-1,$size,$i x $size);
 return $seq;
}
 
sub countMask
{
 my $seq = @_;
 #my @count = ($seq =~ m/(\d+)/g);
 #my $num = scalar (@count);
 $seq=~m/\d/g;
 return $seq;
}


sub getBLAST8
{
 my $record = shift;
 my @line = split(/\t/,$record);
 return ($line[2],$line[3],$line[6],$line[7],$line[8],$line[9],$line[10],$line[1]);
}

#identitiy, alignment length, 
#contig start, contig end 
#ref start, ref end 
#alignement e-value 
 

sub parseAFG  # copy from AmosLib.pm by Mihai pop of TIGR
{
    my $record = shift;
    my @lines = split('\n', $record);
    my $type;
    my %fields;
    my @recs;
    # get record type
    $lines[0] =~ /\{(\w+)/;
    if (! defined $1){
        die ("Weird start of record: $record\n");
    }
    $type = $1;
    if ($lines[$#lines] !~ /^\s*\}/){
      die ("Weird end of record: $record\n");
    }
    my $level = 0;
    my $fieldname;
    for (my $i = 1; $i < $#lines; $i++){
      if ($lines[$i] =~ /^(\w+):(.+)$/){   # simple field
          $fields{$1} = $2;
      } # simple field
      if ($lines[$i] =~ /^(\w+):$/){ # complex field
          $fieldname = $1;
          $fields{$fieldname} = "";
          $i++;
          while ($i < $#lines && ($lines[$i] !~ /^\.$/)){
            $fields{$fieldname} .= "$lines[$i]\n";
            $i++;
          }
      } # complex field
      if ($lines[$i] =~ /^\s*\{/){ # subrecord
          my $level = 1;
          my $thisrec = ++$#recs;
          $recs[$thisrec] .= "$lines[$i]\n";
          $i++;
          while ($level > 0 && $i < $#lines){
            if ($lines[$i] =~ /^\s*\{/){
                $level++;
            }
            if ($lines[$i] =~ /^\s*\}/){
                $level--;
            }
            $recs[$thisrec] .= "$lines[$i]\n";
            if ($level == 0){
                last;
            } else {
                $i++;
            }
          }
          if ($level != 0){
            die ("Error parsing sub_record in: $record\n");
          }
      } # subrecord
    } # for $i...
    return ($type, \%fields, \@recs);
}

sub getOabed {	
	my $str = shift;
	chomp $str;
	my %oabed;
	my @oabed = split (/\t/,$str);
	$oabed{ACHR}=$oabed[0];
	$oabed{ASTART}=$oabed[1];
	$oabed{AEND}=$oabed[2];
	$oabed{AID}=$oabed[3];
	$oabed{ASEQ}=$oabed[4];
	$oabed{AQUAL}=$oabed[5];
	$oabed{AMQ}=$oabed[6];
	$oabed{ASTR}=$oabed[7];
	$oabed{ACIGAR}=$oabed[8];
	$oabed{OCHR}=$oabed[9];
	$oabed{OSTART}=$oabed[10];
	$oabed{OEND}=$oabed[11];
	$oabed{OSEQ}=$oabed[12];
	$oabed{OQUAL}=$oabed[13];
	$oabed{OMQ}=$oabed[14];
	$oabed{OSTR}=$oabed[15];
	$oabed{OSIZE}=$oabed[16];
	$oabed{OCIGAR}=$oabed[17];
	$oabed{OCLASS}=$oabed[18];
	$oabed{OCLUSTER}=$oabed[19];
	return \%oabed;
}

sub localAssembly {
	my ($ref,$pid,$sup,$ind,$kmer,$inslen,$minlen) = @_; 
	#input1 = oabed_stored_hash (key=odd:anchor_part // even:overhang_part) -> value=oabed_hash
	#input2 = cluster_id 
	#input3 = supporting_reads_number
	#input4 = output_file_name
	#input5 = kmer
	#input6 = library_insert_length
	#input7 = contig_minmum_length
	my %store = %$ref; 
	if (scalar (keys %store) >= $sup){
		my $out1=Bio::SeqIO->new(
			-file=>">$ind.$pid.anchor.fa", 
			-format=>'fasta');
		my $out2=Bio::SeqIO->new(
			-file=>">$ind.$pid.overhang.fa",
			-format=>'fasta');

		foreach my $key (sort {$a<=>$b} keys %store){
			if($key %2 == 1) { ##Anchor_part
			my $seq1 = Bio::Seq->new(
				-seq=>$store{$key}->{ASEQ},
				-id=>$key);
				$out1->write_seq($seq1);
			}else{ ##Overhang_part
			my $seq2 = Bio::Seq->new(
				-seq=>$store{$key}->{OSEQ},
				-id=>$key);
				$out2->write_seq($seq2);
			}	
		}
	}
	system ("velveth $ind.Assem $kmer -shortPaired -separate -fsata $ind.$pid.anchor.fa $ind.$pid.overhang.fa > $ind.$pid.log");
	system ("velvetg $ind.Assem -ins_length $inslen -exp_cov auto -cov_cutoff auto -min_contig_lgth $minlen -read_trkg yes -amos_file yes > $ind.$pid.log");
	#system ("rm $ind.$pid.log");
	#system ("rm $ind.$pid.anchor.fa");
	#system ("rm $ind.$pid.overhang.fa");
}

1;

