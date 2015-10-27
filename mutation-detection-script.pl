##Author: Mustafa Adel
##mustafa.adel@aucegypt.edu

#!/usr/bin/perl
use strict;
use Bio::SearchIO;

print "\nUsage: $0 [input mutations list file] [input blast text report]\n\n";

unless(open(LIST, $ARGV[0])) {print STDERR "missing or inaccessible input mutations list\n";exit;}

unless(-e $ARGV[1]) {print STDERR "missing or inaccessible input blast text report\n";exit;}
my $in = new Bio::SearchIO(-format => 'blast', -file => $ARGV[1]);

#if(-e $ARGV[1].".tsv") {print STDERR "output file $ARGV[1].tsv already exists\n";exit;}
#unless(open(OUT,">".$ARGV[1].".tsv")) {print STDERR "inaccessible output file\n";exit;}

my %list; my ($gene,$pos,$mt);
while (<LIST>)
{
	chomp;
	$gene = (split(/\t/,$_))[0];
	$pos = (split(/\t/,$_))[1];
	$mt = (split(/\t/,$_))[2];
	$list{$gene}{$pos}{$mt} = '';
}

my $q = 0;
while(my $r = $in->next_result)
{
	$q++;
	my $h = 0;
	while (my $hit = $r->next_hit)
	{
		$h++;
		my $gene = $hit->name;			
		if (exists $list{$gene})
		{
			while(my $hsp = $hit->next_hsp)
			{
				for my $pos (keys %{$list{$gene}})
				{
					my ($p_start,$p_end); if($pos=~/-/){$p_start=(split("-",$pos))[0];$p_end=(split("-",$pos))[1];} else {$p_start=$pos;$p_end=$pos;}
					if ( ($p_start <= $hsp->end('hit')) && ($p_end >= $hsp->start('hit')) )
					{
						for my $mt (keys %{$list{$gene}{$pos}})
						{
								my $offset=$p_start-$hsp->start('hit');	my $length=$p_end-$p_start+1;
								if ($p_start < $hsp->start('hit')) {$offset = 0; $length= $length-($hsp->start('hit')-$p_start);}
								my $g = 0; my $g_l = 0; if ($hsp->gaps('hit') > 0)
								{
									my @g = $hsp->seq_inds('hit', 'gap');
									foreach (@g)
									{
										$g++ if $_ < $p_start;
										if (($_ >= $p_start) && ($_ <= $p_end) ) {$g_l++;};
									}
								}
								my $wt_s = substr($hsp->hit_string,$offset+$g,$length+$g_l);
								my $mt_s = substr($hsp->query_string,$offset+$g,$length+$g_l);
								my $hm_s = substr($hsp->homology_string,$offset+$g,$length+$g_l);
								unless ($wt_s eq $mt_s)
								{
									if ($mt eq $mt_s)
									{
										print $r->query_name," # ",$hit->name," | hsp: ",$hsp->rank," (",$hsp->start('hit'),"-",$hsp->end('hit'),")<-(",$p_start,"-",$p_end,") $wt_s=>$mt_s\n\n";
										my $x;
										if ( length($hsp->start('query')) >= length($hsp->start('hit')) )
										{
											$x=length($hsp->start('query'));
											print "Query: ", $hsp->start('query'), " ", $hsp->query_string, " ", $hsp->end('query'),"\n";
											print "       "; for (1..$x+1) { print " ";}; print $hsp->homology_string, " (", $hsp->length, ")\n";
											print "Sbjct: ", $hsp->start('hit'); for ( 1..($x+1 - length($hsp->start('hit')) )) { print " ";}; print $hsp->hit_string, " ", $hsp->end('hit'),"\n\n";
										}
										else
										{
											$x=length($hsp->start('hit'));
											print "Query: ", $hsp->start('query'); for ( 1..($x+1 - length($hsp->start('query')) )) { print " ";}; print $hsp->query_string, " ", $hsp->end('query'),"\n";
											print "       "; for (1..$x+1) { print " ";}; print $hsp->homology_string, " (", $hsp->length, ")\n";
											print "Sbjct: ", $hsp->start('hit'), " ", $hsp->hit_string, " ", $hsp->end('hit'),"\n\n";
										}
									}
								}
						}
					}
					else
					{
						next;
					}
				}
			}
		}		
		else
		{
			next;
		}
	}
}
