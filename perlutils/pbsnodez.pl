#!/usr/bin/perl
use Time::Local;
use Term::ANSIColor;
use Data::Dumper;


$output = `pbsnodes`;
$zstatout = `zstat.pl batch`;
foreach my $line (split '\n', $zstatout)
{
	chomp $line;
	$line =~ s/[\000-\037]\[(\d|;)+m//g;;
	$line =~ s/^\s+//;
	$line =~ s/\s+/|/g;
	my @z = split(/\|/, $line);
	$zstat{$z[0]} = [@z];
}

@nodes = split "\n\n", $output;

foreach  $n (@nodes)
{
	my @lines = split "\n", $n;
	$node = $lines[0];
	chomp $node;
	for $l (@lines)
	{
		if( $l =~ /\s*(\S*)\s+\=\s+(.+)$/)
		{
			$nodeprops{$node}{$1} = $2;

		}
	}

}




foreach $k (sort keys %nodeprops)
{
	my %nodejobs;
	#my %detailedProps = map { $_ =~ /(.+)\=(.+)/,  2} split /\,/, $nodeprops{$k}{status};
	my %detailedProps = split /[,=]/,  $nodeprops{$k}{status};
	
	#print Dumper(%detailedProps);
	$color = $nodeprops{$k}{jobs} =~ /s+/ ? "red" : "green";
	print color($color) . "$k    " . color("reset");
	print "$nodeprops{$k}{np}core\t";
	$nodeprops{$k}{jobs} =~ s/\-\d+.hpc-pbs.hpcc.usc.edu//g;
	$nodeprops{$k}{jobs} =~ s/.hpc-pbs.hpcc.usc.edu//g;
	$nodeprops{$k}{jobs} =~ s/\d+\///g;
	$nodeprops{$k}{jobs} =~ s/ //g;
	$nodejobs{$_}++ for split ",",$nodeprops{$k}{jobs};
	print "[" . color("blue") . "$_" . color("reset") . ":$nodejobs{$_}cores:" . color("yellow") . "$zstat{$_}[2]" . color("reset") . "]"  for  (keys %nodejobs);
	print color("green") . "FREE" .  color("reset") if !%nodejobs && $nodeprops{$k}{state} ne "down";
	print color("red") . "NODE IS OFFLINE" .  color("reset") if !%nodejobs && $nodeprops{$k}{state} eq "down";
	print color("magenta") . "\tload=$detailedProps{loadave} Node is idle!" . color("reset") if %nodejobs && $detailedProps{loadave} < 0.02;  
	print color("reset") . "\tload=$detailedProps{loadave}" . color("reset") if %nodejobs && $detailedProps{loadave} >= 0.02;  
	print "\n";



}
