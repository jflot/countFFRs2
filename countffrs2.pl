#!/usr/bin/perl --
#Copyright (c) 2015, JF Flot <jflot@gwdg.de>
#
#Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is
#hereby granted, provided that the above copyright notice and this permission notice appear in all copies.
#
#THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
#INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
#USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

warnings;
use strict;
use Getopt::Std;
use List::Util 'shuffle';

my %opts;
getopts("i:r:vfglR", \%opts);
unless (defined($opts{i})) {print "Usage: perl countffrs2.pl -i <input file in sequential FASTA format>  -r <number of replicate samplings (default is 1)> -v <verbiose> -f <save test dataset as FASTA> -g <save test dataset as gzipped FASTA> -l generate list files -R <generate Roehl output file for Network>\nThe name of each sequence in the FASTA should be of the form Txxx_yyy, where xxx is the species ID and yyy is the sequence ID.\n"; exit};
unless (defined($opts{r})) {$opts{r}=1};
#processing alignment
open(DATA, $opts{i}) || die "Error: $opts{i} cannot be opened!\n";
my@align=<DATA>;

    #remove all comments and endline characters
for (my $i=0; $i<=$#align; $i++)
	{
	$align[$i]=~ s/'//g;	
	chomp($align[$i]);
	$/ = "\r";
	chomp($align[$i]);
	$/ = "\n";
	if (substr($align[$i],0,1) eq ';') {splice(@align,$i, 1); $i=$i-1}
	};
    #make a list of species
#print "@align\n\n";

    #makes a list of sequences
my @allsequences;
for (my $i=0; $i<($#align+1)/2; $i++)
    {my @t=split('_', $align[2*$i]);
     $t[0]=~ s/>T//g;	
     $allsequences[$i][0]=$t[0];
     $allsequences[$i][1]=$t[1];
     $allsequences[$i][2]=$align[2*$i+1]};

 #for debugging
    #for (my $i=0; $i<=$#allsequences; $i++) {print "Here: $allsequences[$i][2]\n"};

    #repeat procedure several r times to make stats about the number of FFRs
for (my $k=0; $k<$opts{r}; $k++) {	
my $fk=sprintf("%03d",$k);
print "Replicate $k:\n";
# rearrange the simulated sequences into a table of genotypes
my@table;
my$speciesnumber=1+$allsequences[($#align-1)/2][0];
for (my $s=1; $s<$speciesnumber; $s++) { #print "Species $s\n";
	my@sequences=@allsequences;
	for (my $i=0; $i<=$#sequences; $i++) { if ($sequences[$i][0]!=$s) {splice(@sequences,$i,1); $i=$i-1}}
;
if ($opts{v}) { print "There are ",1+$#sequences, " sequences for species $s; ";
if (0 == $#sequences % 2) {print "this is an odd number, please check input alignment."; exit}}
else {if (0 == $#sequences % 2) {print "There is an odd number of sequences for one species, please check input alignment."; exit}}

 #for debugging
    #for (my $i=0;  $i<=$#sequences; $i++) {print "$sequences[$i][1]\n"};

    

	    #for each simulated species, form diploid individuals
    my @shuffled_indexes; my @shuffled_sequences; my @t;
    @shuffled_indexes=shuffle(0..$#sequences); #print "shuffled: @shuffled_indexes\n";
    @shuffled_sequences=@sequences[@shuffled_indexes];

 #for debugging
   #for (my $i=0;  $i<=$#sequences; $i++) {print "$shuffled_sequences[$i][1]\n"};


my$nbindiv=(1+$#sequences)/2; if ($opts{v}) {print "I will form $nbindiv individuals.\n"};
		     
 	for (my $i=0; $i<$nbindiv; $i++) {my $indiv = sprintf("%03d",$i+1);
         $t[2*$i][0]="T".$s."indiv".$indiv."a"; $t[2*$i][1]=$shuffled_sequences[2*$i][2];
         $t[1+2*$i][0]="T".$s."indiv".$indiv."b"; $t[1+2*$i][1]=$shuffled_sequences[1+2*$i][2]}
 	   
@table=(@table,@t);
if ($opts{v}) {print "There are now ",1+$#table," sequences in the table.\n" }};	
	
	    #for debugging
    #for (my $i=0; $i<$nbindiv; $i++) {print "!$table[2*$i][1]\n$table[1+2*$i][1]\n"};
    
#write output dataset as FASTA if requested
	if ($opts{f} or $opts{g}) {
	open (F, "> $opts{i}.ffrs".$fk.".fasta") || die "Cannot write in FASTA file!";
	    for (my $i=0; $i<($#table+1)/2 ; $i++) {print F ">".$table[2*$i][0]."\n".$table[2*$i][1]."\n>".$table[1+2*$i][0]."\n".$table[1+2*$i][1]."\n"};
	    close F};
    if ($opts{g} and not $opts{f}) {system("gzip -v9 $opts{i}.ffrs".$fk.".fasta") && die "unable to gzip output FASTA\n"} ;
    if ($opts{g} and $opts{f}) {system("gzip -v9 -c $opts{i}.ffrs".$fk.".fasta > $opts{i}.ffrs".$fk.".fasta.gz") && die "unable to gzip output FASTA\n"} ;

	
	
	#all sequences are simulated so have supposedly the same length, but better check
    my $maxseqlength=length($table[0][1]);
	for (my$i=0; $i<=$#table; $i++)
		{my $l=length($table[$i][1]); # print "$l\n";
                 unless ($l==$maxseqlength) {print "Not all sequences in the input dataset $opts{i} have the same length!"; exit}};
	
	#make hash of haplotypes
	my %haplotypes;
	for (my$i=0; $i<=$#table; $i++) {push(@{$haplotypes{$table[$i][1]}}, $table[$i][0])};
	my $number_of_haplotypes = keys %haplotypes;
	
	if ($opts{v}) {
	print "There are $number_of_haplotypes haplotypes in your dataset of ", (1+$#table)/2," individuals.\n";
	}
	
	#make table of haplotypes
	my @haplotable; my @haplotable2;
	my$count=0;
	foreach my$key (sort {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}} keys %haplotypes) 
			{@{$haplotable[$count]}=@{$haplotypes{$key}}; $haplotable2[$count]=$key; $count++};
	
	#display list of haplotypes
	
	if ($opts{v}) {
	for (my$i=0; $i<$number_of_haplotypes; $i++)
			{print "Haplotype ", $i+1," is found ", $#{$haplotable[$i]}+1, " time(s) in the alignment, under the name(s): ", join(', ', @{$haplotable[$i]}), ".\n"};
	}
	
	#make heterozygotes table
	my @heterozygotes;
	for (my$i=0; $i<$#table; $i++)
	    {my $previous = $table[$i][0]; my $next = $table[$i+1][0]; my $endprevious=chop($previous); my 
	$endnext=chop($next);
	if (((($endprevious eq 'a') and ($endnext eq 'b')) or (($endprevious eq 'b') and ($endnext eq 'c')) or 
	(($endprevious eq 'c') and ($endnext eq 'd')) or (($endprevious eq 'd') and ($endnext eq 'e')) or 
	(($endprevious eq 'e') and ($endnext eq 'f'))) and (($previous eq $next) and ($table[$i][1] ne 
	$table[$i+1][1]))) {push(@heterozygotes, $previous) unless grep(/$previous/,@heterozygotes)}};
	
	#display list of heterozygotes
	
	if ($opts{v}) {
	print "A total of ", scalar(@heterozygotes), " individuals contain more than one haplotype: ", join(', ',@heterozygotes),".\n";
	}
	
	#make hash of connections
	my %connect;
	foreach my$indiv (@heterozygotes) 
		{my $ka; my $kb; my $kc; my $kd; my $ke; my $kf; for (my$i=0; $i<=$#haplotable; $i++)
			{if (grep(/$indiv+a/, @{$haplotable[$i]})) {$ka=$i+1}; if (grep(/$indiv+b/, 
	@{$haplotable[$i]})) {$kb=$i+1};
			 if (grep(/$indiv+c/, @{$haplotable[$i]})) {$kc=$i+1}; if (grep(/$indiv+d/, 
	@{$haplotable[$i]})) {$kd=$i+1}
			 if (grep(/$indiv+e/, @{$haplotable[$i]})) {$ke=$i+1}; if (grep(/$indiv+f/, 
	@{$haplotable[$i]})) {$kf=$i+1}};
	
			($ka,$kb,$kc,$kd,$ke,$kf)=sort {$a<=>$b} ($ka,$kb,$kc,$kd,$ke,$kf);
			$connect{"$ke and $kf"}{$indiv}='OK';
			if ($kd) {$connect{"$kd and $ke"}{$indiv}='OK'; $connect{"$kd and $kf"}{$indiv}='OK'};
			if ($kc) {$connect{"$kc and $kd"}{$indiv}='OK'; $connect{"$kc and $ke"}{$indiv}='OK'; 
	$connect{"$kc and $kf"}{$indiv}='OK'};
			if ($kb) {$connect{"$kb and $kc"}{$indiv}='OK'; $connect{"$kb and $kd"}{$indiv}='OK'; 
	$connect{"$kb and $ke"}{$indiv}='OK'; $connect{"$kb and $kf"}{$indiv}='OK'};
			if ($ka) {$connect{"$ka and $kb"}{$indiv}='OK'; $connect{"$ka and $kc"}{$indiv}='OK'; 
	$connect{"$ka and $kd"}{$indiv}='OK'; $connect{"$ka and $ke"}{$indiv}='OK'; $connect{"$ka and 
	$kf"}{$indiv}}};
	
	#remove spurious connections between one haplotype and itself (useful when some sites have been removed from the alignment, in which case individuals may have several times the same haplotype)
	for (my$i=0; $i<=$#haplotable; $i++) {delete($connect{"$i and $i"})};
	
	#display list of connections
	
	if ($opts{v}) {
	foreach my$key (sort {scalar(keys %{$connect{$b}}) <=> scalar(keys %{$connect{$a}})} keys %connect)
	{my @connectingindividuals = sort {$#{$connect{$a}} <=>$#{$connect{$b}}} keys %{$connect{$key}};
	print "Haplotypes $key co-occur in ", scalar(@connectingindividuals), " individual(s), namely: ", 
	join(', ', sort @connectingindividuals),".\n"};
	}
	
	#build FFRs
	foreach my$hetero(@heterozygotes) {my $keya; my $keyb; my $keyc; my $keyd; my $keye; my $keyf; my 
	@reunion;
	while (my ($k, $v) = each (%haplotypes)) {if (grep(/$hetero+a/, @{$v})) {$keya=$k}; if 
	(grep(/$hetero+b/, @{$v})) {$keyb=$k};
	if (grep(/$hetero+c/, @{$v})) {$keyc=$k}; if (grep(/$hetero+d/, @{$v})) {$keyd=$k}; if 
	(grep(/$hetero+e/, @{$v})) {$keye=$k};
	 if (grep(/$hetero+f/, @{$v})) {$keyf=$k}};
	if (($keyb) and ($keya ne $keyb)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyb}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyb}) unless (($keyb eq $keyc) or ($keyb eq $keyd) 
	or ($keyb eq $keye) or ($keyb eq $keyf))};
	if (($keyc) and ($keya ne $keyc)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyc}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyc}) unless (($keyc eq $keyd) or ($keyc eq $keye) 
	or ($keyc eq $keyf))}; 
	if (($keyd) and ($keya ne $keyd)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyd}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyd}) unless (($keyd eq $keye) or ($keyd eq 
	$keyf))}; 
	if (($keye) and ($keya ne $keye)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keye}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keye}) unless ($keye eq $keyf)};
	if (($keyf) and ($keya ne $keyf)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyf}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyf})};
	};
	
	#remove duplicate sequence names in haplotypes hash
	foreach my $k (keys %haplotypes) {my %hash = map { $_, 1 } @{$haplotypes{$k}}; @{$haplotypes{$k}}= keys 
	%hash};
	
	
	#sort sequence names in FFRs
	foreach my$key (keys %haplotypes) {@{$haplotypes{$key}} = sort @{$haplotypes{$key}} };
	
	
	#display FFRs in terms of sequences
	
	if ($opts{v}) {
	if (scalar(keys %haplotypes) == 1) {print "There is a single FFR in your dataset.\n"}
	else {print "There are ", scalar(keys %haplotypes), " FFRs in your dataset.\n"}
	};
	
	if ($opts{v}) {
	foreach my$key (sort {$#{$haplotypes{$b}} <=> $#{$haplotypes{$a}}} keys %haplotypes) {print "One FFR contains the following ", scalar(@{$haplotypes{$key}})," sequence(s): ", join(', ',@{$haplotypes{$key}},".\n")}; 
	}
	my $total_number_sequences=0;
	foreach my$key (keys %haplotypes) {$total_number_sequences+= scalar(@{$haplotypes{$key}})};
	if ($opts{v}) {
	print "Total number of sequences in the dataset: $total_number_sequences.\n";
	}
	#display FFRs in terms of haplotypes
	if ($opts{v}) {print "Or, in terms of haplotypes:\n"};
	my %haplonumbers;
	for (my $i=0; $i<=$#haplotable; $i++) {foreach my $seq (@{$haplotable[$i]}) 
	{$haplonumbers{$seq}=$i+1}};
	
	my %haplotypes2;
	foreach my$key (keys %haplotypes) {foreach my $seq (@{$haplotypes{$key}}) 
	{$haplotypes2{$haplonumbers{$seq}}=$key}};
	
	my %haplotypes3;
	foreach my$key (keys %haplotypes2) {push(@{$haplotypes3{$haplotypes2{$key}}},,$key)};
	
	if ($opts{v}) {
	foreach my$key (sort {$#{$haplotypes3{$b}} <=>$#{$haplotypes3{$a}}} keys %haplotypes3) 
	{print "One FFR contains the following ", scalar(@{$haplotypes3{$key}})," haplotype(s): ", join(', ', 
	sort {$a <=> $b} @{$haplotypes3{$key}}),".\n"}; 
	}
	
	#display FFRs in terms of individuals	
	if ($opts{v}) {
	print "Or, in terms of individuals:\n";
	}
	
	foreach my$key (keys %haplotypes) {
	for (my$i=0; $i<$#{$haplotypes{$key}}; $i++)
	    {my $previous = ${$haplotypes{$key}}[$i]; my $next = ${$haplotypes{$key}}[$i+1]; my $nextnext = 
	${$haplotypes{$key}}[$i+2];
	my $nextnextnext = ${$haplotypes{$key}}[$i+3]; my $nextnextnextnext = ${$haplotypes{$key}}[$i+4]; 
	my $nextnextnextnextnext = ${$haplotypes{$key}}[$i+5];
	my $endprevious=chop($previous); my $endnext=chop($next); my $endnextnext=chop($nextnext); my 
	$endnextnextnext=chop($nextnextnext);
	my $endnextnextnextnext=chop($nextnextnextnext); my 
	$endnextnextnextnextnext=chop($nextnextnextnextnext);
	if (($endprevious eq 'a') and ($endnext eq 'b') and ($previous eq $next)) 
		{${$haplotypes{$key}}[$i]=$previous; splice (@{$haplotypes{$key}},$i+1,1);
		if (($endnextnext eq 'c') and ($previous eq $nextnext)) 
			{splice (@{$haplotypes{$key}},$i+1,1);
			if (($endnextnextnext eq 'd') and ($previous eq $nextnextnext)) 
				{splice (@{$haplotypes{$key}},$i+1,1);
				if (($endnextnextnextnext eq 'e') and ($previous eq $nextnextnextnext))
					{splice (@{$haplotypes{$key}},$i+1,1);
					if (($endnextnextnextnextnext eq 'e') and ($previous eq 
	$nextnextnextnextnext))
						{splice (@{$haplotypes{$key}},$i+1,1);
	}}}}}}};
	
	if ($opts{v}) {
	foreach my$key (sort {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}} keys %haplotypes) {print "One FFR contains the following ", scalar(@{$haplotypes{$key}})," individual(s): ", (join(', ',@{$haplotypes{$key}})), "\n"};
	} 
    else { foreach my$key (sort {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}} keys %haplotypes) {print join(' ',@{$haplotypes{$key}}).";"}; print "\n"};
	my $total_number_individuals=0;
	foreach my$key (keys %haplotypes) {$total_number_individuals+= scalar(@{$haplotypes{$key}})};
	
	if ($opts{v}) {
	print "Total number of individuals in the dataset: $total_number_individuals.\n\n";
	};
	
	if ($opts{l}) {
	#generate list files as output
	if (!open (F, "> $opts{i}.ffrs".$fk.".list")) {die "Writing Error : $!"};
	my $ffr=1;
	my %listfile;
	foreach my$key (sort  {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}}  keys %haplotypes) {foreach my$ind (@{$haplotypes{$key}}) {$listfile{$ind}=$ffr}; $ffr++}; 
	foreach my$key (sort keys %listfile) {print F "$key\t$listfile{$key}\n"};
	close F};

	if ($opts{R}) {
	#generate Roehl output

	if (!open (F, ">$opts{i}.ffrs".$fk.".rdf")){die "Writing error: $!"};
	my @haplotable3; my@varpos=();
	for (my$i=0; $i<=$#haplotable2; $i++) {@{$haplotable3[$i]}=split('', $haplotable2[$i])};
	
	POSITION: for (my$j=0; $j<=$maxseqlength; $j++) {
		for (my$i=0; $i<$#haplotable3; $i++) { if ($haplotable3[$i][$j] ne $haplotable3[$i+1][$j]) 
	{push(@varpos, $j+1); next POSITION}}};
	
	my@varposdigits; for (my$i=0; $i<=$#varpos; $i++) {@{$varposdigits[$i]}=split('',$varpos[$i])};
	for (my$i=0; $i<=$#varpos; $i++) {for (my $j=1+$#{$varposdigits[$i]}; $j<=$#{$varposdigits[$#varpos]}; 
	$j++) {$varposdigits[$i][$j]=' '}};
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][0]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][1]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][2]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][3]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][4]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][5]}; print F "\r\n";
	
	my @haplonumbers; for (my$i=0; $i<$number_of_haplotypes; $i++) {$haplonumbers[$i]=1+$i};
	my @haplonumbersdigits;
	for (my$i=0; $i<=$#haplonumbers; $i++) {@{$haplonumbersdigits[$i]}=split('',$haplonumbers[$i])};
	for (my$i=0; $i<=$#haplonumbers; $i++) {for (my $j=1+$#{$haplonumbersdigits[$i]}; $j<=6; $j++) 
	{$haplonumbersdigits[$i][$j]=' '}};
	
	my @samplesizes; for (my$i=0; $i<=$#haplonumbers; $i++) {$samplesizes[$i]=$#{$haplotable[$i]}+1};
	my @samplesizesdigits;
	for (my$i=0; $i<=$#samplesizes; $i++) {@{$samplesizesdigits[$i]}=split('',$samplesizes[$i])};
	for (my$i=0; $i<=$#samplesizes; $i++) {if ($samplesizes[$i]>999) {print "One haplotype is present more than 999 times. As a result, no Network input file can be generated."; close F; exit} elsif ($samplesizes[$i]<10) {unshift(@{$samplesizesdigits[$i]},'  ')} elsif ($samplesizes[$i]<100) {unshift(@{$samplesizesdigits[$i]},' ')}};
	
	for (my$i=0; $i<=$#haplonumbers; $i++) {print F @{$haplonumbersdigits[$i]}; 
		for (my $j=0; $j<=$#varpos; $j++) {print F "$haplotable3[$i][$varpos[$j]-1]"};
		print F @{$samplesizesdigits[$i]}; print F "\r\n"};
	print F "\r\n";
	my $nb=$#varpos; $nb++;
	while ($nb>125) { print F "1010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010\r\n"; $nb=$nb-125}; 
	my $string='10' x $nb; print F $string;
close F;	}
}


print "\n";



