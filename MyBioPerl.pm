sub info{
	print "just for study.\n";
	print "report the various statistics on DNA/Protein sequences.\n";
}

#-------------------------------------------------------------------------#
#       report the various statistics on DNA sequences                    # 
#-------------------------------------------------------------------------#

#######################################################################
# nucleotides.
#######################################################################
my @nucleotides = ('A','T','C','G'); 

#######################################################################
# iub2character.
#######################################################################
my %iub2character = (
	A => 'A',
	C => 'C',
	G => 'G',
	T => 'T',
	R => '[GA]',
	Y => '[CT]',
	M => '[AC]',
	K => '[GT]',
	S => '[GC]',
	W => '[AT]',
	B => '[CGT]',
	D => '[AGT]',
	H => '[ACT]',
	V => '[ACG]',
	N => '[ACGT]',
);

#######################################################################
# IUB to regexp.
#######################################################################
sub iub2regexp{
	my($iub) = @_;
	my $regular_expression = '';

	#remove the ^ signs from the recognition sites
	$iub =~ s/\^//g;

	for(my $i = 0; $i < length($iub); ++$i){
		$regular_expression.= $iub2character{substr($iub,$i,1)};
	}
	return $regular_expression;
}

#######################################################################
# IUB to complement.
#######################################################################
sub complementIUB{
	my($seq) = @_;
	(my $com = $seq) =~ tr[ACGTRYMKSWBDHVNacgtrymkswbdhvn]
							[TGCAYRKMWSVHDBNtgcayrkmwsvhdbn];
	return $com;
}

#######################################################################
# IUB to reverse complement.
#######################################################################
sub revcomIUB{
	my($seq) = @_;
	my $revcom = reverse complementIUB($seq);
	return $revcom;
}

#######################################################################
# Parse REBASE bionet file
# input : bionet file
# output: a hash
#	  key = restriction enzyme name
#	  value = whitespace-seperated recognition site and regular expression
#######################################################################
sub parseREBASE{
	my($rebasefile) = @_;

	my $name;
	my $restriction_site;
	my $regexp;

	my %rebase_hash = ();

	my $rebase_file_handle =open_file($rebasefile);

	while(<$rebase_file_handle>){
		# discard header lines
		# the range operator(..),(1 .. /Rich Roberts/) are interesting.
		# strange:
		# if openfile to @lines then foreach @lines,it doesn't work
		( 1 .. /Rich Roberts/ ) and next;

		# discard blank lines
		/^\s*$/ and next;

		# split the two(or three if includes parenthesized name) fields
		my @fields = split(" ", $_);

		# get and store the name and the recognition site
		# remove parenthesized names , for simplicity's sake
		# just the first and the last
		$name = shift @fields;
		$restriction_site = pop @fields;

		# translate the recognition sites to regular expression
		$regexp = iub2regexp($restriction_site);

		# store the data into the hash
		$rebase_hash{$name} = "$restriction_site $regexp";
	}

	# return the hash containing the reformatted REBASE data
	return %rebase_hash;
}

#######################################################################
# Initialize annotation.
# make a blank string of same length as the given sequence string.
#######################################################################
sub initialize_annotation{
	my($seq) = @_;
	return ' ' x length($seq);
}

#######################################################################
# add annotations
# add annotation to an annotation string array.
#######################################################################
sub add_annotation{
    # $array is a reference to an array of annotations.
    my($array, $enz, @positions) = @_;
    
    # put the labels for the enzyme name at the correct positions
    # in the annotation .
    foreach my $location (@positions){
	# loop through all the annotations
	for( my $i = 0; $i < @$array; ++$i){
	    if( substr($$array[$i], $location-1, length($enz))
	    	eq (' ' x length($enz)) ){
	        substr($$array[$i], $location-1, length($enz)) = $enz;
	        last;

	    # if the annotation does collide,
	    # add it to the next annotation string
	    # on the next iteration of the 'for' loop.
	    #
	    # if there is not another annotation string, make one.
	    }elsif($i == (@$array -1)){
	    	push(@$array,initialize_annotation($$array[0]));	
	    }
	}
    }
}

#######################################################################
# format restriction map
# interleave a sequence with annotation lines showing the locations
# of restriction enzymes.
#######################################################################
sub format_restriction_map{
	my($line_length,$seq,@annotations) = @_;
	my @output = ();

	# split strings into lines of $line_length
	for( my $pos = 0; $pos < length($seq); $pos += $line_length ){
		foreach my $string ( reverse ($seq, @annotations) ){
			# discard blank lines
			if( substr($string, $pos, $line_length) 
				ne (' ' x $line_length)){

			my $length =length(substr($string,$pos,$line_length));
				
				if( substr($string, $pos, $line_length)
					ne ' 'x $length){ 
					push(@output, 
					substr($string, $pos, $line_length)
					."\n");
				}
			}
		}
		# separate the lines
		push(@output,"\n");	# yes, it is pretty good if separated.
	 }
	return @output;
}

#######################################################################
# Make restriction map from user queries on names of 
# restriction enzymes.
#######################################################################
sub make_restriction_map_from_user_queries{

	my %rebase_hash = ();
	my @dna_file_data = ();

	my $query = '';
	my $dna = '';
	my $site = '';		# recognition_site
	my $regexp = '';
	my @location = ();

	# read in the DNA file "sample.dna"
	@dna_file_data = get_file_data("sample.dna");
	
	# extract the DNA sequence data from the contents of file "sample.dan"
	$dna = extract_sequence_from_fasta_data(@dna_file_data);
	
	%rebase_hash = parseREBASE('bionet');
	
	# prompt user for restriction enzyme names , create restriction map.
	do{
		print "search for what restriction site or sites (or quit )?:";
		$query = <>;
		chomp $query;
	
		# there are two quit : one is here ,
		# another one in the until().
		# infact ,the until()'s quit process won't quit --
		#	it will quit here.
		if( ($query =~ /^\s*$/) or ($query =~ /^\s*quit/i) ){
			exit;
		}

		# Make an array of restriction site names.
		my @queries = split(' ', $query);
		
		# Initialize the annotation array.
		my @annotations = ();
		push(@annotations, initialize_annotation($dna));
	
		foreach my $query (@queries){

			# perform the search in the DNA sequence

			# change the if(exists..) to if(not exists..)
			# if requested restriction enzyme $query is not known.
			if(not exists $rebase_hash{$query}){
				print "the restriction enzyme $query "
					."is not known!\n";
				next;
			}

			# obtain the regular expression of the restriction site.
			($site, $regexp) = split(' ',$rebase_hash{$query});
				
			#create the restriction map
			@locations = match_positions($regexp , $dna);

			# If the restriction sites do NOT exist in the sequence
			if( not @locations){
				print "A restriction site for $query is"
					."NOT in the sequence\n";
				next;
			}
	
			# report the locations of the restriction map and 
			# the presence of the restriction site to the user.
			
			print "search for $query $site $regexp\n";
			print "A site for $query at location:\n";
			print join(" " , @locations), "\n";

			# create the annotations of the restriction sites
			add_annotation( \@annotations, $query, @locations);

			print "\n";
		}

		# Format the restriction map as sequence and annotations
		my @output = format_restriction_map(50,$dna,@annotations);

		# the last step : print it.
		print @output;

	}until($query =~ /quit/);
} # looks like it is not a good idea.

#######################################################################
# get restriction digest.
# the fragments of DNA left after performing a restriction reaction
# parse the REBASE bionet in a different manner,
# ignore restriction enzymes that are not given with a ^ 
# ( ^ indicating a cut site.)
#######################################################################
sub get_restriction_digest{
    my %rebase_hash = ();
    
    my @file_data = ();
    
    my $query ='';
    my $dna ='';
    my $recognition_site ='';
    my $regexp ='';
    
    my @locations =();
    my @digest =();

    # if there is a command-line argument , assume it is a DNA
    if(@ARGV){
	$dna = $ARGV[0];
    }else{
	#print "Input a FASTA filename:";
	#my $filename = <>;
	my $filename = "sample.dna";
	chomp $filename;

	die "$filename doesn't exist!\n" unless(-e $filename);
	@file_data = get_file_data($filename);
	$dna = extract_sequence_from_fasta_data(@file_data);
    }

    # Get the REBASE data into a hash, from file 'bionet.110'
    %rebase_hash = parseREBASE('bionet.110');

    #prompt user for restriction enzyme names, create restriction map.
    do{
	print "search for what restriction site (or quit)?:";
	$query = <>;
	chomp $query;

	# exit if empty query
	if( $query =~ /^\s*$/ ){
	    exit;
	}

	#perform the search in the DNA sequence
	if( exists $rebase_hash{$query} ){
	    ($recognition_site, $regexp) = split(" ",$rebase_hash{$query});

	    # find the offset of the cut site in the recognition site
	    my $cut_site_offset = index($recognition_site,'^');

	    # create the restriction map (?)
	    @locations = match_positions($regexp, $dna);

	    # add the cut site offset to the @locations
	    @locations = grep($_ += $cut_site_offset, @locations);

	    # calculate and report the restriction digest to the user.
	    if(@locations){
		print "Search for $query $recognition_site $regexp\n";
		print "A restriction digest for $query cutting at locations:\n";
		print join(" ",@locations),"\n";

		# extract the DNA fragments between the cut sites:
		# the restriction digest
		for( my $i = 1, my $j = shift(@locations);@locations;$i=$j,$j=shift(@locations)){
		    push(@digest,substr($dna, $i-1, $j-$i));
		}

		print "the resulting restriction digest:\n";
		print join("\n",@digest),"\n";


	    }else{
		print "A restriction enzyme $query is not in the DNA:\n";
	    }
	}
	print "\n";
    }until($query =~ /quit/);
}

#######################################################################
# parseREBASE2 , parse REBASE bionet file, version 2
# A subroutine to return a hash where
# 	key = restriction enzyme name
# 	value = whitespace-separated recognition sites
# 	( the regular expressions will be calculated on the fly)
# Version2 handles multiple definition lines for an enzyme name
# Version2 also handles alternate enzyme names on a line
#######################################################################
sub parseREBASE2{
	my($rebasefile) = @_;

	my @rebasefile = ();
	my %rebase_hash = ();
	my $site;
	my $regexp;

	# read in the REBASE file
	my $rebase_filehandle = open_file($rebasefile);
	
	while(<$rebase_filehandle>){
		my @names = ();

		# discard header lines
		(1 .. /Rich Roberts/) and next;

		# discard blank lines
		/^\s*$/ and next;

		# split the two( or three if includes parenthesized name) fields
		my @fields = split(" ", $_);

		# get and store the recognition site
		$site = pop @fields;

		# for the purposes of this exercise, we'll ignore cut sites (^)
		# this is not something you'd want to do in general.
		$site =~ s/\^//g;

		# Get and store the name and the recognition site.
		# and alternate (parenthesized) names
		# from the middle field, if any
		foreach my $name (@fields){
			if($name =~ /\(.*\)/){
				$name =~ s/\((.*)\)/$1/;	# the s is different to m
			}
			push @names, $name;
		}

		# Store the data into the hash, avoiding duplicates
		# (ignoring ^ cut sites)
		# and ignoring reverse complements
		foreach my $name (@names){

			# Add new enzyme definition
			if(not defined $rebase_hash{$name}){
				$rebase_hash{$name} = "$site";
				next;
			}
			my(@defined_sites) = split(" ",$rebase_hash{$name});

			# Omit already defined sites
			if(grep {$site eq $_} @defined_sites){
				next;
				
			# Omit reverse complements of already defined sites
			}elsif(grep {revcomIUB($site) eq $_} @defined_sites){
				next;

			# Add the additional site
			}else{
				$rebase_hash{$name} .= " $site";
			}
		}
	}

	# Return the hash containing the reformatted REBASE data.
	return %rebase_hash;
}

#######################################################################
# get nonpalindromic recognition site
# Extend the restriction map software to take into account the opposite
# strand for nonpalindromic recognition sites.
#######################################################################
sub get_nonpalindromic_recognition_site{
	my %rebase_hash = ();
	my @file_data = ();
	my $query = '';
	my $dna = '';
	my @recognition_sites = ();
	my $recognition_site = '';
	my $regexp = '';

	# if there is a command-line argument, assume it's dna
	if(@ARGV){
		$dna = $ARGV[0];
	}else{
		#print "input a FASTA filename:";
		#my $filename = <>;
		my $filename = 'sample.dna';
		chomp $filename;

		die "$filename doesn't exist!\n" unless(-e $filename);

		@file_data = get_file_data($filename);

		$dna = extract_sequence_from_fasta_data(@file_data);	
	}

	# new version of parseREBASE
	%rebase_hash = parseREBASE2('bionet.110');
	#my $i=0;
	#foreach my $key (keys %rebase_hash){
	#	print "$key -> $rebase_hash{$key}\n";
	#	$i++;
	#}
	#print "total:$i\n";

	#prompt user for restriction enzyme names, create restriction map.
	do{
		print "search for what restriction site (or quit)?: ";
		$query = <>;
		chomp $query;

		if($query =~ /^\s*$/ ){
			exit;
		}

		# translating the recognition sites to regular expressions
		if ( exists $rebase_hash{$query} ){
	my @locations = ();
			@recognition_sites = split(" ", $rebase_hash{$query});

			foreach $recognition_site (@recognition_sites){
			
				$regexp = iub2regexp($recognition_site);

				# create the restriction map
				# store these positons with a leading + sign, like "+324"
				push @locations, map($_ = "+$_", match_positions($regexp,$dna));
				
				# for non-palindromic recognition sites,
				# search for the complement in the sequence
				if($recognition_site ne revcomIUB($recognition_site)){

					# calculate the regular expression for the complement
					(my $complement = $recognition_site)=~ tr/ACGTRYMKSWBDHVNacgtrymkswbdhvn/TGCAYRKMWSVHDBNtgcayrkmwsvhdbn/;
					
					my $regular_expression_com = iub2regexp($complement);

					# Get the matching positions for the complement
					# Store these in @positions with a leading-sign,
					# like "-324"
					push @locations,map($_="-$_",match_positions($regular_expression_com,$dna));

				}

				# Sort the locations, ignoring the leading '+' and '-' sign
				@locations = sort {
					my($A,$B); 
					($A=$a)=~s/^.//;
					($B=$b)=~s/^.//;
					$A<=>$B;
				}@locations;

				# report the restriction map to the user
				if(@locations){
					print "searching for $query $recognition_site $regexp\n";
					print "A restriction site for $query at locations:\n";
					print join(" ",@locations),"\n";
				}else{
					print "A restriction enzyme $query is not in the DNA:\n";
				}
			}
		}else{
			print "$query is not in hash!\n";
		}

		print "\n";
	}until( $query =~ /quit/ );

}

#######################################################################
# match positions
#######################################################################
sub match_positions{
	my($regexp , $sequence) = @_;

	my @positions = ();
	while( $sequence =~ /$regexp/ig){
		push(@positions, pos($sequence) - length($&) +1);
	}
	return @positions;
}

#######################################################################
# print IUB table.
#######################################################################
sub print_iub_table{
	foreach my $key (keys %iub2character){
		print "$key => $iub2character{$key}\n";
	}
}

#######################################################################
# random position.
#######################################################################
sub random_position{
	my($string) = @_;
	return int rand length $string;	#maybe is better use it direct
}

#######################################################################
# random position.
#######################################################################
sub random_codon{
	my($codon) = @_;
	my $position = random_position($codon);
	
	$new_nucleotide = $nucleotides[rand @nucleotides];
	substr($codon , $position ,1, $new_nucleotide);
	return $codon;
}

#######################################################################
# report how GC-rich some sequence is.
# (just give the precentage of G and C in the DNA)
#######################################################################
sub cg_percentage{
	$dna = shift;
	return (int (100 * ($dna =~ tr/CGcg//)/(length $dna)))."%";
}

#######################################################################
# presence or absence of poly-T sequences
# (long stretches of mostly T's at the 5' (left) end of many DNA sequences)
#######################################################################
sub dna_has_poly_t{
	$dna = shift;

	# we will need at least 10 T's out of the last 12 bases
	$head = substr($dna,0,12);

	$countT = ($head =~ tr/Tt//);
	if($countT >=10){
		return $countT;
	}else{
		return 0;
	}
}

##################################################
# DNA to RNA
##################################################
sub dna2rna{
	$dna = shift;
	($rna = $dna) =~ tr/Tt/Uu/;
	#($rna = $dna) =~ s/Tt/Uu/g;	# dosen't wokr at 5.8.8.
	return $rna;
}

##################################################
# RNA to DNA
##################################################
sub rna2dna{
	$rna = shift;
	($dna = $rna) =~ tr/Uu/Tt/;
	#($dna = $rna) =~ s/Uu/Tt/g;	# just work after 5.010 (?)
	return $dna;
}

################################################
# convert to its reverse complement.
# in other words,Calculate the reverse complement 
# of an strand of DNA.
################################################
sub reverse_complement{
	$dna = shift;
	($revcom = $dna) =~ tr/ATCGatcg/TAGCtagc/;
	return reverse $revcom;
}

################################################
# two DNA are reverse complement or not
################################################
sub is_reverse_complement{
	#check argument number
	if(scalar(@_) ne 2){
		print "Error:Give 2 strings as arguments "
			."to check for reverse complements.\n";
		exit;
	}

	#change to array
	my @a = split('' , $_[0]);
	my @b = split(// , $_[1]);	#for split ,the '' and // are the same.

	#check argument length
	if(length $_[0] ne scalar(@b) ){
		print "Error:two strings' length is not same.\n";
		exit;
	}
	
	while(@a){
		my $base_a = shift @a;
		my $base_b = pop @b;
		#change $base_a to complement reverse:
		$base_a =~ tr/ATCGatcg/TAGCtagc/;
		if($base_a ne $base_b){
			print "not equal,exit .";
			exit;
		}
	}
	print "they are reverse complement.\n";

}

############################################################### 
# check a string and return true if it is a DNA sequence.
############################################################### 
sub is_DNA{
	my($sequence) = @_;

	# if there is a non-ATCG character ,return false.
	# 'i' means case-insensitive.
	# '^' means any character other than the ones listed.
	if($sequence =~ /[^ATCG]/i){
		return 0;
	}
	return 1;	# otherwise , return true;
}

############################################################### 
# switches two bases in a DNA string at specified positions.
# there have two way to do it : use substr or splice.
# substr editon ( 1st of 3 edition , first positon is 0)
############################################################### 
sub switch_two_bases_by_substr{
	my($dna,$position1,$position2) = @_;
	print $dna . " $position1" . " $position2"."\n";

	#the bases at two position
	$base1 = substr($dna,$position1,1);
	$base2 = substr($dna,$position2,1);

	#exchange
	substr($dna,$position1,1) = $base2;
	substr($dna,$position2,1) = $base1;

	return $dna;

}

############################################################### 
# switches two bases in a DNA string at specified positions.
# splice edition ( 2nd of 3 edition, first position is 0.)
############################################################### 
sub switch_two_bases_by_splice{
	my($str , $position1 , $position2) = @_;
	print $str . " $position1" . " $position2"."\n";

	# for use splice , must first explode str into array.
	@str = split('',$str);

	# get base by array
	$base1 = $str[$position1];
	$base2 = $str[$position2];

	# splice them.
	splice(@str , $position1 , 1 , $base2);
	splice(@str , $position2 , 1 , $base1);

	#print it for test.
	#print @str;
	#print "\n";
	# if use print @str."\n";   will get @str 's length .
	# if use Print "@str";		will get @str 's string seperated by ' '.
	
	#change array to string and return .
	$str = join('',@str);
	#print $str;

	return $str;
}

############################################################### 
# switches two bases in a DNA string at specified positions.
# array edition ( 3rd of 3 edition, first position is 0.)
############################################################### 
sub switch_two_bases_by_array{
	my($str, $position1, $position2) = @_;
	print $str . " $position1 ". $position2. "\n";
	@str = split('',$str);
	@str[$position1, $position2] = @str[$position2,$position1];
	print @str;
}

############################################################### 
# String operation , to lower case
############################################################### 
sub to_lower_string{
	$str = shift;
	$str = "\L$str";
	return $str;
}

############################################################### 
# String operation , to upper case
############################################################### 
sub to_upper_string{
	$str = shift;
	$str = "\U$str";	#or use 'uc'
	return $str;
}

############################################################### 
# the frequency of nucleotides
# report the percentage of each nucleotide .
############################################################### 
sub frequency_of_nucleotides{
	my($filename) = @_;
	open(FILE,$filename) or die "Can't open file:$filename\n";
	my @file = <FILE>;
	close FILE;
	
	#print @file;
	#print "\n";
	#print scalar @file."\n";	# just counted the lines.

	$dna = join('',@file);
	$dna =~ s/\s//g;		# remove space.
	#print $dna."\n";
	$length = length $dna;
	print "total length:".$length."\n";

	# check for the length. it is illegal to have a length of 0.
	if($length == 0){
		print "length is 0 . exit.\n";
		return 0;
	}

	# now we are initialize the counts .
	$count_of_a = 0;
	$count_of_t = 0;
	$count_of_c = 0;
	$count_of_g = 0;
	$count_of_errors = 0;

	# and then , there is 5 way to do it.
	# 1. 
	$count_of_a = nucleotide_percentage($dna,'a');
	$count_of_t = nucleotide_percentage($dna,'t');
	$count_of_c = nucleotide_percentage($dna,"c");
	$count_of_g = nucleotide_percentage($dna,"g");
}

############################################################### 
# count nucleotide
# and return the percentage in the DNA
############################################################### 
sub nucleotide_percentage{
	my($dna, $nucleotide) =@_;
	
	my $count = 0;

	# to interpolate the variable $nucleotide in the tr/// function
	# we have to use an 'eval' on the pattern match.
	#
	# and to avoid eval's string being 'strict' by 'use strict',
	# we can use 'no stirct',(but can NOT do it in module).
	# {
	# no strict;
	# $count2 = eval "$dna =~ tr/$nucleotide//";	# 'c'+'t' = 477
	# 
	# #$count = ($dna =~ tr/$nucleotide//);		#same 4 value:477.
	# print "use tr dirctly2:$count2.\n";
	# }
	#
	# so use while and m// instead.

	while($dna =~ /$nucleotide/ig){
		$count++;
	}
	print "count of $nucleotide:".$count.",";

	$percentage =( $count/(length $dna) *100);
	print "percentage :$percentage %\n";

	return $count;
}


############################################################### 
# search by name for a gene in an unsorted array.
# this is same as search_gene_name ,but using foreach. 
############################################################### 
sub search_gene_name_by_foreach{
	my @genes = qw( xgene1 ygene2 zgene3 agene4 bgene5);
	
	print "enter gene name:";
	my $gene_name = <>;
	chomp $gene_name;

	my $found_flag = 0;

	foreach(@genes){
		print $_,"\n";

		#look for the gene_name.
		if(/^$gene_name$/){
			print "gene $gene_name is known.\n";
			$found_flag = 1;

			#exit the 'foreach' loop.
			last;
		}
	}

	unless($found_flag){
		print "gege $gene_name is not known.\n";
	}

}

############################################################### 
# search name by binary search
############################################################### 
sub search_gene_name_by_binarysearch{
	my($genes,$gene_name) = @_;	#genes is array: @$genes
	my @sorted_genes = sort @$genes;

	if(binary_search_recursive($gene_name,@sorted_genes)){
		print "Gene $gene_name is known.\n";
	}else{
		print "Gene $gene_name is NOT known.\n";
	}
}


############################################################### 
# binary search
############################################################### 
sub binary_search{
	my($query,@array) = @_;

	#if the array is empty , failure. (do I need print sth?)
	@array or return 0;

	my $begin = 0;

	# $#array is the offset of the last element of @array.
	# the same number could be obtained by scalar(@array)-1
	my $end = $#array;

	while(($end - $begin) >= 0){

		my $mid_point = int (($end -$begin)/2 + $begin);
		print $array[$mid_point] , ":" , "@array[$begin..$end]" , "\n";

		# " a cmp b" return 0 if operands are equal
		# -1 if a < b   ( like 'a' cmp 'b');
		# 1  if a > b   ( like 'z' cmp 'b');
		my $result = $query cmp $array[$mid_point];

		if($result == 0){
			return $query;
		}elsif($result == -1){
			$end = $mid_point - 1;
		}else{	# reuslt == 1
			$begin = $mid_point + 1;
		}
	}
	return 0;

}

############################################################### 
# binary search recursive
############################################################### 
sub binary_search_recursive{
	my($query,@array) = @_;

	# if the array is empty ,failure (now we will use it)
	@array or return 0;

	my $mid_point = int (@array/2);
	print "Mid point $mid_point with value:",$array[$mid_point]."\n";

	my $result = $query cmp $array[$mid_point];

	if($result == 0){
		return $query;
	}elsif($result == -1){
		return binary_search_recursive($query,@array[0..$mid_point -1]);
	}else{
		return binary_search_recursive($query,@array[$mid_point+1 .. $#array]);
	}
}

############################################################### 
# search by name for a gene in an unsorted array.
# this version use grep
# version 1 : argument , scalar first , array second.
# (easy things should be easy)
############################################################### 
sub search_gene_name{
	#my @genes = qw(xgene1  ygene2  zgene3  agene4  bgene5);

	#print "enter gene name:";
	#my $gene_name = <>;
	#chomp $gene_name;

	my($gene_name,@genes) = @_;	# ok

	#my($genes,$gene_name) = @_;
	
	print "genes:@genes    \n";			# use @genes direct
	print "gene_name:$gene_name\n";

	#/^ $gene_name $/ , beginning and ending .
	if(grep(/^$gene_name$/,@genes)){		# use @genes direct too.
		print "gene $gene_name is known.\n";
	}else{
		print "gene $gene_name is not known.\n";
	}
}

############################################################### 
# search by name for a gene in an unsorted array.
# this version use grep
# version 2: arguments , array first , scalar second.
# the caller must use '\' before array.
# example : search_gene_name_2(\@genes,gene_name);
############################################################### 
sub search_gene_name_2{
	#my @genes = qw(xgene1  ygene2  zgene3  agene4  bgene5);

	#print "enter gene name:";
	#my $gene_name = <>;
	#chomp $gene_name;

	#my($gene_name,@genes) = @_;	# ok

	my($genes,$gene_name) = @_;
	
	print "genes:@$genes    \n";	# the array bypassed are @$genes.
	print "gene_name:$gene_name\n";

	#/^ $gene_name $/ , beginning and ending .
	if(grep(/^$gene_name$/,@$genes)){			#@$genes too.
		print "gene $gene_name is known!\n";
	}else{
		print "gene $gene_name is not known!\n";
	}
}

############################################################### 
# inserts an element into a sorted array.
############################################################### 
sub insert_an_element_into_sorted_array{
	my($array,$element) = @_;
	my @sorted_array = sort @$array;

	my $i;
	for($i = 0 ; 
	  (@sorted_array > $i) and (($element cmp $sorted_array[$i]) != -1);
	  ++$i){
		;
	}
	splice(@sorted_array , $i , 0 , $element);
	print "new sorted array :@sorted_array\n";
	return @sorted_array;
}

#---------------------------------------------------------------------------#
# report the various statistics on Protein sequences                        #	
#---------------------------------------------------------------------------#

############################################################### 
# check a string and return true if it is protein sequence. 
############################################################### 
sub is_protein{
	my($sequence) = @_;

	# if there's a non-amino acid character ,return false;
	if($sequence =~ /[^ACDEFGHIKLMNPQRSTVWXYZ]/i){
		return 0;
	}

	#otherwise , return 'true'
	return 1;
}

############################################################################ 
# report on the percentage of hydrophobic amino acids in a protein sequence.
# hydrophobic amino acids are the amino acids with nonpolar side chains.
# also called 'hydrophobic residues'.
############################################################################    
sub count_hydrophobic{
	$protein = shift;
	$hydrophobic = ($protein =~ tr/GAVLIPFMWC//);
	return $hydrophobic;
}

############################################################### 
# return the Percentage of hydrophobic (with '%')
############################################################### 
sub percentage_hydrophobic{
	$protein = shift;
	$hydrophobic = ($protein =~ tr/GAVLIPFMWC//);
	$percentage = 100 * ($hydrophobic / length $protein)."%";
	return $percentage;
}

############################################################### 
# Get Amino Acid
# Three Letter AA to one Letter AA
############################################################### 
sub get_amino_acid_table{
	%three_to_one = (
		ALA => A, CYS => C, ASP => D, GLU => E,
		PHE => F, GLY => G, HIS => H, ILE => I,
		LYS => K, LEU => L, MET => M, ASN => N,
		PRO => P, GLN => Q, ARG => R, SER => S,
		THR => T, VAL => V, TRP => W, TYR => Y
	);
	return %three_to_one;
}

############################################################### 
# display Amino Acid
############################################################### 
sub display_amino_acid{
	my %three_to_one = get_amino_acid_table();

	$i=0;
	while ( ($key, $value) = each %three_to_one){
		$i++;
		print "$key => $value\t ";
		if( $i % 4 == 0){
			print "\n";
		}
	}

}

############################################################### 
# display Amino Acid
############################################################### 
sub display_amino_acid_sort_by_3{
	my %three_to_one = get_amino_acid_table();

	#the default value of $i is 0 ,so we can ignore its define here.

	foreach $key (sort keys %three_to_one){
		$i++;
		print "$key => $three_to_one{$key} ,";
		if($i%4 == 0){
			print "\n";
		}
	}
}


############################################################### 
# display Amino Acid
############################################################### 
sub display_amino_acid_sort_by_1{
	my %three_to_one = get_amino_acid_table();
	%one_to_three = reverse %three_to_one;

	#the default value of $i is 0 ,so we can ignore its define here.

	foreach $key (sort keys %one_to_three){
		$i++;
		print "$key => $one_to_three{$key} ,";
		if($i%4 == 0){
			print "\n";
		}
	}
}

############################################################### 
# Extract annotation and sequence from GeneBank file
############################################################### 
sub extract_annotation_and_sequence_from_genebank_file{
	my @annotation =();
	my @sequence ='';
	my $filename = 'record.gb';

	parse1(\@annotation, \$sequence, $filename);

	print "\nannotation:\n";
	print @annotation;

	print "\nsequence:\n";
	print_sequence($sequence, 50);
}

############################################################### 
# open file 
############################################################### 
sub parse1{
	my($annotation, $dna, $filename) = @_;

	# $annotation 	--reference to array
	# $dna 		--reference to scalar
	# $filename 	--scalar

	# declare and initialize variables
	my $in_sequence = 0;
	my @genebank_file =();

	# Get the genebank data into an array from a file
	@genebank_file = get_file_data($filename);

	# Extract all the sequence lines
	foreach my $line (@genebank_file){

		if($line =~ /^\/\/\n/){
			last;
		}elsif($in_sequence){
			$$dna .= $line;
		}elsif($line =~ /^ORIGIN/){
			$in_sequence = 1;
		}else{
			push( @$annotation, $line);
		}
	}
	$$dna =~ s/[\s0-9]//g;
}

############################################################### 
# open file 
############################################################### 
sub open_file{
	my($filename) = @_;
	my $fh;

	open($fh,$filename) or die("can't open $filename\n");
	return $fh;
}

############################################################### 
# get file data
############################################################### 
sub get_file_data{
	my($filename) = @_;
	my @filedata = ();
	die("Can't open:$filename\n") unless( open(GET_FILE_DATA, $filename) );
	@filedata = <GET_FILE_DATA>;
	close GET_FILE_DATA;
	return @filedata;
}

############################################################### 
# extract sequence from fasta data
# a subroutine to extract FASTA sequence data from an array
############################################################### 
sub extract_sequence_from_fasta_data{
	my(@fasta_file_data) = @_;
	my $sequence ='';
	foreach my $line (@fasta_file_data){
		if($line =~ /^\s*$/){
			next;
		}elsif($line =~ /^\s*#/){
			next;
		}elsif($line =~ /^>/){
			next;
		}else{
			$sequence.=$line;
		}
	}
	# remove non-sequence data (whitespace) from sequence data
	$sequence =~ s/\s//g;
	return $sequence;
}

sub print_sequence{
	my($sequence, $length) = @_;
	for( my $pos = 0; $pos < length($sequence); $pos +=$length){
		print substr($sequence, $pos, $length),"\n";
	}
}

############################################################### 
# checks an array of data if it's in 
############################################################### 
sub is_fasta{
	my(@contents) = @_;
	my $line = shift @contents;
	unless($line =~ /^>/ ){
		print "error: not begin with > \n";
		return 0;
	}

	# dose is work ?
	#@contents = uc @contents;

	#check if DNA (standard IUB/IUPAC nucleic acid codes
	#if(grep(/[^ACGTUMRWSYKVHDBNacgtumrwsykvhdbn]/,@contents)){
	if(grep(/[^acgt]/,@contents)){
		print "error: not dna \n";
		return 0;
	}

	#check if Protein  
	#            B ?!
#	if(grep(/[^-ABCDEFGHIKLMNPQRSTVWXYZ\*]/,@contents)){
#		return 0;
#	}

	return 1;
}

############################################################### 
# get genetic codons
#
############################################################### 
sub get_genetic_codons{

	my(%genetic_code) = (
		'TCA' => 'S',	#Serine
		'TCC' => 'S',	#Serine
		'TCG' => 'S',	#Serine
		'TCT' => 'S',	#Serine

		'TTC' => 'F',	#Phenylalanine
		'TTT' => 'F',	#Phenylalanine

		'TTA' => 'L',	#Leucine
		'TTG' => 'L',	#Leucine

		'TAC' => 'Y',	#Tyrosine
		'TAT' => 'Y',	#Tyrosine
		
		'TAA' => '_',	#Stop
		'TAG' => '_',	#Stop

		'TGC' => 'C',	#Cysteine
		'TGT' => 'C',	#Cysteine

		'TGA' => '_',	#Stop

		'TGG' => 'W',	#Tryptophan

		'CTA' => 'L',	#Leucine
		'CTC' => 'L',	#Leucine
		'CTG' => 'L',	#Leucine
		'CTT' => 'L',	#Leucine

		'CCA' => 'P',	#Proline
		'CCC' => 'P',	#Proline
		'CCG' => 'P',	#Proline
		'CCT' => 'P',	#Proline

		'CAC' => 'H',	#Histidine
		'CAT' => 'H',	#Histidine

		'CAA' => 'Q',	#Glutamine
		'CAG' => 'Q',	#Glutamine

		'CGA' => 'R',	#Arginine
		'CGC' => 'R',	#Arginine
		'CGG' => 'R',	#Arginine
		'CGT' => 'R',	#Arginine

		'ATA' => 'I',	#Isoleucine
		'ATC' => 'I',	#Isoleucine
		'ATT' => 'I',	#Isoleucine

		'ATG' => 'M',	#Methionine

		'ACA' => 'T',	#Threonine
		'ACC' => 'T',	#Threonine
		'ACG' => 'T',	#Threonine
		'ACT' => 'T',	#Threonine

		'AAC' => 'N',	#Asparagine
		'AAT' => 'N',	#Asparagine

		'AAA' => 'K',	#Lysine
		'AAG' => 'K',	#Lysine

		'AGC' => 'S',	#Serine
		'AGT' => 'S',	#Serine

		'AGA' => 'R',	#Arginine
		'AGG' => 'R',	#Arginine

		'GTA' => 'V',	#Valine
		'GTC' => 'V',	#Valine
		'GTG' => 'V',	#Valine
		'GTT' => 'V',	#Valine
		
		'GCA' => 'A',	#Alanine
		'GCC' => 'A',	#Alanine
		'GCG' => 'A',	#Alanine
		'GCT' => 'A',	#Alanine
		
		'GAC' => 'D',	#Aspartic Acid
		'GAT' => 'D',	#Aspartic Acid

		'GAA' => 'E',	#Glutamic Acid
		'GAG' => 'E',	#Glutamic Acid

		'GGA' => 'G',	#Glycine
		'GGC' => 'G',	#Glycine
		'GGG' => 'G',	#Glycine
		'GGT' => 'G',	#Glycine
	);

	return %genetic_code;
}

############################################################### 
# codon2aa
# tanslate a DNA 3-character codon to an amino acid
# using hash lookup
############################################################### 
sub codon2aa{
	my($codon) = @_;
	$codon = uc $codon;

	my(%genetic_code) = get_genetic_codons();

	if(exists $genetic_code{$codon}){
		return $genetic_code{$codon};
	}else{
		print STDERR "Bad codon \"$codon\"!!\n";
	}
}

############################################################### 
# dna to peptide
############################################################### 
sub dna2peptide{
	my($dna) = @_;
	my $peptide = '';
	for(my $i = 0; $i < (length($dna) - 2); $i += 3){
		$peptide.= codon2aa( substr($dna, $i ,3) );
	}
	return $peptide;
}

############################################################### 
# translate frame
############################################################### 
sub translate_frame{
	my($seq, $start, $end) = @_;
	my $protein;

	unless($end){
		$end = length($seq);
	}

	# calculate and return the translation.
	return dna2peptide( substr($seq, $start -1, $end - $start +1) );
}

############################################################### 
# mutate codons
# given a codon, returns a list of all the amino acids
# that can result from any single mutation in the codon.
# 9 in all : 3 x (4-1) = 9
# input : codon
# output: mutated codon
############################################################### 
sub mutate_codons{
	my($codon) = @_;
	my @bases = ('A','C','G','T');

	# to collect the output mutations
	my @mutations = ();

	# there have two loop:
	# outer loop loops through the 3 positions in the codon,
	# inner loop loops over the possible bases.
	for(my $i = 0; $i < 3; ++$i){
		foreach my $base (@bases){
			# the codon are 'static'
			my $mutation = $codon;

			if($base eq substr($codon, $i ,1)){
				next; #next foreach loop
			}else{
				substr($mutation, $i, 1) = $base;
			}
			push(@mutations, $mutation);
		}
	}
	return @mutations;
}

############################################################### 
# codon mutated and translated
# tanslate the mutated codon .
# the hash statement has the effect of defining each amino acid
# only once , so the return statement for the subroutine doesn't
# give an amion acid twice if there are two mutations that 
# result in that amino acid.
############################################################### 
sub codon_mutated_and_translated{
	my($codon) = @_;
	%translated_codons = ();
	my @mutated_codons = mutate_codons($codon);

	print "the mutated codon of $codon are:@mutated_codons\n";

	foreach my $mutated_codon(@mutated_codons){
		#print "XX $mutated_codon -> ",codon2aa($mutated_codon),"\n";

		$translated_codons{codon2aa($mutated_codon)}++;
		# ^  remember the foolish mistake
	}
	@return_array = keys %translated_codons;
	return keys %translated_codons;
}

sub get_reversed_genetic_codon{
	
	my %reversed_genetic_code = (
		'A' => 'GCA GCC GCG GCT',	# Alanine
		'C' => 'TGC TGT',		# Cysteine
		'D' => 'GAC GAT',		# Aspartic Acid
		'E' => 'GAA GAG',		# Glutamic Acid
		'F' => 'TTC TTT',		# Phenylalamine
		'G' => 'GGA GGC GGG GGT',	# Glycine
		'H' => 'CAC CAT',		# Histidine
		'I' => 'ATA ATC ATT',		# Isoleucine
		'K' => 'AAA AAG',		# Lysine
		'L' => 'CTA CTC CTG CTT TTA TTG',#Leucine
		'M' => 'ATG',			# Methionine
		'N' => 'AAC AAT',		# Asparagine
		'P' => 'CCA CCC CCG CCT',	# Proline
		'Q' => 'CAA CAG',		# Glutamine
		'R' => 'CGA CGC CGG CGT AGA AGG',#Arginine
		'S' => 'TCA TCC TCG TCT AGC AGT',#Serine
		'T' => 'ACA ACC ACG ACT',	# Threonine
		'V' => 'GTA GTC GTG GTT',	# Valine
		'W' => 'TGG',			# Tryptophan
		'Y' => 'TAC TAT',		# Tyrosine
		'_' => 'TAA TAG TGA',		# stop 
	);                                                       
	return %reversed_genetic_code;
}

############################################################### 
# Given one-character code for an amino acid, 
# randomly changes it to one of the amino acids.
# input : amino acid  ( one-character )
# output: amino acid
############################################################### 
sub random_mutate_aa{
	my($aa) = @_;
	$aa = uc $aa;

	#check peptide 
	die "$aa is not a protein sequence \n" unless is_protein($aa);

	my %reverse_genetic_code = get_reversed_genetic_codon();

	# calculate the underlying codon.
	print ">reverse_genetic_code:",(keys %reverse_genetic_code),"\n";
	my @possible_codons = split(' ', $reverse_genetic_code{$aa} );

	#dbg
	print "possible codons : @possible_codons\n";

	# get the list of possible amino acids
	my @possible_aa = ();

	foreach my $codon (@possible_codons){
		push(@possible_aa,codon_mutated_and_translated($codon));
	}

	print "**possible_aa : @possible_aa \n";

	# the @possible_aa dosen't have chance include some amino acids 
	# more than once because we use hash here
	%possible_amino_acid = ();
	foreach my $amino_acid (@possible_aa){
		$possible_amino_acid{($amino_acid)}++;
	}
	@possible_amino_acid = keys %possible_amino_acid;
	print "unique aa:@possible_amino_acid\n";
	print "possible codons:@possible_codons\n";

	# Finnally, randomly select one and return it.
	# warning : use @possible_codons is NOT right
	#		it is better get @possible_amino_acid and use it here.
	#return $possible_aa[rand @possible_codons];
	return $possible_aa[rand @possible_aa];
	print "\n";
}

############################################################### 
# return aa 's probability
############################################################### 
sub aa_probability{
	my($aa) = @_;
	die "$aa is not a protein sequence \n" unless is_protein($aa);

	my(%reverse_genetic_code) = get_reversed_genetic_codon();

	foreach my $key (keys %reverse_genetic_code){
		print "$key -> $reverse_genetic_code{$key}\n";
	}

	$aa = uc $aa;

	# use of implicit split to array is deprecated
	#my $number_of_codons = split(' ', $reverse_genetic_code{$aa});
	
	my @codons = split(' ', $reverse_genetic_code{$aa});
	my $number_of_codons = @codons;
	return $number_of_codons/64;
}


############################################################### 
# takes each codon that encodes the specified amino acid
# and mutates it at the specified position to the specified nucleotide.
# input :an amino acid , position , nucleotide.
# output:a set of amino acids that are encoded by the mutated codons
############################################################### 
sub aa2mutated_aa{
	my($aa , $position, $nucleotide) =@_;
	$aa = uc $aa;
	my(%reverse_genetic_code) = get_reversed_genetic_codon();
	my %mutated_aa = ();
	my @codons = split(' ',$reverse_genetic_code{$aa});
	print "@codons.\n";

	foreach my $codon (@codons){
		print ">$codon";
		substr($codon, $position-1, 1) = $nucleotide;
		print ">>$codon";
		my $new_aa = codon2aa($codon);
		print ">>>$new_aa\n";
		$mutated_aa{$new_aa}++;
	}
	return keys %mutated_aa;
}


############################################################### 
# mutated probability
############################################################### 
sub mutation_probability{
	my($aa1,$aa2) = @_;
	my $count = 0;
	my(%reverse_genetic_code) = get_reversed_genetic_codon();
	
	# the set of codons that code the first aa
	my @codons = split(' ', $reverse_genetic_code{$aa1});
	foreach my $codon (@codons){
		foreach my $mutation (mutate_codons($codon)){
			if(codon2aa($mutation) eq $aa2){
				$count++;
			}
		}
	}

	return $count/(9 * @codons);
}


1;
