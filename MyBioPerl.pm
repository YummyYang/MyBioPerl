sub info{
	print "just for fun.\n";
	print "play with various statistics on DNA/Protein sequences.\n";
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

	    print ">recognition_site:$recognition_site," 
	    	   ."regexp:$regexp,cut_site_offset:$cut_site_offset\n";


	    # create the restriction map (?)
	    @locations = match_positions($regexp, $dna);
	    print ">>locations:@locations\n";

	    # add the cut site offset to the @locations
	    @locations = grep($_ += $cut_site_offset, @locations);
	    print ">>>locations:@locations\n";

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

############################################################### 
# read database
# read in the data and then do a fast lookup
# on the information associated with a gene name.
############################################################### 
sub read_database{
	my($dbfile) =@_;

	# open the database file
	die "can't open $dbfile!\n" unless( open(DATAFILE,"<$dbfile") );

	my %database = ();

	# Set input separator	
	# This allows us to get an entire record as one scalar value
	# A blank line is just a newline following another newline.
	$/ = "\n\n";

	# Read data
	while(my $record = <DATAFILE>){
		# Get the first line as a key, and the remaining lines as the value
		# Trun the scalar $record back into an array of lines to get the key
		my @lines = split(/\n/,$record);
		my $key = shift @lines;

		# And trun the remaining array into a scalar for storing as a value
		# in the hash
		my $value = join("\n", @lines) . "\n";
		$database{$key} = $value;
	}

	close(DATAFILE);

	# Reset input separator to normal default value
	$/ = "\n";

	return %database;
}

############################################################### 
# input new record
############################################################### 
sub input_new_record{
	print "Input a new record (end with a blank line):\n";
	$/ = "\n\n";

	my $new_record = <>;
	$/ = "\n";

	return join("\n", split("\n", $new_record) ),"\n";
}

############################################################### 
# input new record
############################################################### 
sub add_record{
	my($dbfile, @record) = @_;

	die "Can't write to $dbfile!\n" unless(open(DATAFILE,">>$dbfile"));

	print DATAFILE "\n", @record;
	close DATAFILE;

}

############################################################### 
# db demo
############################################################### 
sub db_demo{
	# the ASCII flat file that holds the database
	my $dbfile = 'db.demo';

	# Read in the db from the flat file to a hash
	my %data = read_database($dbfile);

	# Test by examing a key-value pair
	print "See what entry?: ";
	my $ans = <>;
	chomp $ans;
	print $data{$ans};

	# Get a new record for the database from the user
	my @new_record = input_new_record();

	# Add a new record to db
	add_record($dbfile,@new_record);

	# Reload the database from the database file
	%data = read_database($dbfile);

	# Test by examining a key-value pair
	print "See what entry?: ";
	$ans = <>;
	chomp $ans;
	print $data{$ans};
}

############################################################### 
# Print Word frequence
# takes a long DNA sequence as input 
# output the counts of all four-base subsequences,(256 in all)
# sorted by frequency.
############################################################### 
sub print_word_frequence{
	my $fastafile = 'sample.dna';
	my @file_data = get_file_data($fastafile);
	my $sequence = extract_sequence_from_fasta_data(@file_data);

	# get count of tetramers
	my $size = 4;
	my %count = mercount($size, $sequence);

	# sort the keys by the count, and output results
	my @sortedkeys = sort {$count{$b} <=> $count{$a}} keys %count;
	print "total:".@sortedkeys."\n";

	foreach my $key (@sortedkeys){
		print "$key ", $count{$key}, "\n";
	}
}

############################################################### 
# print DNA subsequence frequence in GenBank library
# use mercount() count all 'mers'
############################################################### 
sub print_dna_frequence{
#sub print_dna_sequence_in_genbank_lib{
	my $genbank_lib='library.gb';
	my $fh ;

	die "can't open $genbank_lib" unless(open($fh,$genbank_lib));
	
	my %count = ();
	my $size = 4;

	while(my $record = get_next_fh_record($fh) ){
		# Extract DNA sequence
		# attention the 'undef' !!!
		my(undef,$dna) = get_annotation_and_dna($record);

		# Get count of tetramers
		%count = mercount($size, $dna);
	}

	close($fh);	

	# sort the count by the keys, and output the result
	my @sortedkeys = sort{ $count{$b} <=> $count{$a} } keys %count;

	# sort the counts by the key 'name' ,and output the result.
	#my @sortedkeys = sort keys %count;
	#print "total:".(keys %count)."\n";

	foreach my $key (@sortedkeys){
		print "$key ", $count{$key}, "\n";
	}
}

############################################################### 
# Mercount
# count all 'mers'
#  --subsequences of specified size -- in a sequence
############################################################### 
sub mercount{
	my($size, $seq) = @_;
	my %count = ();

	# iterate through each subsequence
	for(my $i = 0; $i < length($seq) -3; ++$i){
		my $mer = substr($seq, $i, $size);

		if(defined $count{$mer}){
			$count{$mer}++;
		}else{
			$count{$mer} = 1;
		}
	}
	return %count;
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
# parse 1
# extract the annotation and sequence sections
# from the first record of a Genbank library
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
# Parse annotation 
############################################################### 
sub parse_annotation{
	my($annotation) = @_;
	my(%results) = ();

	while($annotation =~ /^[A-Z].*\n(^\s.*\n)*/gm ){
#		print "-".$&."\n";
		
		my $value = $&;
		(my $key = $value) =~ s/^([A-Z]+).*/$1/s;

#		print "--".$1."\n";
#		print "key:$key\n";
#		print "value:$value\n";

		$results{$key} = $value;
	}
	return %results;
}

############################################################### 
# Parse features
############################################################### 
sub parse_features{
	my($features) = @_;
	my @features = ();

	#extract the features
	while($features =~ /^ {5}\S.*\n(^ {21}\S.*\n)*/gm){
		my $feature = $&;
		push @features, $feature;
	}
	return @features;
}

############################################################### 
# Parse features
# input : features from lib
# output : translations in a scalar
############################################################### 
sub parse_translation{
	my(@features) = @_;
	my $translation;

	foreach my $feature (@features){
		if($feature =~ /^.*\/translation=\"(.*)\"/s){
			$trans = $1; 
			$trans =~ s/\s*//g;
			$translation = $trans;
		}
	}
	#print "translation : \n$translation\n";
	return $translation;
}


############################################################### 
# extract annotation and sequence section from a genebank record
# extract the annotation and sequence section
# from the first record of a Genbank library.
# using regular expression to get two sections.
############################################################### 
sub extract_annotation_and_sequence_section{
	my $annotation = '';
	my $dna = '';
	my $record = '';
	my $filename = 'record.gb';
	my $save_input_separator = $/;

	# Open GenBank library file
	die "Can't open GenBank file:$filename!\n" unless(open(GBFILE,$filename));

	# Set input separator to "//\n" and read in a record to a scalar.
	$/ = "//\n";

	# Well,not <> anymore.
	$record = <GBFILE>;

	# Reset input separator
	$/ = $save_input_separator;

	# Separator the annotation from the sequence data.
	# maybe use # instead / is better.
	($annotation, $dna) = ($record =~ /^(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n/s);

	# print the two pieces
	print "annotation:\n$annotation \ndna:\n$dna\n";

}
############################################################### 
# get annotation and dna
# input  : filehandle
# output : annotation and dna
############################################################### 
sub get_annotation_and_dna{
	my($record) = @_;

	my($annotation) = '';
	my($dna) = '';

	#print $record;

	# separate the annotation from the sequence data
	# foolish mistake..................................................|.
	#($annotation, $dna) = ($record =~ /^(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n\s);

	($annotation, $dna) = ($record =~ /^(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n/s);

	# Clean the sequence of any whitespace or / characters	# will there have '/' ?
	#$dna =~ s/[\s\/]//g;	#'s///g'	# this is mistake ? 
	$dna =~ s/[\s0-9]//g;	#'s///g' 

	return ($annotation,$dna);
}

############################################################### 
# Search sequence
# - search sequence with regular expression
############################################################### 
sub search_sequence{
	my($sequence, $regexp) = @_;
	my(@locations) = ();
	while($sequence =~ /$regexp/ig){
		push(@locations , pos);		# pos.
	}
	return @locations;
}

############################################################### 
# Search annotation
# 	- search annotation with regular expression
############################################################### 
sub search_annotation{
	my($annotation, $regexp) = @_;
	my(@locations) = ();

	# note the /s modifier --matches any characters including newline
	while($annotation =~ /$regexp/isg){
		push(@locations, pos);
	}
	return @locations;
}

############################################################### 
# Now with library
############################################################### 
sub search_genebank_library{
	my $fh;
	my $record;
	my $dna;
	my $annotation;
	my $offset;
	my $library = 'library.gb';

	# Open
	$fh = open_file($library);

	$offset = tell($fh);


	while( $record = get_next_fh_record($fh) ){
		($annotation, $dna) = get_annotation_and_dna($record);

		if(search_sequence($dna, 'AAA[CG].')){
			print "Sequence found in record at offset $offset\n";
		}
		if( search_annotation($annotation, 'homo sapiens')){
			print "Annotation found in record at offset $offset\n";
		}

		$offset = tell($fh);
	}

	# remember close fh
	close($fh);
}
############################################################### 
# parse annotation's keys.
############################################################### 
sub get_annotation_keys{
	my $fh;
	my $record;
	my $dna;
	my $annotation;
	my %fields;
	my $library = 'library.gb';
	my @features;

	# open library 
	$fh = open_file($library);

	$record = get_next_fh_record($fh);

	# Parse the sequence and annotation
	($annotation, $dna) = get_annotation_and_dna($record);

	# Extract the fields of the annotation
	print "\n---parse 'ANNOTATION'----\n";
	%fields = parse_annotation($annotation);

	foreach my $key (keys %fields){
		print "**$key**\n";
		print $fields{$key};
	}

	print "\n---parse 'FEATURES'----\n";
	@features = parse_features($fields{'FEATURES'});

	#print the feature
	foreach my $feature (@features){
		# Extract the name of the feature 
		my($feature_name) = ($feature =~ /^ {5}(\S+)/);

		print "-------$feature_name--------\n";
		print $feature;
	}
}

############################################################### 
# get next record from filehandle
# it is simple to handle this:
# just changethe 'separator',get record, and change it back.
############################################################### 
sub get_next_fh_record{
	my($fh) = @_;
	my($offset);
	my $record = '';

	my $saved_input_separator = $/;

	$/ = "//\n";
	$record = <$fh>;
	$/ = $saved_input_separator;

	return $record;
}

############################################################### 
# make a DBM index of a GeneBank library
# and demostrate its use interactively
############################################################### 
sub db_demo_with_genebank{
	my $fh;
	my $record;
	my $dna;
	my $annotation;

	my %fields;
	my %dbm;
	my $ans;
	my $offset;
	my $library = 'library.gb';

	# open DBM, create if necessary
	die "Can't open DBM GB with 0644\n" unless(dbmopen(%dbm, 'GB',0644));

	# Parse GeneBank library , saving accession number and offset in DBM.
	$fh = open_file($library);
	$offset = tell($fh);

	while($record = get_next_fh_record($fh)){
		# Get accession field for this record.
		($annotation, $dna) = get_annotation_and_dna($record);

		%fields = parse_annotation($annotation); 

		my $accession = $fields{'ACCESSION'};

		# Extract just accession number 
		$accession =~ s/^ACCESSION\s*//;
		$accession =~ s/\s*$//;

		# Store the key/value of accession/offset
		$dbm{$accession} = $offset;

		# Get offset for next record
		$offset = tell($fh);
	}

	# Now interactively query the DBM database with accession numbers
	# to see associated records

	print "Here are the available accession numbers:\n";
	print join("\n", keys %dbm),"\n";

	print "Enter accession number :";

	while ($ans = <>){
		chomp $ans;
		if($ans =~ /^\s*q/){
			last;
		}
		$offset = $dbm{$ans};
		
		if(defined $offset){
			seek($fh, $offset, 0);
			$record = get_next_fh_record($fh);
			print $record;
		}else{
			print "do NOT have entry for accession number $ans\n";
		}

		print "Here are the available accession numbers:\n";
		print join("\n", keys %dbm),"\n";
		print "Enter accession number :";
	}

	# Remember close DBM
	dbmclose(%dbm);

	# and Remember close fh
	close($fh);
}

############################################################### 
# Parsing Genbank annotation using array
############################################################### 
sub parsing_genebank_annotation{
	my @genebank = ();
	my $locus = '';
	my $accession = '';
	my $organism = '';
	my $flag = 0;	# use flag to...

	# get genbank filedata.
	@genebank = get_file_data('record.gb');

	# Get some of the identifying information.
	foreach my $line (@genebank){
		if($line =~ /^LOCUS/){
			$line =~ s/^LOCUS\s*//;
			$locus = $line;
		}elsif($line =~ /^DEFINITION/){
			$line =~ s/^DEFINITION\s*//;
			$definition = $line;
			$flag = 1;
		}elsif($line =~ /^ACCESSION/){
			$line =~ s/ACCESSION\s*//;
			$accession = $line;
			$flag = 0;
		}elsif($flag){
			chomp($definition);
			$line =~ s/^\s*//;	
			$definition .= $line;
		}elsif($line =~ /^  ORGANISM/){
			$line =~ s/\s*ORGANISM\s*//;
			$organism = $line;
		}
	}

	print "*LOCUS*\n";
	print $locus;
	print "*DEFINITION*\n";
	print $definition;
	print "*ACCESSION*\n";
	print $accession;
	print "*ORGANISM*\n";
	print $organism;

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

############################################################### 
# Print sequence with length.
############################################################### 
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
		# add exit . for exception.
		exit;
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

############################################################### 
# find the frequency of occurrence of the adjacent amino acids
# coded in a DNA sequence; or in a GenBank library.
############################################################### 
sub print_amino_acid_frequence{
	#amino_acid_frequence_in_dna_sequence();
	amino_acid_frequence_in_genebank_lib();
}

sub amino_acid_frequence_in_dna_sequence{
	$fasta_file = 'sample.dna';
	@file_data = get_file_data($fasta_file);
	$dna = extract_sequence_from_fasta_data(@file_data);
	
	$peptide = translate_frame($dna,1);

	while(my $aa = get_user_input("Count neighbors of what amino acid?:")){
		my %count_adjacent_aa = ();
		while($peptide =~ /(.)$aa(.)/g){
			# save adjacent amino acids as $aa1 & $aa2
			my($aa1, $aa2) = ($1, $2);	

			# store $aa1
			if(defined $count_adjacent_aa{$aa1}){
				$count_adjacent_aa{$aa1}++;
			}else{
				$count_adjacent_aa{$aa1} = 1;
			}
			
			# store $aa2
		if(defined $count_adjacent_aa{$aa2}){
				$count_adjacent_aa{$aa2}++;
			}else{
				$count_adjacent_aa{$aa2} = 1;
			}
		}

		print "in this sequence , the neighbors of the amino acid $aa \n"
			." have the following frequency:\n";

		# sort the keys by the count ,and output results
		my @sortedkeys = sort {$count_adjacent_aa{$b} <=> $count_adjacent_aa{$a}} keys %count_adjacent_aa;
		foreach my $key (@sortedkeys) {
			print "$key " , $count_adjacent_aa{$key}, "\n";
		}
	}
}

sub amino_acid_frequence_in_genebank_lib{
	my $filename = 'library.gb';
	my $fh = open_file($filename);
	my $peptide;

	# well , considered the frame window of the dna2peptide, 
	# it is better to get ranslation from the annotation.
	# let do a sub for that.

	# get all dna from lib (the amount maybe very large)
	while( $record = get_next_fh_record($fh) ){
		($annotation, undef) = get_annotation_and_dna($record);
		#parse_translation($annotation);

		%fields = parse_annotation($annotation);
		@features = parse_features($fields{'FEATURES'});
		$peptide .= parse_translation(@features);
	}

	my %count;
	while(my $aa = get_user_input("Count neighbors of what amino acid?:")){
		%count = count_adjacent_amino_acid($peptide,$aa);
		print_sorted_keys(%count);
	}
}

############################################################### 
# sort and print hash
# input : hash
# output: sorted keys
############################################################### 
sub print_sorted_keys{
	my(%count) = @_;
	my @sortedkeys = sort{$count{$b} <=> $count{$a}} keys %count;
	foreach my $key (@sortedkeys){
		print "$key , $count{$key}\n";
	}
}

############################################################### 
# count adjacent amino acid
# input : peptide
# output: count_adjacent_amino_acid hash
############################################################### 
sub count_adjacent_amino_acid{
	my($peptide,$aa) = @_;
	my %count_adjacent_aa = ();

	while($peptide =~ /(.)$aa(.)/g){
		# save adjacent amino acids as $aa1 & $aa2
		my($aa1, $aa2) = ($1, $2);
		
		# store $aa1
		if(defined $count_adjacent_aa{$aa1}){
			$count_adjacent_aa{$aa1}++;
		}else{
			$count_adjacent_aa{$aa1} = 1;
		}

		# store $aa2
		if(defined $count_adjacent_aa{$aa2}){
			$count_adjacent_aa{$aa2}++;
		}else{
			$count_adjacent_aa{$aa2} = 1;
		}
	}
	
	return %count_adjacent_aa;
}

############################################################### 
# Extract all the words from the annotation of a library
# of GenBank records.
# For each word found , add the offset of the GenBank record
# in the library to a DBM file that has keys equal to the words
# and values that are strings with offset separated by spaces.
# one key can have a space-separated list of offsets for a value  
############################################################### 
sub extract_words_from_genbank_lib{
	my $filename= 'library.gb';
	my $fh;

	die "can't open library $filename" unless(open($fh,$filename));
	
	# A hash to store the byte offsets of records.
	# keyed by a word appearing in the record
	my %words = ();

	# A clearly inadequate list of words to ignore
	my %ignore = (
		'and'=> 1,
		'or' => 1,
		'organism' => 1,
		'id' => 1,
		'if' => 1,
		'species' => 1,
	);

	# Get the GenBank records and their byte offsets
	for( my $byte_offset = tell($fh); my $record = get_next_fh_record($fh);
		$byte_offset = tell($fh)){
		# Extract annotation
		my( $annotation, undef) = get_annotation_and_dna($record);
		
		# Only need to handle each word once   per GenBank record
		my %seen = ();

		# Extract words, saving the byte offset for found words
		# what's the definition of a word for GenBank annotations?
		while($annotation =~ /(\w[\w'-]*)/g){	#\w and \w'-   ???
			# Store everything in lowercase
			my $word = lc $1;
			
			# Skip unwanted words , or words already found in this record
			defined $ignore{$word} and next;
			defined $seen{$word} and next;

			# Mark this new word as seen
			$seen{$word} = 1;

			# Add the byte offset of this record to the value for this word in hash
			if(defined $words{$word}){
				$words{$word} .= " $byte_offset";
			}else{
				$words{$word} = $byte_offset;
			}
		}
	}
	
	# Interact with the user, asking for words and showing the GenBank records
	# containing them
	while(my $query = lc get_user_input("search for what word? (space for quit): ")){
		# ask again if the requested word doesn't appare in the lib
		unless(defined $words{$query}){
			print "The word \"$query\" is not in the lib $filename\n";
			next;
		}

		# Make an array of the byte offsets of the GenBank records
		# containing the query word
		my @offsets = split " ", $words{$query};
			
		print "there are ". @offsets ." record. ,  >>offsets : @offsets\n";

		# If the user wants to see any of the records for the requested work.
		if(get_user_input("Display the records containing that word? (y or n) ")
			=~ /^\s*y/i){
			# Display each GenBank record beginning at a saved byte offset in the lib
			do{
				my $offset = shift @offsets;
				
				# Point the filehandle at the offset
				seek($fh, $offset, 0);

				# Print the record at that offset
				print get_next_fh_record($fh);
				if(@offsets > 0){
					print "there are ". @offsets ." record,  >>offsets : @offsets\n";
				}
			# If there is another ,and if ther user wnat to see it ,loop agin
			}while(@offsets and (get_user_input("Would u like see the next record? (y or n): ") =~ /^\s*y/i));
		}
	}
	close($fh);

}

############################################################### 
# Extract sequence chains from PDB file
############################################################### 
sub extract_sequence_chains_from_pdb_file{
	# Read in PDB file ( some file may very large! )
	# and attention the '1' and 'l'
	my @file = get_file_data('pdb/c1/pdb1c1f.ent');
	
	# Parse the record types of the PDB file
	my %record_types = parse_pdb_record_types(@file);
	
	# Extract the amino acid sequences of all chains in the protein
	my @chains = extract_seqres($record_types{'SEQRES'});
	
	# Translate the 3-character codes to 1-character codes,and print.
	foreach my $chain (@chains){
	#	print "***chain $chain***\n\n";
	#	print "$chain\n";
		print iub3to1($chain),"\n";
	}
}

############################################################### 
# Extract secondary structure information contained in the
# HELIX, SHEET, TURN record types of PDB files.
# Print out the secondary structure and the primary sequence together,
# So that is't easy to see by what secondary structure a given
# residue is included.
# ( Consider using a special alphabet for secondary structure,
#   so that every residue in a helix is represented by H.
#   The main complication involves keeping track of chain names
#   for these lines and for the SEQRES lines: we also have to
#   modify the parseSEQRES subroutine for this.
#   Note that if there is only one chain, the chain name is often
#   left blank.)
############################################################### 
sub extract_secondary_structure_from_pdb_file{
	my $pdb_file = 'pdb/44/pdb244l.ent';
	
	# Get the file data
	my @pdb_file_data = get_file_data($pdb_file);

	# Parse the record types
	my %record_types = parse_pdb_record_types(@pdb_file_data);

	# Extract the name and the sequence of each chain
	my %chains = extract_seqres_2($record_types{'SEQRES'});

	foreach my $key (keys %chains){
		# Change the sequence from 3-letter to 1-letter amino acid codes
		$chains{$key} = iub3to1($chains{$key});

		# for dbg
		print "chains $key:$chains{$key}\n";
	}

	# Parse the HELIX , SHEET, and TURN record types.
	# The following three statement use the ternary operator .
	$record_types{'HELIX'}?my %helix = extractHELIX($record_types{'HELIX'}):();
	$record_types{'SHEET'}?my %sheet = extractSHEET($record_types{'SHEET'}):();
	$record_types{'TURN'} ?my %turn  = extractTURN ($record_types{'TURN'}) :();

	# for dbg
	print "helix:\n";
	foreach my $key (sort keys %helix){
		print "$key:$helix{$key}\n";
	}
	print "sheet:\n";
	foreach my $key (sort keys %sheet){
		print "$key:$sheet{$key}\n";
	}
	print "turn:\n";
	foreach my $key (sort keys %turn){
		print "$key:$turn{$key}\n";
	}

	# Now make a annotation strings that contain the helix, sheet and turn.
	my %annotation = ();
	foreach my $chain_name (sort keys %chains){
		# Make a blank annotation string, 
		# same length as the sequence of the chain
		$annotation{$chain_name} = ' ' x length($chains{$chain_name});

		if( defined $helix{$chain_name}){
			foreach my $structure (split /:/, $helix{$chain_name}){
				print "structure1:$structure\n";
				
				my($structure, $position) = split /;/, $structure;

				print "structure2:$structure\n";
				#dbg
				print "annotation:$annotation{$chain_name}".
				"position :".($position-1).
				"length :".length($structure)."\n";

				substr($annotation{$chain_name}, $position-1, 
					length($structure)) = $structure;
				print "helix ok\n";
				
			}
		}
		if( defined $sheet{$chain_name}){
			foreach my $structure (split /:/, $sheet{$chain_name}){
				my($structure, $position) = split /;/, $structure;
				substr($annotation{$chain_name}, $position-1, 
					length($structure)) = $structure;
				
			}
		}
		if( defined $turn{$chain_name}){
			foreach my $structure (split /:/, $turn{$chain_name}){
				my($structure, $position) = split /;/, $structure;
				substr($annotation{$chain_name}, $position-1, 
					length($structure)) = $structure;
				
			}
		}

		print "Chain name =\" $chain_name\"\n";
		print $chains{$chain_name},"\n";
		print "$annotation{$chain_name}\n";
	}


}
############################################################### 
# Extract HELIX
############################################################### 
sub extractHELIX{
	my($helix) = @_;

	# make array of lines
	my @record = split( /\n/, $helix);

	# hash to store chains
	my %chain_hash = ();

	foreach my $line (@record){
print "$line\n";
		# Chain is in column 20, starting position in column 22-25
		# length is in column 72-76
		# "strip_space" removes leading and trailing spaces.
		# IF there's only one chain in the PDB entry,the chain name maybe blank.

		my($this_chain) = strip_space(substr($line,19,1));
		my($start) 	= strip_space(substr($line,21,4));
		my($length) = strip_space(substr($line,71,5));
print "this_chain:$this_chain,start:$start,length:$length\n";

		if(defined $chain_hash{$this_chain}){
			$chain_hash{$this_chain} .= ':'. 'H' x $length. ";$start";
		}else{
			$chain_hash{$this_chain}  =      'H' x $length. ";$start";
		}
	}
	return %chain_hash;
}

############################################################### 
# Extract SHEET
############################################################### 
sub extractSHEET{
	my($sheet) = @_;

	# make array of lines
	my @record = split( /\n/, $sheet);

	# hash to store chains
	my %chain_hash = ();

	foreach my $line (@record){
		# For sheet,
		# Chain is in column 22, starting position in column 23-26
		# end position is in column 34-37
		# "strip_space" removes leading and trailing spaces.
		# IF there's only one chain in the PDB entry,the chain name maybe blank.

		my($this_chain) = strip_space(substr($line,21,1));
		my($start) 			= strip_space(substr($line,22,4));
		my($end) 				= strip_space(substr($line,33,4));
		my($length) 		= $end - $start +1;

		if(defined $chain_hash{$this_chain}){
			$chain_hash{$this_chain} .= ':'. 'S' x $length. ";$start";
		}else{
			$chain_hash{$this_chain}  =      'S' x $length. ";$start";
		}
	}
	return %chain_hash;
}

############################################################### 
# Extract TURN
############################################################### 
sub extractTURN{
	my($turn) = @_;

	# make array of lines
	my @record = split( /\n/, $turn);

	# hash to store chains
	my %chain_hash = ();

	foreach my $line (@record){
		# Chain is in column 20, starting position in column 21-24
		# length is in column 32-35
		# "strip_space" removes leading and trailing spaces.
		# IF there's only one chain in the PDB entry,the chain name maybe blank.

		my($this_chain) = strip_space(substr($line,19,1));
		my($start) 			= strip_space(substr($line,20,4));
		my($end) 		= strip_space(substr($line,31,4));
		my($length) = $end - $start +1;

		if(defined $chain_hash{$this_chain}){
			$chain_hash{$this_chain} .= ':'. 'T' x $length. ";$start";
		}else{
			$chain_hash{$this_chain}  =      'T' x $length. ";$start";
		}
	}
	return %chain_hash;
}


############################################################### 
# Extract sequence chains from PDB file
############################################################### 
sub extract_atomic_coordinates_from_pdb_file{
	# Read in PDB file ( some file may very large! )
	# and attention the '1' and 'l'
	my @file = get_file_data('pdb/c1/pdb1c1f.ent');
	
	# Parse the record types of the PDB file
	my %record_types = parse_pdb_record_types(@file);
	
	# Extract the amino acid sequences of all chains in the protein
	my %atoms = parse_atom($record_types{'ATOM'});

	# Print out a couple of the atoms.
	print $atoms{'1'},"\n";
	print $atoms{'1078'},"\n";
	
}

############################################################### 
# Extract sequence chains from PDB file
############################################################### 
sub parse_atom{
	my($atom_record) = @_;
	my %results = ();
	
	my(@atom_record) = split(/\n/,$atom_record);

	foreach my $record (@atom_record){
		my $number = substr($record, 6, 5);	# columns 7-11
		my $x = substr($record, 30, 8);	# columns 31-38
		my $y = substr($record, 38, 8);	# columns 39-46
		my $z = substr($record, 46, 8);	# columns 47-54
		my $element  = substr($record, 76, 2);	# columns 77-78
		
		# $number and $element may have leading spaces,strip them
		$number =~ s/^\s*//;
		$element =~ s/^\s*//;

		# Store information in hash
		$results{$number} = "$x $y $z $element";
	}
	return %results;
}

############################################################### 
# Extract sequence chains from PDB file
# Input : an array of a PDB file,
# output: a hash with 
#			key = record type names
#			value = scalar containing lines for that record type
############################################################### 
sub parse_pdb_record_types{
	my @file = @_;
	my %record_types = ();

	foreach my $line (@file){
		# Get the record type name which begins at the
		# start of the line and ends at the first space
		my($record_type) = ($line =~ /^(\S+)/);
		
		if(defined $record_types{$record_type}){
			$record_types{$record_type} .= $line;
		}else{
			$record_types{$record_type} = $line;
		}
	}
	return %record_types;
}

############################################################### 
# Parse out the record types of a PDB file using regular expressions
# instead of iterating through an array of input lines.
# 
# This maybe not the best way to do this job.
# Since each line of PDB files have the record type as first field.
# It is very easy to just iterate through the lines and save the types
# by looking at the field.
# It is also space-efficient, since you only have to read in one line
# at a time from the file, and some  PDB  files are very large.
# Nevertheless, as an exercise, this approach given here maybe useful.
############################################################### 
sub parse_pdb_record_types_by_re{
	# Open a PDB file and read it into a scalar variable
	#my $pdbfile = 'pdb/pdb1a4o.ent';
	my $pdbfile = 'pdb/44/pdb244l.ent';

	open(PDBFILE, "$pdbfile") or die("can't open pdbfile:$pdbfile\n");

	# It's nice to warn users if the program will seen unresponsive for a while
	print "Reading PDB file $pdbfile ... this may take a minute ...";

	my $pdb_data = join('', <PDBFILE>);

	print "\n";		# smart ~

	# Parse the record types into a hash data structure
	my(%record_types) = ();
	
	# Note the regular expression:
	# ([A-Z]+) finds the name of the record type(without trailing digits if any)
	# [^\n]*\n find any number of non-newlines, up to the first newline

	while( $pdb_data =~ /([A-Z]+)[^\n]*\n/gs){

		# $1 matches the parenthesized part of our regular expression
		# $& matches the entire matched string of our regular expression

		if( defined $record_types{$1} ){
			$record_types{$1} .= $&;
		}else{
			$record_types{$1} = $&;
		}
	}

	# Interact with the user, asking for the record types and showing the lines
	# comprising them.
	# 
	# This while loop has two statements separated by a comma.
	# such a list of statement will return the value of the last item 
	# in the list.
	while(print("The record types for this file are:\n",
				join(" ", sort keys %record_types), "\n"),
			my $query = get_user_input("Show which record type?: ")){

			if(defined $record_types{$query}){
				print $record_types{$query},"\n";
			}else{
				print "The record type \"$query\" is not in the PDB files\n";
			}
	}

}


############################################################### 
# the function parse_and_search_pdb_files will use global variable
############################################################### 
my @global_pdb_files=();

############################################################### 
# Parse HEADER, TITLE, KEYWORDS record types of all PDB files.
# Make a hash with key as a word from those record types,
# and value as a list of filenames that contained that word.
############################################################### 
sub parse_and_search_pdb_files{

# Use File::Find to locate all PDB files.
use File::Find;

	# Search the directory 'pdb' for pdb files, 
	# saving their names in the @pdbfiles array.
	my @pdb_files = ();

	find ( \&pdbfiles, ('pdb')   );

#print "pdb_files >>> :@global_pdb_files\n";
	
	# A hash to store the filenames of records
	# keyed by a word appearing in the record
	my %words = ();

	# A clearly inadequate list of words to ignore.
	my %ignore = (
		'and' 	=> 1,
		'or' 	=> 1,
		'organism'=>1,
		'id'	=>1,
		'if'	=>1,
		'species'=>1,
		'header'=>1,
		'title' =>1,
		'keywords'=>1,
	);

	foreach my $pdffile (@global_pdb_files){
		# Read in PDB files.
		my @file = get_file_data($pdffile);

		# Parse the record types of the PDB file
		my %record_types = parse_pdb_record_types(@file);

		# Put all the three desired record types into one scalar
		my $annotation = '';
		(defined $record_types{'HEADER'}) and ($annotation .= $record_types{'HEADER'});
		(defined $record_types{'TITLE'}) and ($annotation .= $record_types{'TITLE'});
		(defined $record_types{'KEYWORDS'}) and ($annotation .= $record_types{'KEYWORDS'});

		# Only need to handle each word once per PDB record
		my %seen = ();

		# Extract words, saving the byte offset for found words.
		# !! what \w'- mean ??
		while($annotation =~ /(\w[\w'-]*)/g){
			# Store everything in lowercase
			my $word = lc $1;

			# skip unwanted words, or words already founded
			defined $ignore{$word} and next;
			defined $seen{$word} and next;

			# Mark this new word and seen
			$seen{$word} = 1;

			# Add the filename of this record
			# to the value for this word in the hash
			if(defined $words{$word}){
				$words{$word}.=" $pdffile";
			}else{
				$words{$word} = $pdffile;
			}
		}
	}

	# Interact with the user,
	# asking for words and showing the names of PDB files contain them
	while(my $query = lc get_user_input("Find the PDB files' word?: ")){
		if(defined $words{$query}){
			print $words{$query},"\n";
		}else{
			print "The word \"$query\" is not in the PDB files\n";
		}
	}
}

############################################################### 
# Find all files whoes name begin with 'pdb' and end with '.ent'
# Warning :this will use a global
# output: 0
############################################################### 
sub pdbfiles{
	# Ignore files that aren't ASCII text files or aren't readable
	-T and -r or return 0;

	/^pdb.*\.ent$/ and push(@global_pdb_files,"$File::Find::name");

	#print "pdb_files : @pdb_files\n";
	return 0;
}

############################################################### 
# Extract SEQRES
# Input : scalar containing SEQRES lines.
# output: an array containing the chains of the sequence.
############################################################### 
sub extract_seqres{
	my($seqres) =@_;
	my $last_chain = '';
	my $sequence = '';
	my @results = ();

	# Make array of lines
	my @records = split(/\n/, $seqres);
		
	foreach my $line (@records){
		# Chain is in column 12,(A,B,C,D) residues start in column 20.
		my $this_chain = substr($line, 11, 1);
		my $residues = substr($line, 19, 52);	# space at end
		
		# Check if new chain, or continuation of previous chain
		if("$last_chain" eq ""){											# new chain
			$sequence = $residues;	# new chain
		}elsif("$this_chain" eq "$last_chain"){				# same chain
			$sequence .= $residues;
		}elsif($sequence){				# Finish gathering previous chain 
								# if the have one and just have one chian ,won't use this.
			push(@results, $sequence);
			$sequence = $residues;
		}
		$last_chain = $this_chain;
	}
	# Save last chain.
	push(@results, $sequence);
	return @results;
}
############################################################### 
# stripe space.
############################################################### 
sub strip_space{
	my($string) = @_;
	$string =~ s/^\s*//;
	$string =~ s/\s*$//;
	return $string;
}

############################################################### 
# Extract SEQRES 2
# Modified version that reports the chain name plus the sequence.
# returning the values as a hash.
#
# Input : scalar containing SEQRES lines.
# output: an array containing the chains of the sequence.
############################################################### 
sub extract_seqres_2{
	my($seqres) =@_;
	my $last_chain ;
	my $sequence = '';
	my %results = ();			# result changed to hash.

	# Make array of lines
	my @records = split(/\n/, $seqres);
		
	foreach my $line (@records){
		# Chain is in column 12,(A,B,C,D) residues start in column 20.
		my $this_chain = strip_space(substr($line, 11, 1));	# need stripe_sapce?
		my $residues = substr($line, 19, 52);	# space at end
		
		# Check if new chain, or continuation of previous chain
		if(not defined $last_chain ){											# if NOT defined lastchain
			$sequence = $residues;											# new chain
		}elsif("$this_chain" eq "$last_chain"){				# same chain
			$sequence .= $residues;
		}elsif($sequence){				# Finish gathering previous chain
															# (unless first chain) 
								# if the have one and just have one chian ,won't use this.
			$results{$last_chain}= $sequence;
			$sequence = $residues;
		}
		$last_chain = $this_chain;
	}

	# Save last chain.
	$results{$last_chain}= $sequence;
	return %results;
}

############################################################### 
# IUB 3 to 1
# change string of 3-character IUB amino acid codes
# (white space separated) into a string of 1-character aa.
############################################################### 
sub iub3to1{
	my($input) = @_;
	my %three2one =(
		'ALA' => 'A',
		'VAL' => 'V',
		'LEU' => 'L',
		'ILE' => 'I',
		'PRO' => 'P',
		'TRP' => 'W',
		'PHE' => 'F',
		'MET' => 'M',
		'GLY' => 'G',
		'SER' => 'S',
		'THR' => 'T',
		'TYR' => 'Y',
		'CYS' => 'C',
		'ASN' => 'N',
		'GLN' => 'Q',
		'LYS' => 'K',
		'ARG' => 'R',
		'HIS' => 'H',
		'ASP' => 'D',
		'GLU' => 'E',
	);

	# clean up the input
	$input =~ s/\n/ /g;

	my $seq = '';

	# this use of split separates on any contiguous whitespace
	my @code3 = split(' ', $input);

	foreach my $code (@code3){
		# A little erroe checking
		if(not defined $three2one{$code}){
			print "Code $code not defined.\n";
			next;
		}
		$seq .= $three2one{$code};
	}
	return $seq;
}

############################################################### 
# call stride for secondary structure prediction
############################################################### 
sub call_stride_for_secondary_structure_prediction{
	# Call "stride" on a file, collect the report
	my(@stride_output) = call_stride('pdb/c1/pdb1c1f.ent');
	
	# Parse the stride report into primary sequence, 
	# and secondary structure prediction
	my($sequence, $structure) = parse_stride(@stride_output);

	# Print out the beginnings of the sequence 
	# and the secondary structure
	print substr($sequence, 0 ,80), "\n";
	print substr($structure, 0 ,80), "\n";
}

############################################################### 
# call stride
# input : a PDB filename,
# output: output from the stride 
#   	stride : secondary sturcture prediction program
############################################################### 
sub call_stride{
	my($filename) = @_;

	# the stride program options
	my($stride) = '/usr/local/bin/stride';	# modify to urs.
	my $options = '';
	my @results = ();

	# Check for presence of PDB file
	die "File \" $filename\" doesn't exist\n" unless(-e $filename);

	# Start up the program ,capture and return the output
	@results = `$stride $options $filename`;
	return @results;
}

############################################################### 
# parse stride
# extract the primary sequence and the secondary structure prediction
# input : stride's  output , 
# output: two-element array.
############################################################### 
sub parse_stride{
	my(@stride_report) = @_;
	my $seq = '';
	my $str = '';
	my $length;

	# Extract the line of interest
	my(@seq) = grep(/^SEQ/,@stride_report);
	my(@str) = grep(/^STR/,@stride_report);

	# Process those lines to discard all 
	# but the sequence or structure information
	for(@seq) { $_ = substr($_, 10, 50)};
	for(@str) { $_ = substr($_, 10, 50)};

	# Return the information as an array of two strings
	$seq = join('', @seq);
	$str = join('', @str);

	# Delete unwanted spaces from the ends of the strings.
	# ($seq has no spaces that are wanted, but $str may)
	$seq =~ s/(\s+)$//;
	$length = length($1);
	$str =~ s/\s{$length}$//;
	
	return( ($seq,$str));

}

############################################################### 
# get user input
############################################################### 
sub get_user_input{
	my($prompt)= @_;
	print $prompt;

	my $ans = <>;
	chomp $ans;

	if($ans =~ /^\s*$/ or  $ans =~ /^\s*quit\s*$/i){
		return ''; 	# return null if response is empyt , q , or quit.
				#                                    ^
				# not include q,because q is an amino acie."
	}else{
		$ans =~ s/^\s*//;
		$ans =~ s/\s*$//;
		return $ans;
	}	
}

#---------------------------------------------------------------------------#
# File operation							    #	
#---------------------------------------------------------------------------#

############################################################### 
# Demonstrationg how to open a floder and list its contents.
############################################################### 
sub list_contents_of_folder{
	my @files = ();
	my($folder)= @_;

	# Open folder
	die "can't open folder $folder\n" unless(opendir(FOLDER, $folder));	
	
	# Read the contents of the folder ( the files and subfolders)
	@files = readdir(FOLDER);

	# Close the folder
	closedir(FOLDER);

	# Use grep to filter out the array entries '.' and '..'
	@files = grep(!/^\.\.?$/, @files);

	# Or use the 'combinded two lines' statement.
	# @files = grep(!/^\.\.?$/, readdir(FOLDER));

	# Print single level folder's files.
	#	print join("\n", @files), "\n";

	# If file , print filename
	# If folder , print its name and contents
	# Notice we need to prepend the folder name.
	foreach my $file (@files){
		if(-f "$folder/$file"){
			print "$folder/$file\n";
		}elsif(-d "$folder/$file"){
			
			my $folder = "$folder/$file";

			# Open the subfolder and list its contents
			list_contents_of_folder($folder);
		}
	}

}

############################################################### 
# find the oldest and largest file on the hard drive
############################################################### 
sub find_oldest_and_largest{
use File::Find;
	print "home : ". $ENV{'HOME'}."\n";
	find( \&find_by_file_find, $ENV{'HOME'} );
}

############################################################### 
# Report all files greater than 10 Mb and older than a year
############################################################### 
sub find_by_file_find{
	-f
	and ( -s > 10000000)
	and ( -A > 365)
	and ( print $File::Find::name," ", -s, " bytes ", -A, " days old\n\n");
}

############################################################### 
# find all perl files within File:Find
############################################################### 
sub find_all_perl_files{

# Use File::Find again
use File::Find;

	find( \&isperl, ($ENV{'HOME'}) );
	#find( \&find_perl_files_by_pl, ($ENV{'HOME'}) );
}

############################################################### 
# report all perl files within File:Find
############################################################### 
sub isperl{
	# Here is a method that finds command interpretation lines.
	# like '#!/usr/bin/perl'

	# Ignore files that aren't ASCII text files or aren't readable
	-T and -r or return 0;

	# Open the file and see if the first line is a command interpreter line
	open(THISFILE,$_) 
	or (print "$File::Find::name :can't open\n") 
	and return 0;
	
	my $firstline = <THISFILE>;
	close THISFILE ;

	$firstline or return 0;

	($firstline =~ /^#\!.*perl/) 
	and (print $File::Find::name,"\n") 
	and (return 1);

	return 0;
}

############################################################### 
# report all perl files by pl
############################################################### 
sub find_perl_files_by_pl{
	# if doesn't use command interpretation on os
	# this method will looks for the filename extention ".pl"

	# Ignore files that aren't ASCII text files or aren't readable
	-T and -r or return 0;

	/\.pl$/ and (print $File::Find::name,"\n") and (return 1);

	return 0;

}




1
