package Xquest_Digest;
use strict; 
#---------------------------------------------------------------------------
# Module: Xquest_Digest.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for petide digest.
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Licence
# This software is licensed under the Apache License, Version 2.0.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# This software and associated documentation is provided "AS IS",
# without warranties or conditions of any kind.
# See the License for the specific language governing permissions and
# limitations under the License.
#---------------------------------------------------------------------------
use Data::Dumper;

sub getpeps_gen{
my $MSTAB   = shift;
my $PARAMS  = shift;
my $ENZ     = shift;
my $seqobj  = shift;
my $verbose = shift;
my $enzymenumber=shift;

### Modify the Params Hash
$PARAMS->{'enzyme_num'}=$enzymenumber;

my $seqs = getpeps($MSTAB, $PARAMS, $ENZ, $seqobj, 1);
return $seqs;
}

sub getpeps {
	my $MSTAB   = shift;
	my $PARAMS  = shift;
	my $ENZ     = shift;
	my $seqobj  = shift;
	my $verbose = shift;

	my ( $i, $j, @sequences, $mincuts, @cuts );
	my $missed_cleav = $PARAMS->{'missed_cleavages'};

	my @enzymes = split /,/, $PARAMS->{'enzyme_num'};

	@cuts = ();
	if ( $PARAMS->{'define_enzyme'} ) {
		push @cuts, _definedEnzyme( $seqobj, $PARAMS );
	}
	else {
		foreach my $enzymenum (@enzymes) {
			if ( $enzymenum == 1 ) {
				push @cuts, _digestTryps( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 2 ) {
				push @cuts, _digestChymoTryps( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 3 ) {
				push @cuts, _digestGluC( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 9 ) {
				push @cuts, _digestTrypsR( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 14 ) {
				push @cuts, _digestLysC( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 10 ) {
				push @cuts, _digestAspN( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 15 ) {
				push @cuts, _digestTrypsAspN( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 16 ) {
				push @cuts, _digestChymoTryps_lowspec( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 17 ) {
				push @cuts, _digestTryps_lowspec( $seqobj, $ENZ );
			}
			elsif ( $enzymenum == 18 ) {
				push @cuts, _digestTryps_highspec( $seqobj, $ENZ );
			}elsif ( $enzymenum == 19 ){
				push @cuts, _digest_semiTryp($seqobj, $ENZ);
			}elsif ( $enzymenum == 20 ){
				push @cuts, _digest_LysN($seqobj, $ENZ);
			}elsif ($enzymenum == 21){
				push @cuts, _digest_mixed($seqobj, $ENZ);				
			}elsif ($enzymenum == 22){
				push @cuts, _digest_chymotrypsin3($seqobj, $ENZ);				
			}elsif ($enzymenum == 23){
				push @cuts, _digest_chymotrypsin4($seqobj, $ENZ);				
			}elsif($enzymenum == 99){
			    push @cuts, _digestTryps_silac($seqobj, $ENZ);
			}
			
			else {
				die "no digest method defined for enzyme ",
				  $PARAMS->{'enzyme_num'}, "$!";
			}
		}
	}

	push @cuts, $seqobj->length + 1;
	push @cuts, 1;

	my %seen = ();
	@cuts = grep { !$seen{$_}++ } @cuts;
	@cuts = sort { $a <=> $b } @cuts;
	my $allpeps = $PARAMS->{'reportimpossiblepeptides'};
	$mincuts = $PARAMS->{'requiredmissed_cleavages'};
	my $singlytryptic = $PARAMS->{'allow1misscleavage'};
	
	if ( $PARAMS->{'printdigestpeps'} ) {
		print "@cuts\n";
	}
	
	my $k = 1;
	my %seenseq;
	for $i ( 0 .. $#cuts - 1 ) {
		my $lastcut = 1 + $mincuts;
		for $j ( 1 + $mincuts .. $missed_cleav + 1 ) {

			if ( ( $i + $j ) <= $#cuts ) {

				if ( $PARAMS->{'mindigestlength'} <=
					   ( ( $cuts[ $i + $j ] - 1 ) - $cuts[$i] )
					&& ( ( $cuts[ $i + $j ] - 1 ) - $cuts[$i] ) <=
					$PARAMS->{'maxdigestlength'} )
				{
					$lastcut = $j;
					my $id = join '::', $seqobj->id, $cuts[$i], $cuts[ $i + $j ] - 1;
					
					#print "Id: $id\n";
					
					#### CONSIDERING MODIFIED NTERMINUS 
					#### ADDED 29/3/2010, by wathomas
					my $sequence_modified;
					if (($PARAMS->{'ntermxlinkable'}) && ($i==0) || $PARAMS->{'ntermmodified'} ){
					#print "Ntermlinkable option! modify nterminus.\n";
					$sequence_modified = getmodified_seq($seqobj->subseq( $cuts[$i], $cuts[ $i + $j ] - 1 ), $PARAMS, $id,1 );    #expand into variable modifications
					}else{
					$sequence_modified = getmodified_seq($seqobj->subseq( $cuts[$i], $cuts[ $i + $j ] - 1 ), $PARAMS, $id );    #expand into variable modifications	
					}
					
					foreach my $xids ( keys %$sequence_modified ) {
						my $sequence = $sequence_modified->{$xids}->{'seq'};
						if ( test_xlinkrequirements( $sequence, $PARAMS, $xids )
							|| $allpeps )
						{
							if ( $PARAMS->{'printdigestpeps'} ) {
								print $k++, " ", "$xids $sequence\n";
							}
							my $desc = join '::', $cuts[$i],
							  $cuts[ $i + $j ] - 1;
							my $pepobj = PepObj->new( $sequence, $xids, $desc, $MSTAB, $PARAMS, $verbose, );
							push @sequences, $pepobj;
							$seenseq{$id}++;
						}

					}

				}
			}
		}
	}
	if ($singlytryptic) {
		for $i ( 0 .. $#cuts - 1 ) {
			for $j ( 1 + $mincuts .. $missed_cleav + 1 ) {
				my $nend = $cuts[$i];
				my $cend = $cuts[ $i + $j ] - 1;

				#print "$nend ... $cend\n";

				while (
					   ( $cend > $nend )
					&& ( $cend - $nend ) >= $PARAMS->{'mindigestlength'}

				  )
				{

					$cend--;
					if ( ( $cend - $nend ) <= $PARAMS->{'maxdigestlength'} ) {
						my $coorid = join '::', $seqobj->id, $nend, $cend;
						my $id = join '::', $seqobj->id, $nend, $cend, "nt";
						unless ( $seenseq{$coorid}++ ) {

							my $sequence_modified =
							  getmodified_seq( $seqobj->subseq( $nend, $cend ),
								$PARAMS, $id )
							  ;    #expand into variable modifications

							foreach my $xids ( keys %$sequence_modified ) {
								my $sequence =
								  $sequence_modified->{$xids}->{'seq'};
								if (
									test_xlinkrequirements( $sequence, $PARAMS,
										$xids )
									|| $allpeps
								  )
								{
									if ( $PARAMS->{'printdigestpeps'} ) {
										print $k++, " nt: ",
										  "$id $coorid $sequence\n";
									}
									my $desc = join '::', $nend, $cend;
									my $pepobj = PepObj->new(
										$sequence, $xids,
										$desc,     $MSTAB,
										$PARAMS,   $verbose,

									);

									push @sequences, $pepobj;
								}
							}

						}

					}
				}    #end while ($cend ...
				$cend = $cuts[ $i + $j ] - 1;
				$nend = $cuts[$i];
				while (
					   ( $cend > $nend )
					&& ( $cend - $nend ) >= $PARAMS->{'mindigestlength'}

				  )
				{

					$nend++;
					if ( ( $cend - $nend ) <= $PARAMS->{'maxdigestlength'} ) {
						my $coorid = join '::', $seqobj->id, $nend, $cend;
						my $id = join '::', $seqobj->id, $nend, $cend, "nt";
						unless ( $seenseq{$coorid}++ ) {
							my $sequence_modified =
							  getmodified_seq( $seqobj->subseq( $nend, $cend ),
								$PARAMS, $id )
							  ;    #expand into variable modifications

							foreach my $xids ( keys %$sequence_modified ) {
								my $sequence =
								  $sequence_modified->{$xids}->{'seq'};
								if (
									test_xlinkrequirements( $sequence, $PARAMS,
										$xids )
									|| $allpeps
								  )
								{

									if ( $PARAMS->{'printdigestpeps'} ) {
										print $k++, " nt: ",
										  "$xids $coorid $sequence\n";
									}
									my $desc = join '::', $nend, $cend;
									my $pepobj = PepObj->new(
										$sequence, $xids,
										$desc,     $MSTAB,
										$PARAMS,   $verbose,

									);

									push @sequences, $pepobj;
								}
							}

						}

					}
				}    #end while ($cend ...

			}
		}
	}
	return \@sequences;
}

sub test_xlinkrequirements {
	my $sequence   = shift;
	my $PARAMS     = shift;
	my $xids       = shift;
	my @xlinksites = ();

	while ( $sequence =~ /$PARAMS->{'AArequired'}/gi ) {
		push @xlinksites, pos($sequence);
	}

	my $nrxlinksites = $#xlinksites + 1;

	if ( $nrxlinksites == 0 ) {

		#	$superverbose && print "$xids $sequence -> no xlinksite discarded\n";
		return 0;
	}
	elsif ( $nrxlinksites > 1 ) {

 #	$superverbose &&  print "$xids $sequence nr of xlinks = ",$nrxlinksites,"\n";
		return 1;
	}

	elsif (( $nrxlinksites == 1 )
		&& ( $xlinksites[0] < length($sequence) ) )
	{    #only one K but in the middle of peptide
		 #	print "$xids $sequence nr sites only ",$nrxlinksites," but in the middle at:",$xlinksites[0],"\n";
		return 1;
	}
	elsif ($PARAMS->{'nocutatxlink'}
		&& ( $nrxlinksites == 1 )
		&& ( $xlinksites[0] == length($sequence) ) )
	{    #only one K at the end of sequence
		 #	print "$xids $sequence nr sites only ",$nrxlinksites," at the end at ",$xlinksites[0]," discarded\n";
		return 0;
	}
	elsif (!( $PARAMS->{'nocutatxlink'} )
		&& ( $nrxlinksites == 1 )
		&& ( $xlinksites[0] == length($sequence) ) )
	{    #only one K but in the middle of peptide
		 #	print "$xids $sequence nr sites only ",$nrxlinksites," at the end at ",$xlinksites[0]," accepted\n";
		return 1;
	}
	else {
		warn "problem with $sequence $!";
		return 0;
	}

}

sub getmodified_seq {
	my $subseq  = shift;
	my $PARAMS  = shift;
	my $basicid = shift;
	my $modnterm=shift;
	
	#print "Modnterm: $modnterm";
	
	
	my %seqcombinations;
	my $i;

	if ( $PARAMS->{'variable_mod'} ) {
		my ( $modAA, $shift ) = split /,|:/, $PARAMS->{'variable_mod'};
		my $nvariablemods = $PARAMS->{'nvariable_mod'};

		$seqcombinations{$basicid}->{'seq'} = $subseq;
		my @xlinksites = ();
		while ( $subseq =~ /$modAA/gi ) {    #all possible sites of modification
			push @xlinksites, pos($subseq) - 1;
		}

		my @seq = split //, $subseq;

		#print "\n$subseq\n";
		for $i ( 1 .. $nvariablemods )
		{    #get all possible combinations from 1 to nvariablemods
			my $combinat = Math::Combinatorics->new(
				count => $i,
				data  => [@xlinksites],
			);

		#			print "combinations of $i from: " . join( " ", @xlinksites ) . "\n";
		#			print "------------------------"
		#			  . ( "--" x scalar(@xlinksites) ) . "\n";
			while (
				my @combo =
				sort { $a <=> $b } $combinat->next_combination
			  )
			{

				#		print join( ' ', @combo ) . "\n";
				my @tmpseq = @seq;
				foreach (@combo) {
					$tmpseq[$_] = 'X'
					  ;   #replace K with X, X is defined as AA in the masstable
				}
				my $replacedseq = join "", @tmpseq;

				#				print "\n";
				my $id = join "::", $basicid, "X", @combo;
				$seqcombinations{$id}->{'seq'} = $replacedseq;
			}
		}

		#		print "\n";

	}
	elsif ( ($PARAMS->{'variable_mod'}) && ($modnterm)) {  
		
		## modify the peptide and generate all variable modifications
		my ( $modAA, $shift ) = split /,|:/, $PARAMS->{'variable_mod'};
		my $nvariablemods = $PARAMS->{'nvariable_mod'};

		#$seqcombinations{$basicid}->{'seq'} = $subseq;
		## modify the subsequence add Z which is defined in the mass table as nterm
		my $subseq="Z".$subseq;
		
		my @xlinksites = ();
		while ( $subseq =~ /$modAA/gi ) {    #all possible sites of modification
			push @xlinksites, pos($subseq) - 1;
		}

		my @seq = split //, $subseq;

		#print "\n$subseq\n";
		for $i ( 1 .. $nvariablemods )
		{    #get all possible combinations from 1 to nvariablemods
			my $combinat = Math::Combinatorics->new(
				count => $i,
				data  => [@xlinksites],
			);

		#			print "combinations of $i from: " . join( " ", @xlinksites ) . "\n";
		#			print "------------------------"
		#			  . ( "--" x scalar(@xlinksites) ) . "\n";
			while (
				my @combo =
				sort { $a <=> $b } $combinat->next_combination
			  )
			{

				#		print join( ' ', @combo ) . "\n";
				my @tmpseq = @seq;
				foreach (@combo) {
					$tmpseq[$_] = 'X'
					  ;   #replace K with X, X is defined as AA in the masstable
				}
				my $replacedseq = join "", @tmpseq;

				#				print "\n";
				my $id = join "::", $basicid, "Z", @combo;
				$seqcombinations{$id}->{'seq'} = $replacedseq;
			}
		}
	
	
	
	}else{
	 #do nothing with peptide	
	$seqcombinations{$basicid}->{'seq'} = $subseq;
	
	if ($modnterm){
	## mod nterminus if ntermod param is set
	my @idarray=split(/::/,$basicid);
	my @combo;
	push @combo, $idarray[1];
	push @combo, $idarray[2];
	my $modid=join "::",$basicid,"Z",@combo;
	my $subseq2="Z".$subseq;
	$seqcombinations{$modid}->{'seq'} = $subseq2;
	#print "Modifying nterminus: \n";
	#print "Basiq ID: $basicid, basic Sequecne $subseq\n";
	#print "mod ID: $modid, basic Sequecne $subseq2\n";
	}
	}
	return \%seqcombinations;
}
		

sub _digest_mixed {
	my $seqobj   = shift;
	my $ENZ      = shift;
	my $seq      = $seqobj->seq;
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};
	#print "digest LysN\n";
	#exit;
	while ( $seq =~ /(K|R|E|D)/gi ) {
		push @cutsites, pos($seq) + 1;
		#push @cutsites, pos($seq);
	}
	return @cutsites;
}


sub _digest_LysN {
	my $seqobj   = shift;
	my $ENZ      = shift;
	my $seq      = $seqobj->seq;
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};
	#print "digest LysN\n";
	#exit;
	while ( $seq =~ /(K)/gi ) {
		push @cutsites, pos($seq);
	}
	return @cutsites;
}

sub _digestLysC {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	while ( $seq =~ /(K)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}
		#push @cutsites, pos($seq);
	}
	#print "@cutsites\n";
	return @cutsites;
}

sub _digestGluC {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	while ( $seq =~ /(D|E)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}

sub _definedEnzyme {
	my $seqobj = shift;
	my $PARAMS = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	my @cutsites = ();
	my $cutAA    = $PARAMS->{'cutAA'};
	my $notcutAA = $PARAMS->{'notcutAA'};
	my $cutterm  = $PARAMS->{'cutterm'};

	if ( $cutterm =~ /c/i ) {
		while ( $seq =~ /($cutAA)/gi ) {
			unless ( $seq[ pos($seq) ] =~ /$notcutAA/i ) {
				push @cutsites, pos($seq) + 1;
			}
		}
	}
	elsif ( $cutterm =~ /n/i ) {
		while ( $seq =~ /($cutAA)/gi ) {
			unless ( $notcutAA && ($seq[ pos($seq) ] =~ /$notcutAA/i) ) {
				push @cutsites, pos($seq);
			}
		}
	}
	return @cutsites;
}

sub _digestTryps {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#        print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#        while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(K|R)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}




sub _digestTryps_silac {
    my $seqobj = shift;
    my $ENZ    = shift;
    my $seq    = $seqobj->seq;
    my @seq    = split //, $seq;

        #        print "$seq\n";                                                                                                                                                                           
    my @cutsites = ();
    my $cutAA    = $ENZ->{'cutAA'};
    my $notcutAA = $ENZ->{'notcutAA'};

        #        while ( $seq =~ /(K[^P]|R[^P])/gi ) {                                                                                                                                                     
    while ( $seq =~ /(K|R|B|J|U)/gi ) {
	unless ( $seq[ pos($seq) ] =~ /P/i ) {
	    push @cutsites, pos($seq) + 1;
	}

                #push @cutsites, pos($seq);                                                                                                                                                                
    }

        #print "@cutsites\n";                                                                                                                                                                              
    return @cutsites;
}

sub _digestTryps_highspec {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(K|R)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P|C/i ) {
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}

sub _digestTryps_lowspec {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(K|R)/gi ) {
		push @cutsites, pos($seq) + 1;
	}
	return @cutsites;
}

sub _digestChymoTryps {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(F|Y|W)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}



sub _digest_chymotrypsin3 {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	while ( $seq =~ /(F|L|Y|W)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}

	}

	return @cutsites;
}

sub _digest_chymotrypsin4 {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	while ( $seq =~ /(F|Y|W|L|M|H|C|)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}
	}

	return @cutsites;
}




sub _digestChymoTryps_lowspec {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(F|Y|W|M|L)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}

sub _digestTrypsAspN {
	my $seqobj   = shift;
	my $ENZ      = shift;
	my $seq      = $seqobj->seq;
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};
	while ( $seq =~ /(D|K[^P]|R[^P])/gi ) {
		push @cutsites, pos($seq);
	}
	return @cutsites;
}

sub _digestTrypsR {
	my $seqobj = shift;
	my $ENZ    = shift;
	my $seq    = $seqobj->seq;
	my @seq    = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(R)/gi ) {
		unless ( $seq[ pos($seq) ] =~ /P/i ) {
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}

sub _digestAspN {
	my $seqobj   = shift;
	my $ENZ      = shift;
	my $seq      = $seqobj->seq;
	my @cutsites = ();
	while ( $seq =~ /(D|E)/gi ) {
		push @cutsites, pos($seq);
	}
	return @cutsites;
}

1;
