#chr1    1189291 1209234 UBE2J2;UBE2J2;UBE2J2;UBE2J2
#chr1    1215815 1227409 SCNN1D;SCNN1D
#chr1    1227763 1243269 ACAP3;MIR6726

while(<>)	{
	chomp;
	split(/\t/);
	@gene=split(/;/,@_[3]);
	my @unique = do { my %seen; grep { !$seen{$_}++ } @gene };
	print "@_[0]\t@_[1]\t@_[2]\t" . join (";",@unique) . "\n";
}

