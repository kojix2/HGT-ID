#!/usr/bin/perl

print "\nChecking Java version...\n\n";

my $ret = `java -version 2>&1`;
print "$ret\n";

if ( index( $ret, '1.7' ) == -1 ) {
    printf "HGT-ID require Java 1.7 or more to run properly\n\n";
    die "";
}

print "\nChecking SAMtools...\n\n";

$ret = `which samtools 2>&1`;
if ( index( $ret, 'no samtools' ) == -1 ) {
    printf "%-30s\tOK\n\n", 'SAMtools';
}
else {
    printf "%-30s\tnot found\n\n", 'SAMtools';
    die "";
}

print "\nChecking BEDtools...\n\n";

$ret = `which samtools 2>&1`;
if ( index( $ret, 'no bedtools' ) == -1 ) {
    printf "%-30s\tOK\n\n", 'BEDtools';
}
else {
    die "%-30s\tnot found\n\n", 'BEDtools';
}

my @required_modules = (
    "Getopt::Long", "Cwd",        "Data::Dumper", "File::Basename",
    "File::Copy",   "File::Path", "File::Spec",   "File::Temp",
    "FindBin",      "Getopt::Long",
);

print "\nChecking CPAN modules required by HGT-ID...\n\n";
my $count = 0;
for my $module (@required_modules) {

    eval("use $module");
    if ($@) {
        printf "%-30s\tFailed\n", $module;
        $count++;
    }
    else {
        printf "%-30s\tOK\n", $module;
    }
}

if ( $count == 1 ) {
    print "\n\nOne module may not be installed properly.\n\n";
}
elsif ( $count > 1 ) {
    print "\n\n$count modules may not be installed properly.\n\n";
}
else {
    print "\n\nAll CPAN modules checked!\n\n";
}
