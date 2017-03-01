#!/usr/bin/perl -w

use strict;
use EasyImport::Core;

## load parameters from an INI-style config file
my %sections = (
  'FILES' => 	{	'GFF' => 1
        },
  'GENE_DESCRIPTIONS' =>	{ },
  'GENE_STABLE_IDS' =>	{ },
  'TRANSCRIPT_DESCRIPTIONS' =>	{ },
  'TRANSCRIPT_STABLE_IDS' =>	{ }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}
## connect to core database
#my $dbh = core_db_connect($params);

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
foreach my $subsection (sort keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

## get name/description information from fasta files if defined in the config
## use pattern matching and substitution to allow for non-exact correspondence between names in different files
## preferentially select descriptions from sources with lower indices
my %gene_properties;
foreach my $subsection (keys %{$params->{'GENE_DESCRIPTIONS'}}){
	next if $subsection eq 'GFF';
	get_properties(\%{$infiles{$subsection}},'description',\%gene_properties,$params->{'GENE_STABLE_IDS'}{$subsection},$params->{'GENE_DESCRIPTIONS'}{$subsection});
}
foreach my $subsection (keys %{$params->{'GENE_NAMES'}}){
	next if $subsection eq 'GFF';
	get_properties(\%{$infiles{$subsection}},'synonym',\%gene_properties,$params->{'GENE_STABLE_IDS'}{$subsection},$params->{'GENE_NAMES'}{$subsection});
}
my %transcript_properties;
foreach my $subsection (keys %{$params->{'TRANSCRIPT_DESCRIPTIONS'}}){
	next if $subsection eq 'GFF';
	get_properties(\%{$infiles{$subsection}},'description',\%transcript_properties,$params->{'TRANSCRIPT_STABLE_IDS'}{$subsection},$params->{'TRANSCRIPT_DESCRIPTIONS'}{$subsection});
}
foreach my $subsection (keys %{$params->{'TRANSCRIPT_NAMES'}}){
	next if $subsection eq 'GFF';
	get_properties(\%{$infiles{$subsection}},'synonym',\%transcript_properties,$params->{'TRANSCRIPT_STABLE_IDS'}{$subsection},$params->{'TRANSCRIPT_NAMES'}{$subsection});
}

## check that gff file is properly formatted, repair problems that can be fixed and weave in gene/transcript names/descriptions
## TODO: include test for phase, currently using flag

my $gff_file = rewrite_gff($params,\%infiles,\%gene_properties,$params->{'GENE_STABLE_IDS'}{'GFF'},$params->{'GENE_DESCRIPTIONS'}{'GFF'},$params->{'GENE_NAMES'}{'GFF'},\%transcript_properties,$params->{'TRANSCRIPT_STABLE_IDS'}{'GFF'},$params->{'TRANSCRIPT_DESCRIPTIONS'}{'GFF'},$params->{'TRANSCRIPT_NAMES'}{'GFF'},$params->{'TRANSLATION_STABLE_IDS'}{'GFF'});

die "ERROR: no features in rewritten gff file\n" unless $gff_file;

## summarise the rewritten gff file
gff_feature_summary($params,$gff_file);




sub usage {
	return "USAGE: perl /path/to/prepare_gff.pl /pat/to/config_file.ini";
}
