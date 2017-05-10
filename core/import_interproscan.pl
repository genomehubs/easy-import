#!/usr/bin/perl -w

use strict;
use EasyImport::Core;
use EasyImport::Xref;

## load parameters from an INI-style config file
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file);
}

## download/obtain files using methods suggested by file paths and extensions
die "ERROR: an interproscan result file must be specified at [FILES] IPRSCAN\n" unless $params->{'FILES'}{'IPRSCAN'};
my %infiles;
($infiles{'IPRSCAN'}{'name'},$infiles{'IPRSCAN'}{'type'}) = fetch_file($params->{'FILES'}{'IPRSCAN'});

## check that an IPRSCAN file has been specified
die "ERROR: an interproscan result file must be specified at [FILES] IPRSCAN\n" unless $infiles{'IPRSCAN'};

## connect to the db
my $dbh = core_db_connect($params);

## load the IPRSCAN file into the database
push @ARGV,$infiles{'IPRSCAN'}{'name'};
read_iprscan($dbh,$params);

