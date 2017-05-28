#!/usr/bin/perl -w

use strict;
use File::Basename;
use EasyImport::Core;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

## load parameters from an INI-style config file
my %sections = (
  'ORTHOGROUP'  => {
    'PREFIX'       => 1
  }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

#===============================================================================

system "parallel ORTHOGROUPID={/} /ensembl/easy-import/compara/genetree.sh ::: orthogroups/".$params->{'ORTHOGROUP'}{'PREFIX'}."*"

sub usage {
  return "USAGE: perl /path/to/run_genetrees.pl /path/to/config_file.ini";
}
