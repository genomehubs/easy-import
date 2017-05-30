#!/usr/bin/perl -w

use strict;
use EasyImport::Core;
use EasyImport::Compara;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Graph::NewickParser;
use Bio::EnsEMBL::Compara::SpeciesTreeNode;
use Bio::EnsEMBL::Compara::SpeciesTree;

## load parameters from an INI-style config file
my %sections = (
  'DATABASE_COMPARA' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            },
  'DATABASE_TEMPLATE' =>	{ 	'NAME' => 1,
              'URL' => 1
            },
  'SPECIES_SET' => { 'SPECIES_TREE' => 1,
              'TREE_LABEL' => 1
            }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections);
}

#-------------------------------------------------

# create the compara database from a template
# populate ncbi_taxa_node and ncbi_taxa_tree tables

setup_compara_db($params);

add_species_tree($params);

sub usage {
	return "USAGE: perl setup_database.pl ini_file";
}
