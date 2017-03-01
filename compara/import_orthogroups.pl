#!/usr/bin/perl -w

use strict;
use File::Find;
use EasyImport::Core;
use EasyImport::Compara;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Graph::NewickParser;
use Bio::EnsEMBL::Compara::GeneTreeNode;
use Bio::EnsEMBL::Compara::GeneTreeMember;
use Bio::EnsEMBL::Compara::SpeciesTreeNode;
use Bio::EnsEMBL::Compara::SpeciesTree;
use Bio::EnsEMBL::Compara::DBSQL::NCBITaxonAdaptor;

## load parameters from an INI-style config file
my %sections = (
  'DATABASE_COMPARA' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            },
    'TAXA' => {},
    'ORTHOGROUP' => {  'PREFIX' => 1,
  	           'PROTEIN' => 1,
  	           'PROTEIN_ALIGN' => 1,
  	           'PROTEIN_TRIMMED' => 1,
  	           'FNAFILE' => 1,
  	           'BOUNDEDFILE' => 1,
  	           'TREE' => 1
            }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections);
}

my $dbh = compara_db_connect($params);

our $prefix = $params->{'ORTHOGROUP'}{'PREFIX'};
# this is the current version prefix, trim the version digits at the end
$prefix =~ s/\d*$//;

my $taxlist;
foreach my $key (keys %{$params->{'TAXA'}}){
  $taxlist .= $key.'|';
}
chop $taxlist;
$params->{'ORTHOGROUP'}{'TAXA'} = $taxlist;

my ($st_nodes,$core_dbs) = fetch_species_tree_nodes($params,$dbh);

find({wanted => sub {
  my $file = $File::Find::name;
  if (!-d $file){
    if ($file =~ m/$prefix\w+$/){
      warn "importing $file\n";
      load_sequences($dbh,$params,$st_nodes,$core_dbs,$file);
    }
  }
},
no_chdir => 1},'.');

sub usage {
        return "USAGE: perl -I /path/to/dir/containing/Ensembl_Import.pm /path/to/import_homologues.pl ini_file";
}
