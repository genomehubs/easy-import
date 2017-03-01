#!/usr/bin/perl -w

use strict;
use EasyImport::Core;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

## load parameters from an INI-style config file
my %sections = (
  'DATABASE_CORE' =>    {
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1,
              'RO_USER' => 1
            }
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
$params->{'DATABASE_CORE'}{'NAME'} = shift @ARGV;
while (my $ini_file = shift @ARGV){
        load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => $params->{'DATABASE_CORE'}{'HOST'},
    -user => $params->{'DATABASE_CORE'}{'RO_USER'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql',
);

my $dbname = $params->{'DATABASE_CORE'}{'NAME'};

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $params->{'DATABASE_CORE'}{'RW_USER'},
    -pass   => $params->{'DATABASE_CORE'}{'RW_PASS'},
    -dbname => $dbname,
    -host   => $params->{'DATABASE_CORE'}{'HOST'},
    -port   => $params->{'DATABASE_CORE'}{'PORT'},
    -driver => 'mysql'
);


my $meta_container = $dba->get_adaptor("MetaContainer");

my $display_name    = $meta_container->get_display_name();
my $assembly_name   = $meta_container->single_value_by_key('ASSEMBLY.NAME');

open INI,">/import/conf/".$params->{'DATABASE_CORE'}{'NAME'}.".ini";
print INI <<EOF;
[DATABASE_CORE]
  NAME = $params->{'DATABASE_CORE'}{'NAME'}
[META]
  SPECIES.DISPLAY_NAME = $display_name
  ASSEMBLY.DEFAULT = $assembly_name
EOF
close INI;
