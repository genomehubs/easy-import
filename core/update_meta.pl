#!/usr/bin/perl -w

use strict;
use EasyImport::Core;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

## load parameters from an INI-style config file
my %sections = (
  'DATABASE_CORE' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RW_USER' => 1,
              'RW_PASS' => 1
            },
  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

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

# Get/set meta values
my %meta = %{$params->{'META'}};
my %key_exists;
foreach my $key (keys %meta){
  my @values = ref($meta{$key}) eq 'ARRAY' ? @{$meta{$key}} : ($meta{$key});
  $key = lc $key;
  while (my $value = shift @values){
    unless ($meta_container->key_value_exists($key,$value)){
      my @current = @{$meta_container->list_value_by_key($key)};
      if (!@current || @current > 1 || $key_exists{$key}){
        $meta_container->store_key_value($key,$value);
      }
      elsif ($value =~ m/\S/){
        $meta_container->update_key_value($key,$value);
      }
      else {
        $meta_container->delete_key($key);
      }
    }
  }
  $key_exists{$key}++;
}

sub usage {
	return "USAGE: perl /path/to/update_meta.pl /path/to/config_file.ini";
}

__END__
