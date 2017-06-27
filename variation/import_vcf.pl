#!/usr/bin/perl -w

use strict;
use strict;
use File::Basename;

use EasyImport::Core;
use EasyImport::Variation;

use Bio::EnsEMBL::Registry;

$| = 1;

## load parameters from an INI-style config file
my %sections = (
  'DATABASE_CORE' =>	{ 	'NAME' => 1,
              'HOST' => 1,
              'PORT' => 1,
              'RO_USER' => 1
            },
  'DATABASE_VARIATION' =>	{ 	'NAME' => 1
            },
  'META' =>	{ 	'SPECIES.PRODUCTION_NAME' => 1,
              'SPECIES.SCIENTIFIC_NAME' => 1
            },
  'FILES' =>	{ 	'VCF' => 1,
              'PANEL' => 1
            },
  'STUDY' =>	{ 	'SOURCE' => 1
            }

  );
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
	load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

## download/obtain files using methods suggested by file paths and extensions
my %infiles;
foreach my $subsection (sort keys %{$params->{'FILES'}}){
	($infiles{$subsection}{'name'},$infiles{$subsection}{'type'}) = fetch_file($params->{'FILES'}{$subsection});
}

# set database name if not specified
$params->{'DATABASE_VARIATION'}{'NAME'} = $params->{'DATABASE_CORE'}{'NAME'} unless $params->{'DATABASE_VARIATION'}{'NAME'};
$params->{'DATABASE_VARIATION'}{'NAME'} =~ s/_core_/_variation_/;
$params->{'DATABASE_VARIATION'}{'HOST'} = $params->{'DATABASE_CORE'}{'HOST'} unless $params->{'DATABASE_VARIATION'}{'HOST'};
$params->{'DATABASE_VARIATION'}{'RW_USER'} = $params->{'DATABASE_CORE'}{'RW_USER'} unless $params->{'DATABASE_VARIATION'}{'RW_USER'};
$params->{'DATABASE_VARIATION'}{'RW_PASS'} = $params->{'DATABASE_CORE'}{'RW_PASS'} unless $params->{'DATABASE_VARIATION'}{'RW_PASS'};
$params->{'DATABASE_VARIATION'}{'PORT'} = $params->{'DATABASE_CORE'}{'PORT'} unless $params->{'DATABASE_VARIATION'}{'PORT'};

# print a registry file
open REG,">".$params->{'DATABASE_CORE'}{'NAME'}.".registry.conf";
print REG <<EOT;
use strict;

use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host    => '$params->{'DATABASE_CORE'}{'HOST'}',
  -user    => '$params->{'DATABASE_CORE'}{'RO_USER'}',
  -port    => '$params->{'DATABASE_CORE'}{'PORT'}',
  -species => '$params->{'META'}{'SPECIES.PRODUCTION_NAME'}',
  -group   => 'core',
  -dbname  => '$params->{'DATABASE_CORE'}{'NAME'}'
);

new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
  -host    => '$params->{'DATABASE_VARIATION'}{'HOST'}',
  -user    => '$params->{'DATABASE_VARIATION'}{'RW_USER'}',
  -pass    => '$params->{'DATABASE_VARIATION'}{'RW_PASS'}',
  -port    => '$params->{'DATABASE_VARIATION'}{'PORT'}',
  -species => '$params->{'META'}{'SPECIES.PRODUCTION_NAME'}',
  -group   => 'variation',
  -dbname  => '$params->{'DATABASE_VARIATION'}{'NAME'}'
);

my \@aliases = ( \'$params->{'META'}{'SPECIES.SCIENTIFIC_NAME'}\' );

Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
  -species => '$params->{'META'}{'SPECIES.PRODUCTION_NAME'}',
  -alias   => \\\@aliases
);

1;
EOT
close REG;

# make lists of seq_regions with and without transcripts
#load $params->{'DATABASE_CORE'}{'NAME'}.'.registry.conf';
#system "cat ".$params->{'DATABASE_CORE'}{'NAME'}.'.registry.conf';
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($params->{'DATABASE_CORE'}{'NAME'}.'.registry.conf');
my $slice_adaptor = $reg->get_adaptor( $params->{'META'}{'SPECIES.PRODUCTION_NAME'}, 'core', 'Slice' );
my $tr_adaptor    = $reg->get_adaptor( $params->{'META'}{'SPECIES.PRODUCTION_NAME'}, 'core', 'Transcript' );
my @supercontigs  = @{$slice_adaptor->fetch_all('toplevel')};
open WITH, ">WITH.REGIONS";
open WITHOUT, ">WITHOUT.REGIONS";
my (%with,%without);
foreach my $slice (@supercontigs) {
  my @transcripts = @{$tr_adaptor->fetch_all_by_Slice($slice)};
  my $row = $slice->seq_region_name()."\t".$slice->start()."\t".$slice->end()."\n";
  if (scalar(@transcripts) > 0){
    print WITH $row;
  }
  else {
    print WITHOUT $row;
  }
}
close WITH;
close WITHOUT;

# extract samples to keep from panel file
my $panfile = $infiles{'PANEL'}{'name'};
system "cut -f1 $panfile > SAMPLES";

# split the input vcf into with and without files
my $infile = $infiles{'VCF'}{'name'};
warn "indexing vcf file\n";
system 'tabix -f '.$infiles{'VCF'}{'name'};

my $withfile = $infile;
$withfile =~ s/\.vcf\.gz/.with.vcf.gz/;
my $withoutfile = $infile;
$withoutfile =~ s/\.vcf\.gz/.without.vcf.gz/;
my $filter = '';
if ($params->{'VCF'} && $params->{'VCF'}{'FILTER'}){
  $filter = "-i ".$params->{'VCF'}{'FILTER'};
}

system "bcftools view $filter -S SAMPLES -R WITH.REGIONS $infile | bgzip -c > $withfile";
system "bcftools view $filter -S SAMPLES -R WITHOUT.REGIONS $infile | bgzip -c > $withoutfile";

# index the vcf files
system "tabix -f $withfile";
system "tabix -f $withoutfile";

# run the variant effect predictor
if (! -d 'vep'){
  system "perl /ensembl/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl"
         .' --input_file $withfile'
         .' --skip_db_check'
         .' --species "'.$params->{'META'}{'SPECIES.PRODUCTION_NAME'}.'"'
         .' --build all'
         .' --registry "'.$params->{'DATABASE_CORE'}{'NAME'}.'.registry.conf"'
         .' --host "'.$params->{'DATABASE_CORE'}{'HOST'}.'"'
         .' --port "'.$params->{'DATABASE_CORE'}{'PORT'}.'"'
         .' --user "'.$params->{'DATABASE_CORE'}{'RO_USER'}.'"'
         .' --dir "vep"'
         .' --tsl'
         .' --protein'
         ;
  system "mv vep/".$params->{'META'}{'SPECIES.PRODUCTION_NAME'}."/*/* vep/";
  system "rm -rf vep/".$params->{'META'}{'SPECIES.PRODUCTION_NAME'};
}

# create the variation database from a template
setup_variation_db($params);

# loop through vcf files and add to database
import_chunk($withoutfile,$params);
import_chunk($withfile,$params,1);

sub usage {
	return "USAGE: perl import_vcf.pl ini_file";
}

sub import_chunk {
 my ($file,$params,$use_vep) = @_;
 my $cmd = "perl /ensembl/ensembl-variation/scripts/import/import_vcf.pl"
       .' --input_file "'.$file.'"'
       .' --tmpdir "/tmp"'
       .' --species "'.$params->{'META'}{'SPECIES.PRODUCTION_NAME'}.'"'
       .' --source "'.$params->{'STUDY'}{'SOURCE'}.'"'
       .' --panel "'.$panfile.'"'
       .' --registry "'.$params->{'DATABASE_CORE'}{'NAME'}.'.registry.conf"'
       .' --sql "/ensembl/ensembl-variation/sql/table.sql"'
       .' --coord_system "scaffold"'
       .' --var_prefix "var"'
       ;

#  if ($params->{'STUDY'}{'SOURCE_DESCRIPTION'}){
#    $cmd .= ' --source_description "'.$params->{'STUDY'}{'SOURCE_DESCRIPTION'}.'"'
#  }
  if ($use_vep){
    $cmd .= ' --add_tables transcript_variation';
    $cmd .= ' --cache "vep"';
  }
  system $cmd;

}


sub fetch_file {
	my ($type,$location,$new_name);
	$location = shift;
	if (ref $location){
		$type = $location->[0];
		if (defined $location->[2]){
			$new_name = $location->[2];
		}
		$location = $location->[1];
	}
	# work out file name from location
	$location =~ m/.+\/([^\/]+)$/;
	my $filename = $1 ? $1 : $location;
	my $command;
	my $compression = '';
	if ($filename =~ s/\.(gz|gzip|tar\.gz|tgz|zip)$//){
		$compression = ".".$1;
	}
	if ($type eq 'vcf'){
		$filename = $filename.$compression;
		$compression = '';
	}
	if (($new_name && !-e $new_name) || (!$new_name && !-e $filename)){
		if ($location =~ m/^(?:ftp|http|https):/){
			$command = 'wget';
			system "wget \"$location\" -O $filename"."$compression";
		}
		elsif ($location =~ m/:[\/~]/){
			$command = 'scp';
			system "scp $location $filename"."$compression";
		}
		else {
			$command = 'cp';
			system "cp $location $filename"."$compression";
		}
	}
	if ($compression && !-e $filename){
		if (!-e $filename.$compression){
			# something went wrong when copying the file
			die "ERROR: could not $command $location to $filename"."$compression\n";
		}
		else {
			if ($compression =~ m/^\.t/){
				system "tar xf $filename"."$compression $new_name";
			}
			elsif ($compression =~ m/^\.g/){
				system "gunzip $filename"."$compression";
			}
			elsif ($compression =~ m/^\.z/){
				system "unzip $filename"."$compression";
			}
			if (!-e $filename && !-e $new_name){
				# this compression type is not currently supported
				die "ERROR: could not extract $filename"."$compression to $filename\n";
			}
		}
	}
	if (!-e $filename){
		# something went wrong when copying the file
		die "ERROR: could not $command $location to $filename\n";
	}
	if (!$type){
		# TODO: infer type properly
		if ($filename =~ m/\.([^\.])$/){
			$type = $1;
			warn "WARNING: no type specified for $filename from $location, inferring '$type' based on file extension\n";

		}
		else {
			die "ERROR: no type specified for $filename from $location, unable to infer based on file extension\n";
		}
	}
	if (defined $new_name){
		system "mv $filename $new_name";
		if (!-e $new_name){
			die "ERROR: could not move $filename from $location to $new_name\n";
		}
		$filename = $new_name;
	}
	if ($type =~ m/^(?:fa|faa|fas|fasta|fna|fsa)/i  ){
		$type = 'fas';
	}
	elsif ($type =~ m/^gff/i){
		$type = 'gff';
	}
	elsif ($type =~ m/^agp/i){
		$type = 'agp';
	}
	elsif ($type =~ m/^csv/i){
		$type = 'csv';
	}
	elsif ($type =~ m/^tsv/i){
		$type = 'tsv';
	}
	return ($filename,$type);
}
