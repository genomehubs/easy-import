#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Find;
use Module::Load;
use File::Copy;

## find the full path to the directory that this script is executing in
our $dirname;
BEGIN {
  $dirname  = dirname(abs_path($0));
}
use lib "$dirname/../modules";
use lib "$dirname/../gff-parser";
use EasyImport::Core;
use EasyImport::Compara;

## load parameters from an INI-style config file
my %sections = (
    'TAXA' => {},
    'SETUP' => { 
                'FASTA_DIR' => 1,
                'TMP_DIR' => 1,
                'ORTHOFINDER_DIR' => 1,
#                'GENETIC_CODE' => 1,
               },
    'ORTHOGROUP' => {
                     'PREFIX' => 1,
  	             'PROTEIN' => 1,
  	             'PROTEIN_ALIGN' => 1,
  	             'PROTEIN_TRIMMED' => 1,
  	             'FNAFILE' => 1,
  	             'BOUNDEDFILE' => 1,
  	             'TREE' => 1,
                     'HOMOLOG' => 1,
                    }
);
## check that all required parameters have been defined in the config file
die "ERROR: you must specify at least one ini file\n",usage(),"\n" unless $ARGV[0];
my %params;
my $params = \%params;
while (my $ini_file = shift @ARGV){
        load_ini($params,$ini_file,\%sections,scalar(@ARGV));
}

my %sequences;


my $orthologousgroupstxt_filename = $params{SETUP}->{ORTHOFINDER_DIR}.'/Orthogroups/Orthogroups.txt';

my @species_list = keys %{$params{TAXA}};
my $species_regexp = join("|", @species_list);

my $fastadir = $params{SETUP}->{FASTA_DIR};
my $aligndir = $params{SETUP}->{TMP_DIR};
my $trimdir = $params{SETUP}->{ORTHOFINDER_DIR}.'/MultipleSequenceAlignments';
my $treedir = $params{SETUP}->{ORTHOFINDER_DIR}.'/Resolved_Gene_Trees/';
my $seqidfile = $params{SETUP}->{ORTHOFINDER_DIR}.'/WorkingDirectory/SequenceIDs.txt';
my $dupfile = $params{SETUP}->{ORTHOFINDER_DIR}.'/Gene_Duplication_Events/Duplications.tsv';
my $ortholog_dir = $params{SETUP}->{ORTHOFINDER_DIR}.'/Orthologues';
my $orthogroup_dir = $params{SETUP}->{ORTHOGROUP_DIR};
my $genetree_suff = $params{ORTHOGROUP}->{TREE};                
my $align_suff = $params{ORTHOGROUP}->{PROTEIN_ALIGN}[0];                
my $trimmed_suff = $params{ORTHOGROUP}->{PROTEIN_TRIMMED};                
my $homolog_suff = $params{ORTHOGROUP}->{HOMOLOG};                
my $species_tree = $params{SETUP}->{ORTHOFINDER_DIR}.'/Species_Tree/SpeciesTree_rooted_node_labels.txt';                
#my $genetic_code = textfile2hash($params{SETUP}->{GENETIC_CODE});

# parse species tree
open TREE,"<$species_tree";
my $nwk;
while (<TREE>){
  s/\r\n/\n/;
  chomp;
  $nwk .= $_;
}
close TREE;
my $st_parsed = parse_tree($nwk);
my $st_nodes = nodes_by_taxa($st_parsed);


# load all sequences in memory
for my $species (@species_list) {
  $sequences{$species}{faa} = fastafile2hash("$fastadir/canonical_proteins/$species.fa", );
  $sequences{$species}{fna} = fastafile2hash("$fastadir/canonical_cds_translationid/$species.fa");
  $sequences{$species}{fba} = fastafile2hash("$fastadir/canonical_protein_bounded_exon/$species.fa");
}

# load sequence IDs in memory
open ID, "<$seqidfile" or die $!;
my %ids;
while (<ID>){
  s/\r\n/\n/;
  chomp;
  my ($id,$title) = split /[:\s]+/;
  $ids{$id} = $title;
}
close ID;

# load gene duplications
open DUP, "<$dupfile" or die $!;
my %dups;
while (<DUP>){
  s/\r\n/\n/;
  chomp;
  my ($og,$stnode,$gtnode,$support) = split /\s+/;
  if (!$dups{$og}) {
    $dups{$og} = {}
  }
  $dups{$og}{$gtnode} = $support;
}
close DUP;

# load orthologues
my %orthos;
for my $spA (@species_list){
  for my $spB (@species_list){
    my $file = "$ortholog_dir/Orthologues_$spA/$spA\__v__$spB.tsv";
    if (-e $file){
      open ORTH,"<$file" or die $!;
      <ORTH>;
      while (<ORTH>){
        s/\r\n/\n/;
        chomp;
        my ($og,$alist,$blist) = split /\t/;
        my @agenes = split /,\s+/,$alist;
        my @bgenes = split /,\s+/,$blist;
        foreach my $agene (@agenes){
          foreach my $bgene (@bgenes){
            $orthos{$og}{$agene}{$bgene} = 1; 
            $orthos{$og}{$bgene}{$agene} = 1; 
          }
        }
      }
      close ORTH
    }
  }
}


# create folder with sequences for each orthogroup
open OG, "<$orthologousgroupstxt_filename" or die $!;
mkdir "orthogroups";
<OG>;
while (<OG>) {
  s/\r\n/\n/;
  chomp;
  my @tokens = split /\s+/;
  my $orthogroup_id = shift @tokens;
  next if scalar @tokens < 2;
  $orthogroup_id =~ s/:$//;
  print $orthogroup_id,"\n";
  my $genetree_id = $orthogroup_id;
  my $prefix = $params{ORTHOGROUP}->{PREFIX} || 'OG';
  $genetree_id =~ s/^OG/$prefix/;
  my $orthodir = "$orthogroup_dir/$genetree_id";
  if (-e "$treedir/$orthogroup_id\_tree.txt"){
  mkdir $orthodir;
  open IDS, ">$orthodir/$genetree_id" or die $!;
  open FAA, ">$orthodir/$genetree_id.faa" or die $!;
  open FNA, ">$orthodir/$genetree_id.fna" or die $!;
  open FBA, ">$orthodir/$genetree_id.fba" or die $!;
  my @gene_ids;
  my %trimmed_seqs = %{fastafile2hash("$trimdir/$orthogroup_id.fa")};
  for my $full_sequence_id (@tokens) {
    if (exists $trimmed_seqs{$full_sequence_id}){
      if ($full_sequence_id =~ /^(${species_regexp})_(\S+)$/) {
        my $species = $1;
        my $sequence_id = $2;
        if (exists $sequences{$species}{faa}{$full_sequence_id} and 
            exists $sequences{$species}{fna}{$full_sequence_id} and
            exists $sequences{$species}{fba}{$full_sequence_id}){
          push @gene_ids, $full_sequence_id;
          print IDS "$full_sequence_id\n";
          print FAA ">$full_sequence_id\n" . $sequences{$species}{faa}{$full_sequence_id}->{seq} . "\n";
          print FNA ">$full_sequence_id\n" . $sequences{$species}{fna}{$full_sequence_id}->{seq} . "\n";
          print FBA ">$full_sequence_id\n" . $sequences{$species}{fba}{$full_sequence_id}->{seq} . "\n";
        }
        else {
          die "$full_sequence_id does not exist in one of the input fasta files";
        }
      }
      else {
        warn "$full_sequence_id does not match species list";
      }
    }
  }
  close IDS;
  close FAA;
  close FNA;
  close FBA;

 #   my %protalign = %{fastafile2hash("$aligndir/$orthogroup_id.fa_pre_trim" )};
 #   my %nucalign;
 #   foreach my $idA (keys %protalign){
 #     my $nameA = $ids{$idA};
 #     $nameA =~ /^(${species_regexp})_(\S+)$/;
 #     my $speciesA = $1;
 #     $nucalign{$nameA} = prot2nuc($protalign{$idA}->{seq},$sequences{$speciesA}{fna}{$nameA}->{seq});
 #     foreach my $idB (keys %protalign){
 #       next if $idA eq $idB;
 #       my $nameB = $ids{$idB};
 #       $nameB =~ /^(${species_regexp})_(\S+)$/;
 #       unless ($nucalign{$nameB}){
 #         my $speciesB = $1;
 #         $nucalign{$nameB} = prot2nuc($protalign{$idB}->{seq},$sequences{$speciesB}{fna}{$nameB}->{seq});
 #       }
 #       my ($qid,$tid) = pairwise_identity($protalign{$idA}->{seq},$protalign{$idB}->{seq});
 #       my ($dn,$ds) = dn_ds($nucalign{$nameA},$nucalign{$nameB});
 #       print $protalign{$idA}->{seq},"\n";
 #       print $protalign{$idB}->{seq},"\n";
 #       print $qid," - ",$tid,"\n";
 #       #print $dn," : ",$ds,"\n";
 #     }
 #   }

    open ALIN,"<$aligndir/$orthogroup_id.fa_pre_trim" or die $!;
    open ALOUT, ">$orthodir/$genetree_id$align_suff" or die $!;
    while (<ALIN>){
      s/^>(\S+)/>$ids{$1}/;
      print ALOUT;
    }
    close ALIN;
    close ALOUT;

    open TRIN,"<$trimdir/$orthogroup_id.fa" or die $!;
    open TROUT, ">$orthodir/$genetree_id$trimmed_suff" or die $!;
    while (<TRIN>){
      s/(\w+?_)\w+?_/$1/g;
      print TROUT;
    }
    close TRIN;
    close TROUT;

    open GTIN,"<$treedir/$orthogroup_id\_tree.txt" or die $!;
    open GTOUT, ">$orthodir/$genetree_id$genetree_suff" or die $!;
    my $gt_nwk;
    while (<GTIN>){
      s/\r\n/\n/;
      chomp;
      $gt_nwk .= $_;
    }
    close GTIN;
    $gt_nwk =~ s/[A-Z]{6}_([A-Z]{6}_)/$1/g;
    my $parsed = parse_tree($gt_nwk);
    my $gt_nodes = gt_nodes_by_taxa($parsed,$st_nodes,$st_parsed);
        foreach my $gtnode (keys %{$gt_nodes}){
          my $dup = 'N';
          my $spnode = $gt_nodes->{$gtnode};
          my $score = 0;
          if ($dups{$orthogroup_id} && $dups{$orthogroup_id}{$gtnode}){
            $dup = 'Y';
            $score = $dups{$orthogroup_id}{$gtnode};
          }
          if ($gtnode eq 'n0'){
            $gt_nwk =~ s/(\b$gtnode\b)/$1\[&&NHX:S=$spnode:D=$dup:Bl=0:L=$score\]/;
          }
          else {
            $gt_nwk =~ s/(\b$gtnode\b[:\),]*([0-9\.e-]*))/$1\[&&NHX:S=$spnode:D=$dup:Bl=$2:L=$score\]/;
          }
        } 
      print GTOUT $gt_nwk;
    close GTOUT;

    open HTOUT, ">$orthodir/$genetree_id$homolog_suff" or die $!;
    print HTOUT <<EOF;
Homolog Table for: $genetree_id$genetree_suff
P == Paralogous
O == Orthologous
X == Xenologous
. == Genes on X and Y axis are the same.

EOF
    for my $gene_id (@gene_ids){
      print HTOUT "\t$gene_id";
    }
    print HTOUT "\n";
    for my $geneA (@gene_ids){
      print HTOUT "$geneA";
      for my $geneB (@gene_ids){
        if ($geneA eq $geneB){
          print HTOUT "\t.";
        }
        elsif ($orthos{$orthogroup_id}{$geneB}{$geneA}){
          print HTOUT "\tO";
        }
        else {
          print HTOUT "\tP";
        }
      }
      print HTOUT "\n";
    }
    close HTOUT;
  }
}


#############################################################################
sub fastafile2hash
{
  my $fastafile  = shift @_;
  my $changecase = "N";
  my $order      = "S"; # S = same as input, or R = random
  $changecase    = substr(uc(shift @_),0,1) if @_;
  $order         = substr(uc(shift @_),0,1) if @_;
  my %sequences;
  my $fh = &read_fh($fastafile);
  my $seqid;
  my $seq_counter;
  while (<$fh>)
  {
    if (/^>(\S+)(.*)/) {
      $seqid = $1;
      $sequences{$seqid}{desc} = $2;
      $sequences{$seqid}{order} = $order eq "S" ? $seq_counter++ : rand;
    }
    else {
      if (/\d/) {
        chomp($sequences{$seqid}{seq} .= " $_"); # add space to sep qual values
        $sequences{$seqid}{seq} =~ s/^\s+//;
        $sequences{$seqid}{seq} =~ s/\s+$//;
        next;
      }
      chomp($sequences{$seqid}{seq} .= lc($_)) if $changecase eq "L";
      chomp($sequences{$seqid}{seq} .= uc($_)) if $changecase eq "U";
      chomp($sequences{$seqid}{seq} .= $_    ) if $changecase eq "N";
    }
  }
  return \%sequences;
}

#############################################################################

sub read_fh {
  my $filename = shift @_;
  my $filehandle;
  if ($filename =~ /gz$/) {
    open $filehandle, "gunzip -dc $filename |" or die "Problem opening $filename\n$!";
  }
  else {
    open $filehandle, "<$filename" or die "Problem opening $filename\n$!";
  }
  return $filehandle;
}

#############################################################################

sub parse_tree {
  my $nwk = shift;
  print $nwk,"\n";
  my %nodes;
  while ($nwk =~ s/\(([^\(\)]+)\)(\w+)/$2/){
    my $match = $1;
    my $node = $2;
    $nodes{$node} = ();
    while ($match =~ s/\b([A-Zn][^:,\(\)]+)://){
      if ($nodes{$1}){
        push @{$nodes{$node}},@{$nodes{$1}};
      }
      else {
        push @{$nodes{$node}},$1;
      }
    }
  }
  foreach my $node (keys %nodes){
    @{$nodes{$node}} = sort @{$nodes{$node}};
  }
  return \%nodes;
}


sub nodes_by_taxa {
  my $nodes = shift;
  my %taxa;
  foreach my $node (keys %{$nodes}){
    $taxa{join '_',sort @{$nodes->{$node}}} = $node;
  }
  return \%taxa;
}

sub unique_gt_taxa {
  my %arr = map { local $_ = $_; s/^([A-Z]+?)_.+/$1/; $_ => 1 } @_;
  return sort keys %arr;
}

sub gt_nodes_by_taxa {
  my $gt_nodes = shift;
  my %nodes;
  my $st_nodes = shift;
  my $st_parsed = shift;
  foreach my $gt_node (keys %{$gt_nodes}){
    my @unique_taxa = unique_gt_taxa(@{$gt_nodes->{$gt_node}});
    if (scalar @unique_taxa > 1){
      my $taxlist = join '_', @unique_taxa;
      if ($st_nodes->{$taxlist}){
        $nodes{$gt_node} = $st_nodes->{$taxlist};
      }
      else {
        my $shortest = "inf";
        foreach my $st_node (keys %{$st_parsed}){
          my @st_taxa = @{$st_parsed->{$st_node}};
          my $st_len = scalar @st_taxa;
          if ($st_len < $shortest){
            # are all elements in @unique_taxa also in @st_taxa?
            my $ctr = 0;
            for (my $i = 0; $i < scalar @unique_taxa; $i++){
              while (my $taxon = shift @st_taxa){
                if ($unique_taxa[$i] eq $taxon){
                  $ctr++;
                  last;
                }
              }
            }
            if ($ctr == scalar @unique_taxa){
              $shortest = $st_len;
              $nodes{$gt_node} = $st_node;
            }
          }
	}
        #$taxlist = join '\w*?',sort @unique_taxa;
      }
      #$nodes{$gt_node} = $st_nodes->{join '_',sort @unique_taxa};
    }
    else {
      $nodes{$gt_node} = $unique_taxa[0];
    }
    foreach my $gene (@{$gt_nodes->{$gt_node}}){
      $gene =~ m/^([A-Z]+)_/;
      $nodes{$gene} = $1;
    }
  }
  return \%nodes;
}

#############################################################################

sub prot2nuc {
  my $prot = shift;
  my $nuc = shift;
  my @prot = split //,$prot;
  my @nuc = $nuc =~ /(.{1,3})/g;
  my $nucaln = '';
  while (my $aa = shift @prot){
    if ($aa eq '-'){
      $nucaln .= '---';
    }
    else {
      $nucaln .= shift @nuc;
    }
  }
  return $nucaln;
}

#############################################################################

sub textfile2hash {
  my $filename = shift;
  my %hash;
  open TXT,"<$filename" or die $!;
  while (<TXT>){
    chomp;
    my ($key, $value) = split /\s+/;
    if ($key && $value){
      $hash{$key} = $value;
    }
  }
  return \%hash;
}

#############################################################################

sub pairwise_identity {
  my $seqA = shift;
  my $seqB = shift;
  my @seqA = split //,$seqA;
  my @seqB = split //,$seqB;
  my ($idA, $idB, $lenA, $lenB);
  my $id = 0;
  for (my $i = 0; $i < scalar @seqA; $i++){
    if ($seqA[$i] ne '-'){
      $lenA++;
      if ($seqA[$i] eq $seqB[$i]){
        $id++;
      }
    }
    if ($seqB[$i] ne '-'){
      $lenB++;
    }
  }
  $idA = $id / $lenA * 100;
  $idB = $id / $lenB * 100;
  return ($idA, $idB);
}

#############################################################################

#sub dn_ds {
#  my $seqA = shift;
#  my $seqB = shift;
#  my $in = "\nquery $seqA\ntarget $seqB\n";
#  my ($out, $err);  
#  run3 ['perl', '/Snap/SNAP.pl', '-plot', 'N', '-list', 'N', '-'], \$in, \$out, \$err;
#  print $out,"\n";
#  warn $err if $err;
#  my ($dn, $ds);
#  return ($dn, $ds);
#}

#############################################################################

sub usage {
  return "USAGE: perl /path/to/make_orthogroup_files.pl ini_file";
}
