#!/usr/bin/perl -w

use strict;
use DBI;
use Bio::SeqIO;
use File::Basename;
use Data::Dumper;

# fill in/select from
sub load_sequences {
  my ($dbh,$params,$st_nodes,$core_dbs,$fullname) = @_;
  my ($file,$path) = fileparse($fullname);
  my %seqs;

  return if not -f $path.'/'.$file.''.$params->{'ORTHOGROUP'}{'HOMOLOG'};

  # Leave out gene_tree_root function for now, assume supertree of treetype clusterset default exists in database with root_id 1
  # my $root_id = add_gene_tree_root($dbh,'protein','clusterset','default',$mlss_id);

  read_seqs_to_hash(\%seqs,$path.'/'.$file.''.$params->{'ORTHOGROUP'}{'PROTEIN'},$params->{'ORTHOGROUP'}{'TAXA'},'protein');
  read_seqs_to_hash(\%seqs,$path.'/'.$file.''.$params->{'ORTHOGROUP'}{'FNAFILE'},$params->{'ORTHOGROUP'}{'TAXA'},'fna');
  read_seqs_to_hash(\%seqs,$path.'/'.$file.''.$params->{'ORTHOGROUP'}{'BOUNDEDFILE'},$params->{'ORTHOGROUP'}{'TAXA'},'bounded');
  my $cluster_id = $file;
  my $aln_method = $params->{'ORTHOGROUP'}{'PROTEIN_ALIGN'}->[1];
  my $aln_length = read_seqs_to_hash(\%seqs,$path.'/'.$file.''.$params->{'ORTHOGROUP'}{'PROTEIN_ALIGN'}->[0],$params->{'ORTHOGROUP'}{'TAXA'},'aln');


  my $gene_align_id;
  my %aln_seqs;
  $gene_align_id = add_gene_align($dbh,$aln_method,$aln_length);
  foreach my $sp (keys %seqs){
    foreach my $tsl_stable_id (keys %{$seqs{$sp}}){
      next unless $seqs{$sp}{$tsl_stable_id}{'protein'};
      my $seqstr = $seqs{$sp}{$tsl_stable_id}{'protein'};
      my $length = length $seqstr;
      my $taxon_id;
      my $sp_name;
      my $seq_member_id;
      if (my $value = $params->{'TAXA'}{$sp}){
        if (ref $value){
          ($sp_name,$taxon_id) = @$value;
          # TODO: can't ignore these forever...
        }
        elsif (!$core_dbs->{$sp}) {
          fetch_meta_from_core_db($dbh,$core_dbs,$sp,$value,$params);
        }
      }
      if ($core_dbs->{$sp}) {
        $seq_member_id = fetch_gene_from_core_db($dbh,$core_dbs,$sp,$tsl_stable_id,$seqstr,$length);
      }
      if ($seq_member_id){
        # load seqs from from fna and faa.bounded files
        add_other_member_sequence($dbh,$seq_member_id,'cds',$seqs{$sp}{$tsl_stable_id}{'fna'}) if $seqs{$sp}{$tsl_stable_id}{'fna'};
        add_other_member_sequence($dbh,$seq_member_id,'exon_bounded',$seqs{$sp}{$tsl_stable_id}{'bounded'}) if $seqs{$sp}{$tsl_stable_id}{'bounded'};
        add_gene_align_member($dbh,$seq_member_id,$gene_align_id,$seqs{$sp}{$tsl_stable_id}{'aln'});
        $aln_seqs{$tsl_stable_id} = $seqs{$sp}{$tsl_stable_id}{'aln'};
      }
    }
  }
  my $species_set_id = get_species_set_id($dbh,'default');
  my $method_link_id = 401;
  my $mlss_id = add_method_link_species_set($dbh,$method_link_id,$species_set_id,'default','GenomeHubs');
  my $genetree = add_gene_tree ($params,$path.'/'.$file.''.$params->{'ORTHOGROUP'}{'TREE'}, 'protein', 'tree', 'default', $mlss_id, $gene_align_id, $cluster_id, 1, $st_nodes);
#  my $st_nodes = fetch_species_tree_nodes ($params);

  add_homology ($dbh,$params,$genetree,$st_nodes,$path.'/'.$file.''.$params->{'ORTHOGROUP'}{'HOMOLOG'},\%aln_seqs);

  return 1;
}

sub pairwise_identity {
  my $seqA = shift;
  my $seqB = shift;
  my @seqA = split //,$seqA;
  my @seqB = split //,$seqB;
  my ($perc_idA, $perc_covA, $perc_idB, $perc_covB, $lenA, $lenB);
  my $id = 0;
  my $aln = 0;
  for (my $i = 0; $i < scalar @seqA; $i++){
    if ($seqA[$i] ne '-'){
      $lenA++;
      if ($seqA[$i] eq $seqB[$i]){
        $id++;
        $aln++;
      }
      elsif ($seqB[$i] ne '-'){
        $aln++;
      }
    }
    if ($seqB[$i] ne '-'){
      $lenB++;
    }
  }
  $perc_idA = $id / $lenA * 100;
  $perc_covA = $aln / $lenA * 100;
  $perc_idB = $id / $lenB * 100;
  $perc_covB = $aln / $lenB * 100;
  return {perc_idA => $perc_idA, perc_idB => $perc_idB, perc_covA => $perc_covA, perc_covB => $perc_covB};
}

sub add_homology {
  my ($dbh,$params,$genetree, $st_nodes, $notung_homolog_filename, $seqs) = @_;
  my $cdba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host => $params->{'DATABASE_COMPARA'}{'HOST'},
    -user => $params->{'DATABASE_COMPARA'}{'RW_USER'},
    -pass => $params->{'DATABASE_COMPARA'}{'RW_PASS'},
    -port => $params->{'DATABASE_COMPARA'}{'PORT'},
    -dbname => $params->{'DATABASE_COMPARA'}{'NAME'},
  );
  my $gtna   = $cdba->get_adaptor("GeneTreeNode");

  my $homology_pair = parse_notung_homolog ($params,$notung_homolog_filename, 6); #6=ignore_lines

  $genetree = bless $genetree, 'Bio::EnsEMBL::Compara::GeneTree';
  my $gene_tree_root_id = $genetree->root_id();
  my $description = "";
  my $species_tree_node_id;
  my $gene_tree_node_id;

  # get all pairwise gene tree node ancestors and store in hash
  my %gt_nodes;
  my $leaves = $genetree->get_all_sorted_leaves();
  foreach my $g1 (@{$leaves}){
    $g1 = bless $g1, 'Bio::EnsEMBL::Compara::GeneTreeNode';
    $gt_nodes{$g1->name()}{$g1->name()} = $g1->node_id();
    foreach my $g2 (@{$leaves}){
      $g2 = bless $g2, 'Bio::EnsEMBL::Compara::GeneTreeNode';
      my $node_id = $g1->find_first_shared_ancestor($g2)->node_id();
      $gt_nodes{$g1->name()}{$g2->name()} = $node_id;
      $gt_nodes{$g2->name()}{$g2->name()} = $node_id;
    }
  }

  my $taxa = $params->{'ORTHOGROUP'}{'TAXA'};

  for my $seq1 (sort keys %{$homology_pair}) {
    for my $seq2 (sort keys %{$homology_pair->{$seq1}}) {

      my ($tax1,$seqm1) = ($1,$2) if $seq1 =~ /^($taxa).(\S+)$/;
      my ($tax2,$seqm2) = ($1,$2) if $seq2 =~ /^($taxa).(\S+)$/;
      my ($gtn1_id,$gtn2_id);
      my $sp1 = $params->{'TAXA'}{$tax1};
      $sp1 =~ s/_core_.+$//;
      my $sp2 = $params->{'TAXA'}{$tax2};
      $sp2 =~ s/_core_.+$//;
      $description = $homology_pair->{$seq1}{$seq2}{description};

      my $species_set_name = $sp1 eq $sp2 ? $sp1 : join('-',(sort(($sp1,$sp2))));
      my $species_set_id = get_species_set_id($dbh,$species_set_name);
      my $aln1 = $seqs->{$seqm1};
      my $aln2 = $seqs->{$seqm2};
      my $identities = pairwise_identity($aln1, $aln2);
      my ($method_link_id,$mlss_name,$homology_type);
      if ($description =~ m/ortholog/i){
        $method_link_id = 201;
        $mlss_name = $species_set_name.'-Orthologues';
        $homology_type = 'Orthologues';
      }
      else {
        $method_link_id = 202;
        $mlss_name = $species_set_name.'-Paralogues';
        $homology_type = 'Paralogues';
      }

      my $mlss_id = add_method_link_species_set($dbh,$method_link_id,$species_set_id,$mlss_name,'lepbase');

      # get gene_tree_node id
      $gene_tree_node_id = $gt_nodes{$seq1}{$seq2};

      # get species_tree_node id
      $species_tree_node_id = $st_nodes->{$tax1}{$tax2};

      #insert into homology table
      $dbh->do("INSERT INTO homology (method_link_species_set_id,description,is_tree_compliant,species_tree_node_id,gene_tree_node_id,gene_tree_root_id)"
              ." VALUES (".$mlss_id
              .",".$dbh->quote($description)
              .","."1"
              .",".$species_tree_node_id
              .",".$gene_tree_node_id
              .",".$gene_tree_root_id
              .")");
      my $homology_id = $dbh->last_insert_id(undef,undef,undef,undef);

      add_homology_members($dbh,$params,$homology_id,$seqm1,$seqm2,$homology_type,$identities);

    }
  }

}

sub add_homology_members {
  my ($dbh,$params,$homology_id,$seqm1,$seqm2,$homology_type,$identities) = @_;
  my $cdba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host => $params->{'DATABASE_COMPARA'}{'HOST'},
    -user => $params->{'DATABASE_COMPARA'}{'RW_USER'},
    -pass => $params->{'DATABASE_COMPARA'}{'RW_PASS'},
    -port => $params->{'DATABASE_COMPARA'}{'PORT'},
    -dbname => $params->{'DATABASE_COMPARA'}{'NAME'},
  );
  my $seqma = $cdba->get_adaptor("SeqMember");

  my $sm1 = $seqma->fetch_by_stable_id($seqm1);
  my $gm1 = $sm1->gene_member();
  my $seq_member_id_1  = $sm1->dbID();
  my $gene_member_id_1 = $gm1->dbID();
  my $cigar1 = "";
  if ($seq_member_id_1) {
    my $sth = $dbh->prepare("SELECT cigar_line FROM gene_align_member WHERE seq_member_id = $seq_member_id_1");
    $sth->execute();
    if ($sth->rows > 0){
      $cigar1 = $sth->fetchrow_arrayref()->[0]
    } else {
      die "Could not find SeqMember $seqm1 with dbID $seq_member_id_1 in gene_align_member";
    }
  }

  my $sm2 = $seqma->fetch_by_stable_id($seqm2);
  my $gm2 = $sm2->gene_member();
  my $seq_member_id_2  = $sm2->dbID();
  my $gene_member_id_2 = $gm2->dbID();
  my $cigar2 = "";
  if ($seq_member_id_2) {
    my $sth = $dbh->prepare("SELECT cigar_line FROM gene_align_member WHERE seq_member_id = $seq_member_id_2");
    $sth->execute();
    if ($sth->rows > 0){
      $cigar2 = $sth->fetchrow_arrayref()->[0]
    } else {
      die "Could not find SeqMember $seqm2 with dbID $seq_member_id_2 in gene_align_member";
    }
  }

  ($cigar1,$cigar2) = Bio::EnsEMBL::Compara::Utils::Cigars::minimize_cigars($cigar1,$cigar2);

  #insert into homology_member table
  $dbh->do("INSERT INTO homology_member (homology_id,gene_member_id,seq_member_id,cigar_line,perc_id,perc_cov)"
          ." VALUES (".$homology_id
          .",".$gene_member_id_1
          .",".$seq_member_id_1
          .",".$dbh->quote($cigar1)
          .",".$identities->{perc_idA}
          .",".$identities->{perc_covA}
          .")");
  $dbh->do("INSERT INTO homology_member (homology_id,gene_member_id,seq_member_id,cigar_line,perc_id,perc_cov)"
          ." VALUES (".$homology_id
          .",".$gene_member_id_2
          .",".$seq_member_id_2
          .",".$dbh->quote($cigar2)
          .",".$identities->{perc_idB}
          .",".$identities->{perc_covB}
          .")");

  # increment orthologue/paralogue counts
  $dbh->do("UPDATE gene_member_hom_stats SET `$homology_type` = `$homology_type` + 1 WHERE gene_member_id = $gene_member_id_1");
  $dbh->do("UPDATE gene_member_hom_stats SET `$homology_type` = `$homology_type` + 1 WHERE gene_member_id = $gene_member_id_2");
}

sub parse_notung_homolog {
  my ($params,$notung_homolog_filename,$ignore_lines) = @_;
  my %tax_count;
  my %homology_pair;
  my @seq_members;

  use List::Util qw(min max);

  open NOTUNG_HOMOLOG, "<$notung_homolog_filename" or die $!;

  for (1..$ignore_lines) {
    <NOTUNG_HOMOLOG>
  };

  my $all_seqs = <NOTUNG_HOMOLOG>;

  my $taxa = $params->{'ORTHOGROUP'}{'TAXA'};

  while ($all_seqs =~ /($taxa).(\S+)/g) {
    $tax_count{$1}++;
    push @seq_members, "$1_$2";
  }

  while (<NOTUNG_HOMOLOG>) {
    chomp;
    my @row  = split("\t");
    my @seq_members_local = @seq_members;
    my $seq1 = shift @row;
    my ($tax1, $seqm1) = ($1, $2) if $seq1 =~ /^($taxa).(\S+)$/;
    while (my $homolog_value = shift @row) {
      my $seq2 = shift @seq_members_local;
      my ($tax2, $seqm2) = ($1, $2) if $seq2 =~ /^($taxa).(\S+)$/;
      next if $seq1 eq $seq2;
      next if $homology_pair{$seq2}{$seq1};

      if ($homolog_value eq "O") {
        if ($tax_count{$tax1} == 1 and $tax_count{$tax2} == 1) {
          $homology_pair {$seq1}{$seq2}{description} = "ortholog_one2one";
        } elsif ( min($tax_count{$tax1},$tax_count{$tax2}) == 1 and max($tax_count{$tax1},$tax_count{$tax2}) > 1) {
          $homology_pair {$seq1}{$seq2}{description} = "ortholog_one2many";
        } elsif ( min($tax_count{$tax1},$tax_count{$tax2})  > 1 and max($tax_count{$tax1},$tax_count{$tax2}) > 1 ) {
          $homology_pair {$seq1}{$seq2}{description} = "ortholog_many2many";
        }
      } elsif ($homolog_value eq "P") {
        $homology_pair {$seq1}{$seq2}{description} = ($tax1 eq $tax2) ? "within_species_paralog" : "between_species_paralog";
      }
    }
  }

  return \%homology_pair;

  # test hash by printing:
  #   for my $seq1 (sort keys %homology_pair) {
  #     for my $seq2 (sort keys %{$homology_pair{$seq1}}) {
  #       print "$seq1\t$seq2\t$homology_pair{$seq1}{$seq2}{description}\n";
  #     }
  #   }

}

sub read_seqs_to_hash {
  my ($seqs,$file,$taxa,$key) = @_;
  my $seqin = Bio::SeqIO->new(-file => $file, -format => 'fasta');
  my $len = -1;
  while (my $seq = $seqin->next_seq()){
    $seq->display_id() =~ m/^($taxa).(.+)$/;
    my $sp = $1;
    my $tsl_stable_id = $2;
    $seqs->{$sp}{$tsl_stable_id}{$key} = $seq->seq();
    $len = $seq->length() if $seq->length() > $len;
  }
  return $len;
}

sub fetch_species_tree_nodes {
  my ($params,$dbh) = @_;
  my %st_nodes;
  my $cdba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host => $params->{'DATABASE_COMPARA'}{'HOST'},
    -user => $params->{'DATABASE_COMPARA'}{'RW_USER'},
    -pass => $params->{'DATABASE_COMPARA'}{'RW_PASS'},
    -port => $params->{'DATABASE_COMPARA'}{'PORT'},
    -dbname => $params->{'DATABASE_COMPARA'}{'NAME'},
  );
  my $sta = $cdba->get_adaptor("SpeciesTree");
  my $ncbita = $cdba->get_adaptor("NCBITaxon");

  my $speciestree = $sta->fetch_by_root_id(1000);
  my $speciestree_root = $speciestree->root();
  $speciestree_root = bless $speciestree_root, 'Bio::EnsEMBL::Compara::SpeciesTreeNode';

  my $leaves = $speciestree_root->get_all_sorted_leaves();
  my %core_dbs;
  while (my $sp1 = shift @{$leaves}){
    if ($sp1->name()){
      my $production_name = $params->{'TAXA'}{$sp1->name()};
      $production_name =~ s/_core_.+$//;
      my $genome_db_id = fetch_meta_from_core_db($dbh,\%core_dbs,$sp1->name(),$params->{'TAXA'}{$sp1->name()},$params);
      my $taxon_id = $core_dbs{$sp1->name()}{'taxon_id'};
      $sp1->taxon_id($taxon_id);
      $dbh->do("update species_tree_node set genome_db_id = ".$genome_db_id.", taxon_id = ".$taxon_id." where node_name = ".$dbh->quote($sp1->name()));
      $st_nodes{$sp1->name()}{$sp1->name()} = $sp1->node_id();
    }
    foreach my $sp2 (@{$leaves}){
      my $node_id = $sp1->find_first_shared_ancestor($sp2)->node_id();
      $st_nodes{$sp1->name()}{$sp2->name()} = $node_id;
      $st_nodes{$sp2->name()}{$sp1->name()} = $node_id;
    }
  }

  my $nodes = $speciestree_root->get_all_nodes();
  while (my $node = shift @{$nodes}){
    next if $node->is_leaf();
    my @leaf_tax_nodes;
    my $node_leaves = $node->get_all_sorted_leaves();
    for my $node_leaf (@{$node_leaves}) {
      push @leaf_tax_nodes, $ncbita->fetch_node_by_taxon_id($node_leaf->taxon_id());
    }
    my $tax_node = $ncbita->fetch_first_shared_ancestor_indexed(@leaf_tax_nodes);
    $dbh->do("update species_tree_node set taxon_id = ". $tax_node->taxon_id() ." where node_id = ".$dbh->quote($node->node_id()));
  }

  return (\%st_nodes,\%core_dbs);;
}

sub add_species_tree {
  my ($params,$newick_treefile, $mlss_id) = @_;
  my $cdba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host => $params->{'DATABASE_COMPARA'}{'HOST'},
    -user => $params->{'DATABASE_COMPARA'}{'RW_USER'},
    -pass => $params->{'DATABASE_COMPARA'}{'RW_PASS'},
    -port => $params->{'DATABASE_COMPARA'}{'PORT'},
    -dbname => $params->{'DATABASE_COMPARA'}{'NAME'},
  );
  my $sta   = $cdba->get_adaptor("SpeciesTree");
  my $stna  = $cdba->get_adaptor("SpeciesTreeNode");

  my $species_newick = "";
  open SPECIES_NEWICK, $params->{'SPECIES_SET'}{'TREE_FILE'};
  while (<SPECIES_NEWICK>) {
    chomp; $species_newick .= $_;
  }
  my $speciestree = $sta->fetch_by_root_id(1000);
  my $newroot = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($species_newick, "Bio::EnsEMBL::Compara::SpeciesTreeNode");
  $speciestree->root($newroot);
  $newroot->build_leftright_indexing;
  $newroot->distance_to_parent;


  $sta->store($speciestree);
}

sub add_gene_tree {

  my ($params,$newick_treefile, $member_type, $tree_type, $clusterset_id, $mlss_id, $gene_align_id, $stable_id, $version, $st_nodes) = @_;
  my $cdba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host => $params->{'DATABASE_COMPARA'}{'HOST'},
    -user => $params->{'DATABASE_COMPARA'}{'RW_USER'},
    -pass => $params->{'DATABASE_COMPARA'}{'RW_PASS'},
    -port => $params->{'DATABASE_COMPARA'}{'PORT'},
    -dbname => $params->{'DATABASE_COMPARA'}{'NAME'},
  );
  my $sta   = $cdba->get_adaptor("SpeciesTree");
  my $gta   = $cdba->get_adaptor("GeneTree");
  my $gtna  = $cdba->get_adaptor("GeneTreeNode");
  my $seqma = $cdba->get_adaptor("SeqMember");

  my $newtree = new Bio::EnsEMBL::Compara::GeneTree( -member_type => $member_type, -tree_type => $tree_type, -clusterset_id => $clusterset_id,
    -method_link_species_set_id => $mlss_id, -gene_align_id => $gene_align_id, -stable_id => $stable_id, -version => $version);

  $gta->store($newtree);

  open TREEFILE, "<$newick_treefile" or die "Could not open newick treefile $newick_treefile\n";
  chomp(my $newick_tree = <TREEFILE>);
  my $newroot = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($newick_tree, "Bio::EnsEMBL::Compara::GeneTreeNode");
  $newroot->build_leftright_indexing;

  $newroot->node_id($newtree->root_id);
  $newroot->distance_to_parent($newtree->root->distance_to_parent);
  $newroot->adaptor($newtree->root->adaptor);
  $newroot->tree($newtree);
  $newtree->{'_root'} = $newroot;

  my $supertree = $gta->fetch_by_root_id(1);
  $supertree->root->add_child($newroot);
  my $species_tree = $sta->fetch_by_root_id(1000);
  my $species_tree_root = $species_tree->root();
#  $supertree->species_tree($species_tree);
  my $taxa = $params->{'ORTHOGROUP'}{'TAXA'};
  foreach my $node (@{$newroot->get_all_nodes}) {
    #print $node->name . "\n" if defined $node->name;
    if ($node->is_leaf) {
      my $node_name = $node->name();

      $node_name =~ m/^($taxa).(.+)$/;
      my $taxon = $1;
      $node->name($node_name);
      my $seqm = $seqma->fetch_by_stable_id($2);
      $node = bless $node, 'Bio::EnsEMBL::Compara::GeneTreeMember';
      $node->seq_member_id($seqm->seq_member_id) if $seqm;
      $node->add_tag('species_tree_node_id',$st_nodes->{$taxon}{$taxon});
    } else {
      if (defined $node->get_tagvalue("Duplication")) {
        my $is_dup = $node->get_tagvalue("Duplication");
        $node->add_tag('is_dup', $is_dup );
        $node->add_tag('node_type',$is_dup == 1 ? 'duplication' : 'speciation');
        if (defined $node->get_tagvalue("l")) {
          my $dup_score = $node->get_tagvalue("l");
          $node->add_tag('duplication_confidence_score', $dup_score );
        }
      }
      my $leaves = $node->get_all_sorted_leaves();
      my %species;
      while (my $leaf = shift @{$leaves}){
        $species{$leaf->get_tagvalue("s")} = 1;
      }
      my @leaf_list;
      foreach my $sp (keys %species){

        push @leaf_list, $species_tree_root->find_leaf_by_name($sp)
      }
      my $spt_node_id;
      my $first = shift @leaf_list;
      if (scalar @leaf_list > 0){
        $spt_node_id = $first->find_first_shared_ancestor_from_leaves(\@leaf_list)->node_id;
      }
      else {
        $spt_node_id = $first->node_id();
      }
      $node->add_tag('species_tree_node_id',$spt_node_id);
#     if (defined $node->name and $node->name =~ /(\d*)(Y|N)/) {
#       $node->add_tag('bootstrap',$1);
#       $node->add_tag('is_dup',   $2 eq "Y" ? 1 : 0);
#       $node->add_tag('node_type',$2 eq "Y" ? 'duplication' : 'speciation');
#     }
    }
  }

  $gta->store($supertree);

  my $all_nodes = $newtree->get_all_nodes;
  $gtna->_store_all_tags($all_nodes);
  return $newtree;
}

sub add_gene_tree_root {
  my ($dbh,$member_type,$tree_type,$clusterset_id,$mlss_id,$gene_align_id) = @_;
  if ($gene_align_id){
    my $sth = $dbh->prepare("SELECT root_id FROM gene_tree_root WHERE member_type = ".$dbh->quote($member_type)." AND tree_type = ".$dbh->quote($tree_type)." AND clusterset_id = ".$dbh->quote($clusterset_id)." AND method_link_species_set_id = $mlss_id AND gene_align_id = $gene_align_id");
    $sth->execute();
    if ($sth->rows > 0){
      return $sth->fetchrow_arrayref()->[0]
    }
    $dbh->do("INSERT INTO gene_tree_root (member_type,tree_type,clusterset_id,method_link_species_set_id,gene_align_id)"
              ." VALUES (".$dbh->quote($member_type)
              .",".$dbh->quote($tree_type)
              .",".$dbh->quote($clusterset_id)
              .",".$mlss_id
              .",".$gene_align_id
              .")");
    $sth->execute();
    return $sth->fetchrow_arrayref()->[0];
  }
  else {
    my $sth = $dbh->prepare("SELECT root_id FROM gene_tree_root WHERE member_type = ".$dbh->quote($member_type)." AND tree_type = ".$dbh->quote($tree_type)." AND clusterset_id = ".$dbh->quote($clusterset_id)." AND method_link_species_set_id = $mlss_id");
    $sth->execute();
    if ($sth->rows > 0){
      return $sth->fetchrow_arrayref()->[0]
    }
    $dbh->do("INSERT INTO gene_tree_root (member_type,tree_type,clusterset_id,method_link_species_set_id)"
              ." VALUES (".$dbh->quote($member_type)
              .",".$dbh->quote($tree_type)
              .",".$dbh->quote($clusterset_id)
              .",".$mlss_id
              .")");
    $sth->execute();
    return $sth->fetchrow_arrayref()->[0];
  }
  return;
}


sub add_method_link_species_set {
  my ($dbh,$method_link_id,$species_set_id,$name,$source) = @_;
  my $sth = $dbh->prepare("SELECT method_link_species_set_id FROM method_link_species_set WHERE method_link_id = $method_link_id AND species_set_id = $species_set_id");
  $sth->execute();
  if ($sth->rows > 0){
    return $sth->fetchrow_arrayref()->[0]
  }
  $dbh->do("INSERT INTO method_link_species_set (method_link_id,species_set_id,name,source)"
            ." VALUES (".$method_link_id
            .",".$species_set_id
            .",".$dbh->quote($name)
            .",".$dbh->quote($source)
            .")");
  $sth->execute();
  return $sth->fetchrow_arrayref()->[0]
}

sub get_species_set_id {
  my ($dbh,$value) = @_;
  my $species_set_id;
  my $sth = $dbh->prepare("SELECT species_set_id FROM species_set_tag WHERE tag = 'name' AND value = ".$dbh->quote($value));
  $sth->execute();
  if ($sth->rows > 0){
    $species_set_id = $sth->fetchrow_arrayref()->[0]
  }

  return $species_set_id;
}

sub add_gene_align_member {
  my ($dbh,$seq_member_id,$gene_align_id,$seqstr) = @_;
  my @arr = split /(-+)/,$seqstr;
  my $cigar;
  foreach my $chunk (@arr){
    if (my $length = length($chunk)){
      if ($chunk =~ m/-/){
        $cigar .= $length.'D';
      }
      else {
        $cigar .= $length.'M';
      }
    }
  }
  my $sth = $dbh->prepare("SELECT cigar_line FROM gene_align_member WHERE seq_member_id = $seq_member_id AND gene_align_id = $gene_align_id");
  $sth->execute();
  if ($sth->rows > 0){
    return 1;
  }
  $dbh->do("INSERT INTO gene_align_member (seq_member_id,gene_align_id,cigar_line)"
            ." VALUES (".$seq_member_id
            .",".$gene_align_id
            .",".$dbh->quote($cigar)
            .")");
  return 1;
}


sub add_gene_align {
  my ($dbh,$aln_method,$aln_length) = @_;
  $dbh->do("INSERT INTO gene_align (aln_method,aln_length)"
            ." VALUES (".$dbh->quote($aln_method)
            .",".$aln_length
            .")");
  return $dbh->last_insert_id(undef,undef,undef,undef);
}

sub add_other_member_sequence {
  my ($dbh,$seq_member_id,$seq_type,$sequence) = @_;
  my $length = length $sequence;
  my $sth = $dbh->prepare("SELECT sequence,length FROM other_member_sequence WHERE seq_member_id = $seq_member_id AND seq_type = ".$dbh->quote($seq_type));
  $sth->execute();
  if ($sth->rows > 0){
    return 1;
  }
  $dbh->do("INSERT INTO other_member_sequence (seq_member_id,seq_type,length,sequence)"
            ." VALUES (".$seq_member_id
            .",".$dbh->quote($seq_type)
            .",".$length
            .",".$dbh->quote($sequence)
            .")");
  return 1;
}

sub fetch_gene_from_core_db {
  my ($dbh,$core_dbs,$sp,$tsl_stable_id,$seqstr,$length) = @_;
  my $sth = $dbh->prepare("SELECT seq_member_id,gene_member_id FROM seq_member WHERE genome_db_id = $core_dbs->{$sp}{'genome_db_id'} AND stable_id = ".$dbh->quote($tsl_stable_id));
  $sth->execute();
  if ($sth->rows > 0){
    my @arr = $sth->fetchrow_array();
    return @arr;
  }
  #$dbh->do("INSERT INTO sequence (length, sequence) VALUES ($length,".$dbh->quote($seqstr).")");
  #$dbh->do("INSERT INTO seq_member (length, sequence) VALUES ($length,".$dbh->quote($seqstr).")");
  my $cdbh = $core_dbs->{$sp}{'db_handle'};
  # TODO: change this to translation.stable_id
  my $csth = $cdbh->prepare("SELECT cs.name cs_name, sr.name sr_name, sr.length sr_length,
                                    tsc.transcript_id, tsc.seq_region_start tsc_sr_start, tsc.seq_region_end tsc_sr_end, tsc.seq_region_strand tsc_sr_strand, tsc.description tsc_description, tsc.display_xref_id tsc_xref_id,
                                    g.seq_region_start g_sr_start, g.seq_region_end g_sr_end, g.seq_region_strand g_sr_strand, g.stable_id g_stable_id, g.description g_description, g.canonical_transcript_id, g.display_xref_id g_xref_id, g.gene_id,
                                    tsl.seq_start tsl_seq_start, tsl.seq_end tsl_seq_end, tsl.stable_id tsl_stable_id,
                                    se.seq_region_start se_sr_start, se.seq_region_end se_sr_end, ee.seq_region_start ee_sr_start, ee.seq_region_end ee_sr_end
                                  FROM translation tsl
                                  JOIN exon se ON tsl.start_exon_id = se.exon_id
                                  JOIN exon ee ON tsl.end_exon_id = ee.exon_id
                                  JOIN transcript tsc ON tsl.transcript_id = tsc.transcript_id
                                  JOIN gene g ON tsc.gene_id = g.gene_id
                                  JOIN seq_region sr ON tsc.seq_region_id = sr.seq_region_id
                                  JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id
                  WHERE tsl.stable_id LIKE ".$dbh->quote($tsl_stable_id)
                );
  $csth->execute();
  #print $tsl_stable_id,"!\n";
  if (my $ref = $csth->fetchrow_hashref()){
    my $dnafrag_id = add_dnafrag($dbh,$core_dbs->{$sp}{'genome_db_id'},$ref->{'cs_name'},$ref->{'sr_name'},$ref->{'sr_length'});
    #print $dnafrag_id,"\t",$ref->{'tsl_stable_id'},"\n";
    my $display_label = 'NULL'; #TODO: look for a display label (display_xref...)
    my $gene_member_id = add_gene_member($dbh,$core_dbs,$sp,$core_dbs->{$sp}{'genome_db_id'},$ref->{'g_stable_id'},'default',$core_dbs->{$sp}{'taxon_id'},$ref->{'g_description'},$dnafrag_id,$ref->{'g_sr_start'},$ref->{'g_sr_end'},$ref->{'g_sr_strand'},$ref->{'gene_id'},$ref->{'g_xref_id'},$display_label);
    #print $gene_member_id,"\t",$ref->{'g_stable_id'},"\n";
    my $sequence_id = add_sequence($dbh,$seqstr,$length);
    #print $sequence_id," = seq_id\n";
    my ($sr_start,$sr_end);
    if ($ref->{'tsc_sr_strand'} != -1){
      $sr_start = $ref->{'se_sr_start'} + $ref->{'tsl_seq_start'} - 1;
      $sr_end = $ref->{'ee_sr_start'} + $ref->{'tsl_seq_end'} - 1;
    }
    else {
      $sr_start = $ref->{'ee_sr_end'} - $ref->{'tsl_seq_start'} + 1;
      $sr_end = $ref->{'se_sr_end'} - $ref->{'tsl_seq_end'} + 1;
    }
    my $seq_member_id = add_seq_member($dbh,$core_dbs,$sp,$core_dbs->{$sp}{'genome_db_id'},$ref->{'tsl_stable_id'},'ENSEMBLPEP',$core_dbs->{$sp}{'taxon_id'},$sequence_id,$gene_member_id,$ref->{'g_description'},$dnafrag_id,$sr_start,$sr_end,$ref->{'tsc_sr_strand'},$ref->{'transcript_id'},$ref->{'tsc_xref_id'},$display_label);
    #print $seq_member_id,"\t",$ref->{'g_stable_id'},"\n";
    if ($ref->{'canonical_transcript_id'} == $ref->{'transcript_id'} ){
      simple_update($dbh,'gene_member',{'canonical_member_id' => $seq_member_id},{'gene_member_id' => $gene_member_id});
    }
    return $seq_member_id;
  }
}

sub add_seq_member {
  my ($dbh,$core_dbs,$sp,$genome_db_id,$stable_id,$source_name,$taxon_id,$sequence_id,$gene_member_id,$description,$dnafrag_id,$dnafrag_start,$dnafrag_end,$dnafrag_strand,$transcript_id,$display_xref_id,$display_label) = @_;
  my $sth = $dbh->prepare("SELECT seq_member_id from seq_member WHERE genome_db_id = $genome_db_id AND stable_id = ".$dbh->quote($stable_id));
  $sth->execute();
  if ($sth->rows() > 0){
    return $sth->fetchrow_arrayref()->[0];
  }
  $description = $description ? $dbh->quote($description) : 'NULL';
  $display_label = $dbh->quote($display_label) unless $display_label eq 'NULL';
  $dbh->do("INSERT INTO seq_member (stable_id,source_name,taxon_id,genome_db_id,sequence_id,gene_member_id,description,dnafrag_id,dnafrag_start,dnafrag_end,dnafrag_strand,display_label)"
          ." VALUES (".$dbh->quote($stable_id)
          .",".$dbh->quote($source_name)
          .",".$taxon_id
          .",".$genome_db_id
          .",".$sequence_id
          .",".$gene_member_id
          .",".$description
          .",".$dnafrag_id
          .",".$dnafrag_start
          .",".$dnafrag_end
          .",".$dnafrag_strand
          .",".$display_label
          .")");
  my $seq_member_id;
  $sth->execute();
  if ($sth->rows() > 0){
    $seq_member_id = $sth->fetchrow_arrayref()->[0];

    my $cdbh = $core_dbs->{$sp}{'db_handle'};
    my $csth = $cdbh->prepare("SELECT x.xref_id, x.display_label
                            FROM xref x
                            JOIN object_xref o ON o.xref_id = x.xref_id
                            WHERE ensembl_id = $transcript_id
                            AND ensembl_object_type = 'transcript'");
    $csth->execute();
    while (my $ref = $csth->fetchrow_hashref()){
      if ($display_xref_id && $ref->{'xref_id'} == $display_xref_id && $ref->{'display_label'}){
        # update display_label
        simple_update($dbh,'seq_member',{'display_label' => $dbh->quote($ref->{'display_label'})},{'seq_member_id' => $gene_member_id});
      }
    }

  }
  return $seq_member_id;
}

sub add_sequence {
  my ($dbh,$seqstr,$length) = @_;
  $dbh->do("INSERT INTO sequence (length,sequence)"
          ." VALUES ($length"
          .",".$dbh->quote($seqstr)
          .")");
  return $dbh->last_insert_id(undef,undef,undef,undef);
}

sub add_gene_member {
  my ($dbh,$core_dbs,$sp,$genome_db_id,$stable_id,$source_name,$taxon_id,$description,$dnafrag_id,$dnafrag_start,$dnafrag_end,$dnafrag_strand,$gene_id,$display_xref_id,$display_label) = @_;
  my $sth = $dbh->prepare("SELECT gene_member_id from gene_member WHERE genome_db_id = $genome_db_id AND stable_id = ".$dbh->quote($stable_id));
  $sth->execute();
  if ($sth->rows() > 0){
    return $sth->fetchrow_arrayref()->[0];
  }

  $description = $description ? $dbh->quote($description) : 'NULL';
  $display_label = $dbh->quote($display_label) unless $display_label eq 'NULL';
  my $gene_trees = 1;
  # TODO: Update gene_trees, gene_gain_loss_trees, orthologues, paralogues when adding trees
  $dbh->do("INSERT INTO gene_member (stable_id,source_name,taxon_id,genome_db_id,description,dnafrag_id,dnafrag_start,dnafrag_end,dnafrag_strand,display_label)"
          ." VALUES (".$dbh->quote($stable_id)
          .",".$dbh->quote($source_name)
          .",".$taxon_id
          .",".$genome_db_id
          .",".$description
          .",".$dnafrag_id
          .",".$dnafrag_start
          .",".$dnafrag_end
          .",".$dnafrag_strand
          .",".$display_label
          .")");
  my $gene_member_id;
  $sth->execute();
  if ($sth->rows() > 0){
    $gene_member_id = $sth->fetchrow_arrayref()->[0];

    $dbh->do("INSERT INTO gene_member_hom_stats (gene_member_id,collection,gene_trees,gene_gain_loss_trees)"
            ." VALUES (".$gene_member_id
            .",".$dbh->quote($source_name)
            .",1"
            .",1"
            .")");

    # TODO: fetch an xref and fill in display_xref if the gene has a display_xref_id

    my $cdbh = $core_dbs->{$sp}{'db_handle'};
    my $csth = $cdbh->prepare("SELECT x.xref_id, x.dbprimary_acc, x.display_label,
                                      e.db_name db_name, e.db_release, e.status, e.priority, e.db_display_name, e.type, e.secondary_db_name,e.secondary_db_table, e.description e_desc
                            FROM xref x
                            JOIN external_db e ON x.external_db_id = e.external_db_id
                            JOIN object_xref o ON o.xref_id = x.xref_id
                            WHERE ensembl_id = $gene_id
                            AND ensembl_object_type = 'gene'");
    $csth->execute();
    while (my $ref = $csth->fetchrow_hashref()){
      my $external_db_id = add_external_db($dbh,$ref->{'db_name'},$ref->{'db_release'},$ref->{'status'},$ref->{'priority'},$ref->{'db_display_name'},$ref->{'type'},$ref->{'secondary_display_name'},$ref->{'secondary_db_table'},$ref->{'e_desc'});
      add_member_xref($dbh,$gene_member_id,$ref->{'dbprimary_acc'},$external_db_id);
      if ($display_xref_id && $ref->{'xref_id'} == $display_xref_id && $ref->{'display_label'}){
        # update display_label
        simple_update($dbh,'gene_member',{'display_label' => $dbh->quote($ref->{'display_label'})},{'gene_member_id' => $gene_member_id});
      }
    }

  }

  return $gene_member_id;
}

sub add_member_xref {
  my ($dbh,$gene_member_id,$dbprimary_acc,$external_db_id) = @_;
  my $sth = $dbh->prepare("SELECT dbprimary_acc FROM member_xref WHERE gene_member_id = $gene_member_id AND dbprimary_acc = ".$dbh->quote($dbprimary_acc)." AND external_db_id = $external_db_id");
  $sth->execute();
  if ($sth->rows() > 0){
    return $sth->fetchrow_arrayref()->[0];
  }
  $dbh->do("INSERT INTO member_xref (gene_member_id,dbprimary_acc,external_db_id)"
          ." VALUES (".$gene_member_id
          .",".$dbh->quote($dbprimary_acc)
          .",".$external_db_id
          .")");
  return 1;
}

sub add_external_db {
  my ($dbh,$db_name,$db_release,$status,$priority,$db_display_name,$type,$secondary_db_name,$secondary_db_table,$description) = @_;
  #print $db_name,"\n";
  my $sth = $dbh->prepare("SELECT external_db_id FROM external_db WHERE db_name = ".$dbh->quote($db_name));
  $sth->execute();
  if ($sth->rows() > 0){
    return $sth->fetchrow_arrayref()->[0];
  }
  $dbh->do("INSERT INTO external_db (db_name,db_release,status,priority,db_display_name,type,secondary_db_name,secondary_db_table,description)"
          ." VALUES (".$dbh->quote($db_name)
          .",".$dbh->quote($db_release)
          .",".$dbh->quote($status)
          .",".$priority
          .",".$dbh->quote($db_display_name)
          .",".$dbh->quote($type)
          .",".$dbh->quote($secondary_db_name)
          .",".$dbh->quote($secondary_db_table)
          .",".$dbh->quote($description)
          .")");
  return $dbh->last_insert_id(undef,undef,undef,undef);
}


sub add_dnafrag {
  my ($dbh,$genome_db_id,$cs_name,$sr_name,$sr_length) = @_;
  my $sth = $dbh->prepare("SELECT dnafrag_id from dnafrag WHERE genome_db_id = $genome_db_id AND coord_system_name = ".$dbh->quote($cs_name)." AND name = ".$dbh->quote($sr_name));
  $sth->execute();
  if ($sth->rows() > 0){
    return $sth->fetchrow_arrayref()->[0];
  }
  $dbh->do("INSERT INTO dnafrag (length,name,genome_db_id,coord_system_name)"
          ." VALUES ($sr_length"
          .",".$dbh->quote($sr_name)
          .",".$genome_db_id
          .",".$dbh->quote($cs_name)
          .")");
  $sth->execute();
  if ($sth->rows() > 0){
    return $sth->fetchrow_arrayref()->[0];
  }
  return;

}

sub fetch_genome_db_id {
  my ($name,$assembly,$params) = @_;
  my $dbh = compara_db_connect($params);
  my $sth = $dbh->prepare("SELECT genome_db_id FROM genome_db WHERE name = ".$dbh->quote($name)." AND assembly = ".$dbh->quote($assembly));
  $sth->execute();
  if ($sth->rows > 0){
    return $sth->fetchrow_arrayref()->[0];
  }
  return;
}

sub fetch_meta_from_core_db {
  my ($dbh,$core_dbs,$sp,$db_name,$params) = @_;
  $core_dbs->{$sp}{'db_name'} = $db_name;
  my $cdbh = core_db_connect($db_name,$params);
  $core_dbs->{$sp}{'db_handle'} = $cdbh;
  my $csth = $cdbh->prepare("SELECT meta_value FROM meta WHERE meta_key LIKE ?");
  $csth->execute('species.production_name');
  if ($csth->rows > 0){
    $core_dbs->{$sp}{'name'} = lc($csth->fetchrow_arrayref()->[0]);
  }
  else {
    return;
  }
  my $sth = $dbh->prepare("SELECT genome_db_id,taxon_id,assembly,genebuild FROM genome_db WHERE name = ".$dbh->quote($core_dbs->{$sp}{'name'}));
  $sth->execute();
  if ($sth->rows > 0){
    my $ref = $sth->fetchrow_arrayref();
    $core_dbs->{$sp}{'genome_db_id'} = $ref->[0];
    $core_dbs->{$sp}{'taxon_id'} = $ref->[1];
    $core_dbs->{$sp}{'assembly'} = $ref->[2];
    $core_dbs->{$sp}{'genebuild'} = $ref->[3];
    return $core_dbs->{$sp}{'genome_db_id'};
  }
  $csth->execute('species.production_name');
  if ($csth->rows > 0){
    $core_dbs->{$sp}{'name'} = lc($csth->fetchrow_arrayref()->[0]);
  }
  else {
    return;
  }
  $csth->execute('species.taxonomy_id');
  if ($csth->rows > 0){
    $core_dbs->{$sp}{'taxon_id'} = $csth->fetchrow_arrayref()->[0];
  }
  else {
    return;
  }
  $csth->execute('assembly.default');
  if ($csth->rows > 0){
    $core_dbs->{$sp}{'assembly'} = $csth->fetchrow_arrayref()->[0];
  }
  else {
    return;
  }
  $csth->execute('genebuild.start_date');
  if ($csth->rows > 0){
    $core_dbs->{$sp}{'genebuild'} = $csth->fetchrow_arrayref()->[0];
  }
  else {
    return;
  }
  $core_dbs->{$sp}{'genome_db_id'} = fetch_genome_db_id($core_dbs->{$sp}{'name'},$core_dbs->{$sp}{'assembly'},$params);
  if ($core_dbs->{$sp}{'genome_db_id'}){
    $dbh->do("INSERT INTO genome_db (genome_db_id,taxon_id,name,assembly,genebuild) "
        ."VALUES (  ".$core_dbs->{$sp}{'genome_db_id'}
              .",".$core_dbs->{$sp}{'taxon_id'}
              .",".$dbh->quote($core_dbs->{$sp}{'name'})
              .",".$dbh->quote($core_dbs->{$sp}{'assembly'})
              .",".$dbh->quote($core_dbs->{$sp}{'genebuild'})
              .")");
  }
  else {
    $dbh->do("INSERT INTO genome_db (taxon_id,name,assembly,genebuild) "
        ."VALUES (  ".$core_dbs->{$sp}{'taxon_id'}
              .",".$dbh->quote($core_dbs->{$sp}{'name'})
              .",".$dbh->quote($core_dbs->{$sp}{'assembly'})
              .",".$dbh->quote($core_dbs->{$sp}{'genebuild'})
              .")");
  }
  $sth->execute();
  if ($sth->rows > 0){
    my $ref = $sth->fetchrow_arrayref();
    $core_dbs->{$sp}{'genome_db_id'} = $ref->[0];
    $core_dbs->{$sp}{'taxon_id'} = $ref->[1];
    $core_dbs->{$sp}{'assembly'} = $ref->[2];
    $core_dbs->{$sp}{'genebuild'} = $ref->[3];
    add_to_species_set($dbh,$core_dbs,$sp);
    return $core_dbs->{$sp}{'genome_db_id'};
  }
  return;

}

sub add_to_species_set {
  my ($dbh,$core_dbs,$sp) = @_;
  my $msth = $dbh->prepare("SELECT max(species_set_id) FROM species_set");
  $msth->execute();
  my $ss_id;
  if ($msth->rows > 0){
    $ss_id = $msth->fetchrow_arrayref()->[0];
  }
  if (!$ss_id || $ss_id < 1) {
    $ss_id = 1;
    $dbh->do("INSERT INTO species_set_tag (species_set_id,tag,value) "
        ."VALUES (  ".1
              .",'name'"
              .",'default'"
              .")");
    $dbh->do("INSERT INTO species_set_header (species_set_id,name,size) "
        ."VALUES (  ".1
              .",'default'"
              .",0"
              .")");

  }
  $dbh->do("INSERT INTO species_set (species_set_id,genome_db_id) "
        ."VALUES (  ".1
              .",".$core_dbs->{$sp}{'genome_db_id'}
              .")");
  $dbh->do("UPDATE species_set_header SET size = size + 1 WHERE species_set_id = 1");
  $ss_id++;
  # add single sequence and pairwise sequence sets
  my $sth = $dbh->prepare("SELECT genome_db_id,name FROM genome_db WHERE genome_db_id != ".$core_dbs->{$sp}{'genome_db_id'});
  $sth->execute();
  while (my $ref = $sth->fetchrow_arrayref()){
    my $db_id = $ref->[0];
    my $db_name = $ref->[1];
    $dbh->do("INSERT INTO species_set (species_set_id,genome_db_id) "
        ."VALUES (  ".$ss_id
              .",".$core_dbs->{$sp}{'genome_db_id'}
              .")");
    $dbh->do("INSERT INTO species_set (species_set_id,genome_db_id) "
        ."VALUES (  ".$ss_id
              .",".$db_id
              .")");
    $dbh->do("INSERT INTO species_set_tag (species_set_id,tag,value) "
        ."VALUES (  ".$ss_id
              .",'name'"
              .",".$dbh->quote(join("-",sort(($core_dbs->{$sp}{'name'},$db_name))))
              .")");
    $dbh->do("INSERT INTO species_set_header (species_set_id,name,size) "
        ."VALUES (  ".$ss_id
              .",".$dbh->quote(join("-",sort(($core_dbs->{$sp}{'name'},$db_name))))
              .",2"
              .")");
    $ss_id++;
  }
  $dbh->do("INSERT INTO species_set (species_set_id,genome_db_id) "
        ."VALUES (  ".$ss_id
              .",".$core_dbs->{$sp}{'genome_db_id'}
              .")");
  $dbh->do("INSERT INTO species_set_tag (species_set_id,tag,value) "
        ."VALUES (  ".$ss_id
              .",'name'"
              .",".$dbh->quote($core_dbs->{$sp}{'name'})
              .")");
  $dbh->do("INSERT INTO species_set_header (species_set_id,name,size) "
        ."VALUES (  ".$ss_id
              .",".$dbh->quote($core_dbs->{$sp}{'name'})
              .",1"
              .")");
  return $ss_id;

}

# create the compara database from a template
# populate meta, ncbi_taxa_node, ncbi_taxa_tree tables
sub setup_compara_db {
  my $params = shift;
  my $dbh = compara_db_connect($params);
  $dbh->do('DROP DATABASE IF EXISTS '.$params->{'DATABASE_COMPARA'}{'NAME'}) || die "ERROR: unable to drop existing [DATABASE_COMPARA] NAME ".$params->{'DATABASE_COMPARA'}{'NAME'}." using provided settings";
  $dbh->do('CREATE DATABASE '.$params->{'DATABASE_COMPARA'}{'NAME'}) || die "ERROR: unable to create to [DATABASE_COMPARA] NAME ".$params->{'DATABASE_COMPARA'}{'NAME'}." using provided settings";

  my $file = $params->{'ENSEMBL'}{'LOCAL'}.'/ensembl-compara/sql/table.sql';
  my $connection_info =    '-u'.$params->{'DATABASE_COMPARA'}{'RW_USER'}
              .' -h'.$params->{'DATABASE_COMPARA'}{'HOST'}
              .' -P'.$params->{'DATABASE_COMPARA'}{'PORT'}
              .' -D'.$params->{'DATABASE_COMPARA'}{'NAME'}
              .' -p'.$params->{'DATABASE_COMPARA'}{'RW_PASS'};
  system "mysql $connection_info < $file";

  $dbh->do('USE '.$params->{'DATABASE_COMPARA'}{'NAME'});
  ## add method_link, ncbi_taxa_name and ncbi_taxa_node tables

  system "wget ".$params->{'DATABASE_TEMPLATE'}{'URL'}."/".$params->{'DATABASE_TEMPLATE'}{'NAME'}."/method_link.txt.gz";
  system "wget ".$params->{'DATABASE_TEMPLATE'}{'URL'}."/".$params->{'DATABASE_TEMPLATE'}{'NAME'}."/ncbi_taxa_name.txt.gz";
  system "wget ".$params->{'DATABASE_TEMPLATE'}{'URL'}."/".$params->{'DATABASE_TEMPLATE'}{'NAME'}."/ncbi_taxa_node.txt.gz";
  system "gunzip *.txt.gz";
  system 'mysqlimport --fields_escaped_by=\\\\ '
      .' -h '.$params->{'DATABASE_COMPARA'}{'HOST'}
      .' -P '.$params->{'DATABASE_COMPARA'}{'PORT'}
      .' -u '.$params->{'DATABASE_COMPARA'}{'RW_USER'}
      .' -p'.$params->{'DATABASE_COMPARA'}{'RW_PASS'}
      .' '.$params->{'DATABASE_COMPARA'}{'NAME'}
      .' -L method_link.txt ncbi_taxa_name.txt ncbi_taxa_node.txt';
  system "rm method_link.txt ncbi_taxa_name.txt ncbi_taxa_node.txt";

  #initialise gene_tree_root and gene_tree_node
  $dbh->do("INSERT INTO gene_tree_root (root_id,member_type,tree_type,clusterset_id,method_link_species_set_id,gene_align_id,ref_root_id,stable_id,version) "
      ."VALUES (1,'protein','clusterset','default',1,NULL,NULL,NULL,NULL)");
  $dbh->do("INSERT INTO gene_tree_node (node_id,parent_id,root_id,left_index,right_index,distance_to_parent,seq_member_id) "
      ."VALUES (1,NULL,1,0,0,0,NULL)");
  #initialise species_tree_root and species_tree_node
  my $species_newick = "";
  open SPECIES_NEWICK, $params->{'SPECIES_SET'}{'TREE_FILE'};
  while (<SPECIES_NEWICK>) {
    chomp; $species_newick .= $_;
  }
  $dbh->do("INSERT INTO species_tree_root (root_id,method_link_species_set_id,label)"
      ."VALUES (1000,1,".$dbh->quote($params->{'SPECIES_SET'}{'TREE_LABEL'}).")");
  return $dbh;
}

sub compara_db_connect {
  my $params = shift;
  my $dsn = "DBI:mysql:host=$params->{'DATABASE_COMPARA'}{'HOST'};port=$params->{'DATABASE_COMPARA'}{'PORT'}";
  my $dbh = DBI->connect($dsn,"$params->{'DATABASE_COMPARA'}{'RW_USER'}","$params->{'DATABASE_COMPARA'}{'RW_PASS'}") || die "ERROR: unable to connect to [DATABASE_COMPARA] HOST ".$params->{'DATABASE_COMPARA'}{'HOST'}." using provided settings";
  $dbh->do('CREATE DATABASE IF NOT EXISTS '.$params->{'DATABASE_COMPARA'}{'NAME'}) || die "ERROR: unable to create [DATABASE_COMPARA] NAME ".$params->{'DATABASE_COMPARA'}{'NAME'}." using provided settings";
  $dbh->do('USE '.$params->{'DATABASE_COMPARA'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_COMPARA] NAME ".$params->{'DATABASE_COMPARA'}{'NAME'}." using provided settings";
  return $dbh;
}

#sub template_db_connect {
# my $params = shift;
# my $dsn = "DBI:mysql:host=$params->{'DATABASE_TEMPLATE'}{'HOST'};port=$params->{'DATABASE_TEMPLATE'}{'PORT'}";
# my $dbh = DBI->connect($dsn,"$params->{'DATABASE_TEMPLATE'}{'RO_USER'}") || die "ERROR: unable to connect to [DATABASE_TEMPLATE] HOST ".$params->{'DATABASE_TEMPLATE'}{'HOST'}." using provided settings";
# $dbh->do('USE '.$params->{'DATABASE_TEMPLATE'}{'NAME'}) || die "ERROR: unable to connect to [DATABASE_TEMPLATE] NAME ".$params->{'DATABASE_TEMPLATE'}{'NAME'}." using provided settings";
# return $dbh;
#}
#
sub core_db_host_connect {
  my $params = shift;
  my $dsn = "DBI:mysql:host=$params->{'DATABASE_CORE'}{'HOST'};port=$params->{'DATABASE_CORE'}{'PORT'}";
  my $dbh = DBI->connect($dsn,"$params->{'DATABASE_CORE'}{'RO_USER'}") || die "ERROR: unable to connect to [DATABASE_CORE] HOST ".$params->{'DATABASE_CORE'}{'HOST'}." using provided settings";
  return $dbh;
}

sub core_db_connect {
  my ($db_name,$params) = @_;
  my $dbh = core_db_host_connect($params);
  $dbh->do('USE '.$db_name) || die "ERROR: unable to connect to ".$db_name." using provided settings";
  return $dbh;
}

sub make_orthogroup_files {
  my $params = shift @_;
  my @species_list   = keys %{$params->{'TAXA'}};
  my $species_regexp = join("|",@species_list);

  my %sequences;

  # load all sequences in memory
  for my $species (@species_list) {
    $sequences{$species}{faa} = _fasta_file_to_hash($params->{'SETUP'}{'FASTA_DIR'} . "/$species\_-_canonical_proteins.fa");
    $sequences{$species}{fna} = _fasta_file_to_hash($params->{'SETUP'}{'FASTA_DIR'} . "/$species\_-_cds_translationid.fa");
    $sequences{$species}{fba} = _fasta_file_to_hash($params->{'SETUP'}{'FASTA_DIR'} . "/$species\_-_protein_bounded_exon.fa");
  }

  # create folder with sequences for each orthogroup
  open OG, "<" . $params->{'SETUP'}{'ORTHOFINDER_GROUPS'} or die $!;
  mkdir "orthogroups";
  while (<OG>) {
    my @tokens = split /\s+/;
    my $orthogroup_id = shift @tokens;
    $orthogroup_id =~ s/:$//;
    mkdir "orthogroups/$orthogroup_id";
    open IDS, ">orthogroups/$orthogroup_id/$orthogroup_id" or die $!;
    open FAA, ">orthogroups/$orthogroup_id/$orthogroup_id.faa" or die $!;
    open FNA, ">orthogroups/$orthogroup_id/$orthogroup_id.fna" or die $!;
    open FBA, ">orthogroups/$orthogroup_id/$orthogroup_id.fba" or die $!;
    for my $sequence_id (@tokens) {
      if ($sequence_id =~ /^(${species_regexp})_\S+$/) {
        my $species = $1;
        if (exists $sequences{$species}{faa}{$sequence_id} and 
            exists $sequences{$species}{fna}{$sequence_id} and
            exists $sequences{$species}{fba}{$sequence_id}) {
          print IDS "$sequence_id\n";
          print FAA ">$sequence_id\n" . $sequences{$species}{faa}{$sequence_id}->{seq} . "\n";
          print FNA ">$sequence_id\n" . $sequences{$species}{fna}{$sequence_id}->{seq} . "\n";
          print FBA ">$sequence_id\n" . $sequences{$species}{fba}{$sequence_id}->{seq} . "\n";
        }
        else {
          warn "$sequence_id does not exist in one of the input fasta files";
        }
      }
      else {
        warn "$sequence_id does not match species list";
      }
    }
    close IDS;
    close FAA;
    close FNA;
    close FBA;
  }
}

1;
