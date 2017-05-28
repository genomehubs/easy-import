# Change values below only if you are using your own
# cluster infrastructure to run the genetree pipeline

# cd to the folder where you have the orthogroup protein, cds, and bounded exon fasta files:
cd orthogroups/$ORTHOGROUPID

# set paths:
MAFFT=/usr/local/bin/mafft
NOISY=/usr/local/bin/noisy
RAXML=/usr/local/bin/raxmlHPC-PTHREADS-SSE3
NOTUNG="java -jar /Notung-2.9/Notung-2.9.jar"

# set the species tree newick file to reconcile the gene tree with:
NOTUNG_SPECIESTREE=/import/data/speciestree.newick

$MAFFT --treeout --auto --reorder $ORTHOGROUPID.faa > $ORTHOGROUPID.faa.mafft && \
$NOISY --seqtype P $ORTHOGROUPID.faa.mafft && \
$RAXML -f a -x 12345 -# 100 -T 1 -p 12345 -m PROTGAMMAAUTO -s $ORTHOGROUPID.faa_out.fas -n $ORTHOGROUPID && \
rename 's/(RAxML_\S+?)\.(\S+)/$2.$1/' RAxML_*.$ORTHOGROUPID && \
$NOTUNG --treeoutput nhx --root -s $NOTUNG_SPECIESTREE -g $$ORTHOGROUPID.RAxML_bipartitionsBranchLabels && \
$NOTUNG --treeoutput nhx --root -s $NOTUNG_SPECIESTREE -g $$ORTHOGROUPID.RAxML_bipartitionsBranchLabels.rooting.0 \
  --nolosses --treeoutput nhx --homologtabletabs --reconcile --stpruned && \
touch $ORTHOGROUPID.done
