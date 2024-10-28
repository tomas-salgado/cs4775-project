# Create a small test genome and annotations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Create a small test genome
test_genome = SeqRecord(
    Seq("ATGCATGCATGCATGCATGC" * 200),  # 4000bp sequence
    id="test_genome",
    name="test_genome"
)

# Write test genome
with open("test_genome.fasta", "w") as handle:
    SeqIO.write(test_genome, handle, "fasta")

# Create test annotations
test_feature = SeqFeature(
    location=FeatureLocation(2000, 2500),  # Gene in the middle
    type="gene",
    qualifiers={
        "gene_id": ["test_gene_1"],
        "gene": ["TEST1"],
        "strand": 1
    }
)

test_record = SeqRecord(
    Seq(""),
    id="test_genome",
    features=[test_feature]
)

# Write test annotations
with open("test_annotations.gff", "w") as handle:
    GFF.write([test_record], handle)
