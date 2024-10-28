from Bio import SeqIO
from Bio import SeqFeature
from BCBio import GFF
import subprocess
import json
import os

# Parse genome - use test files
genome = SeqIO.read("test_genome.fasta", "fasta")

# Parse annotations 
with open("test_annotations.gff", "r") as handle:
    annotations = list(GFF.parse(handle))

def extract_upstream(genome, gene, upstream_length=2000):
    """Extract upstream region with better boundary handling and strand consideration"""
    if gene.location.strand == 1:  # Forward strand
        start = max(0, gene.location.start - upstream_length)
        end = gene.location.start
        sequence = genome[start:end]
    else:  # Reverse strand
        start = gene.location.end
        end = min(len(genome), gene.location.end + upstream_length)
        sequence = genome[start:end].reverse_complement()
    return sequence

def parse_bprom_output(output):
    """Parse BPROM output into structured format"""
    parsed_results = {
        "TATA_box": [],
        "minus35": [],
        "minus10": [],
        "TSS": []
    }
    
    for line in output.split('\n'):
        if '-10 box at position' in line:
            pos = int(line.split()[-1])
            parsed_results["minus10"].append(pos)
        elif '-35 box at position' in line:
            pos = int(line.split()[-1])
            parsed_results["minus35"].append(pos)
        elif 'TATA box at position' in line:
            pos = int(line.split()[-1])
            parsed_results["TATA_box"].append(pos)
        elif 'TSS at position' in line:
            pos = int(line.split()[-1])
            parsed_results["TSS"].append(pos)
    
    return parsed_results

def test_workflow():
    """Test the entire workflow"""
    # Load test data
    genome = SeqIO.read("test_genome.fasta", "fasta")
    
    with open("test_annotations.gff", "r") as handle:
        annotations = list(GFF.parse(handle))
    
    # Process one gene
    upstream_regions = []
    gene_info = []
    
    for record in annotations:
        for feature in record.features:
            if feature.type == "gene":
                upstream = extract_upstream(genome, feature)
                upstream_regions.append(upstream)
                
                gene_info.append({
                    "gene_id": feature.qualifiers.get("gene_id", ["Unknown"])[0],
                    "name": feature.qualifiers.get("gene", ["Unknown"])[0],
                    "strand": feature.location.strand,
                    "position": (feature.location.start, feature.location.end)
                })
    
    # Run predictions
    promoter_predictions = [run_bprom(seq) for seq in upstream_regions]
    
    # Create promoter bank
    promoter_bank = []
    for i, (seq, pred, gene) in enumerate(zip(upstream_regions, promoter_predictions, gene_info)):
        entry = {
            "id": f"promoter_{i}",
            "sequence": str(seq),
            "gene_info": gene,
            "predictions": pred,
            "length": len(seq)
        }
        promoter_bank.append(entry)
        print(f"Processed promoter {i}:")
        print(f"  Sequence length: {len(seq)}")
        print(f"  Gene: {gene['name']}")
        print(f"  Predictions: {pred}")
    
    return promoter_bank

def run_bprom(sequence):
    with open("temp.fa", "w") as f:
        f.write(f">seq\n{sequence}")
    
    # Use full path to mock_bprom.py
    script_path = os.path.join(os.getcwd(), "mock_bprom.py")
    result = subprocess.run([script_path, "temp.fa"], capture_output=True, text=True)
    return parse_bprom_output(result.stdout)

if __name__ == "__main__":
    results = test_workflow()
