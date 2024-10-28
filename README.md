# Bacterial Promoter Prediction Pipeline

A Python-based pipeline for predicting bacterial promoter elements in genomic sequences using BPROM (or a mock version for testing).

## Overview

This pipeline extracts upstream regions from annotated genes and predicts promoter elements including:
- TATA box
- -35 box
- -10 box
- Transcription Start Site (TSS)

## Installation

1. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # Unix/macOS
venv\Scripts\activate     # Windows
```

2. Install required packages:
```bash
pip install biopython bcbio-gff
```

## Usage

1. Generate test data:
```bash
python test_data.py
```

2. Run the pipeline:
```bash
python main.py
```

## File Structure

- `main.py`: Core pipeline implementation
- `mock_bprom.py`: Mock BPROM predictor for testing
- `test_data.py`: Creates test genome and annotation files

## Implementation Details

### Test Data Generation
```python
# Create a small test genome
test_genome = SeqRecord(
    Seq("ATGCATGCATGCATGCATGC" * 200),  # 4000bp sequence
    id="test_genome",
    name="test_genome"
)
```
Creates a 4000bp test genome with a repeating pattern.

### Upstream Region Extraction
```python
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
```
Extracts upstream regions considering:
- Strand orientation
- Genome boundaries
- Configurable length

### Promoter Prediction
```python
def mock_bprom_output():
    return """
Promoter Prediction Output
-10 box at position 35
-35 box at position 15
TATA box at position 25
TSS at position 45
    """
```
Runs BPROM (or mock version) on upstream sequences to predict promoter elements.

## Output Format

The pipeline generates a promoter bank with entries containing:
- Promoter ID
- Sequence information
- Associated gene details
- Predicted promoter elements
- Sequence length

Example output:
```
Processed promoter 0:
  Sequence length: 1500
  Gene: TEST1
  Predictions: {'TATA_box': [25], 'minus35': [15], 'minus10': [35], 'TSS': [45]}
```
