#!/usr/bin/env python3

import sys

def mock_bprom_output():
    return """
Promoter Prediction Output
-10 box at position 35
-35 box at position 15
TATA box at position 25
TSS at position 45
    """

if __name__ == "__main__":
    input_file = sys.argv[1]
    print(mock_bprom_output())

