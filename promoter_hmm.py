#!/usr/bin/env python3

'''Script for scoring potential promoter sequences.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states (can keep as default)

For prokaryotes
python3 promoter_hmm.py -f promoter.fasta --organism prokaryote

For eukaryotes (default)
python3 promoter_hmm.py -f promoter.fasta
or
python3 promoter_hmm.py -f promoter.fasta --organism eukaryote
'''

import argparse
import numpy as np
import math

def sumLogProbs(a,b):
    if a > b: 
        return a + np.log1p(math.exp(b-a))
    else:
        return b + np.log1p(math.exp(a-b))

def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s

def forward_backward(obs, trans_probs, emiss_probs, init_probs):
    states = list(init_probs.keys())
    N = len(obs)
    F = {state: np.full(N, -np.inf) for state in states}
    B = {state: np.full(N, -np.inf) for state in states}
    R = {state: np.full(N, 0.0) for state in states}

    # Forward algorithm
    for state in states:
        F[state][0] = init_probs[state] + emiss_probs[state][obs[0]]

    for t in range(1, N):
        for curr_state in states:
            sum_log_probs = -np.inf 
            for prev_state in states:
                log_prob = F[prev_state][t - 1] + trans_probs[prev_state][curr_state]
                sum_log_probs = sumLogProbs(sum_log_probs, log_prob)
            F[curr_state][t] = sum_log_probs + emiss_probs[curr_state][obs[t]]

    like_f = -np.inf
    for state in states:
        like_f = sumLogProbs(like_f, F[state][N - 1])

    # Backward algorithm
    for state in states:
        B[state][N - 1] = 0.0

    for t in range(N - 2, -1, -1):
        for curr_state in states:
            sum_log_probs = -np.inf
            for next_state in states:
                log_prob = (trans_probs[curr_state][next_state] +
                          emiss_probs[next_state][obs[t + 1]] +
                          B[next_state][t + 1])
                sum_log_probs = sumLogProbs(sum_log_probs, log_prob)
            B[curr_state][t] = sum_log_probs

    like_b = -np.inf
    for state in states:
        log_prob = (init_probs[state] +
                   emiss_probs[state][obs[0]] +
                   B[state][0])
        like_b = sumLogProbs(like_b, log_prob)

    # Calculate posteriors
    for t in range(N):
        px_t = -np.inf
        for state in states:
            log_prob = F[state][t] + B[state][t]
            px_t = sumLogProbs(px_t, log_prob)

        for state in states:
            log_posterior = F[state][t] + B[state][t] - px_t
            R[state][t] = np.exp(log_posterior)

    return F, like_f, B, like_b, R

def score_promoter(sequence, organism_type='eukaryote'):
    """Score promoter sequences with organism-specific parameters
    Args:
        sequence: DNA sequence string
        organism_type: either 'prokaryote' or 'eukaryote'
    """
    # Define states for promoter elements
    mu = 0.1  # Probability of switching states

    # Define transition probabilities based on organism type
    transition_probabilities = {
        'TATA': {
            'TATA': np.log(0.7 if organism_type == 'prokaryote' else 0.6),
            'GC': np.log(0.1),
            'CAAT': np.log(0.1),
            'BG': np.log(0.1 if organism_type == 'prokaryote' else 0.2)
        },
        'GC': {
            'TATA': np.log(0.1),
            'GC': np.log(0.7 if organism_type == 'prokaryote' else 0.6),
            'CAAT': np.log(0.1),
            'BG': np.log(0.1 if organism_type == 'prokaryote' else 0.2)
        },
        'CAAT': {
            'TATA': np.log(0.1),
            'GC': np.log(0.1),
            'CAAT': np.log(0.7 if organism_type == 'prokaryote' else 0.6),
            'BG': np.log(0.1 if organism_type == 'prokaryote' else 0.2)
        },
        'BG': {
            'TATA': np.log(0.1),
            'GC': np.log(0.1),
            'CAAT': np.log(0.1),
            'BG': np.log(0.7 if organism_type == 'prokaryote' else 0.6)
        }
    }

    # Define emission probabilities based on organism type
    emission_probabilities = {
        'TATA': {
            'A': np.log(0.45 if organism_type == 'prokaryote' else 0.35),
            'T': np.log(0.45 if organism_type == 'prokaryote' else 0.35),
            'G': np.log(0.05 if organism_type == 'prokaryote' else 0.15),
            'C': np.log(0.05 if organism_type == 'prokaryote' else 0.15)
        },
        'GC': {
            'G': np.log(0.45 if organism_type == 'prokaryote' else 0.4),
            'C': np.log(0.45 if organism_type == 'prokaryote' else 0.4),
            'A': np.log(0.05 if organism_type == 'prokaryote' else 0.1),
            'T': np.log(0.05 if organism_type == 'prokaryote' else 0.1)
        },
        'CAAT': {
            'C': np.log(0.35 if organism_type == 'prokaryote' else 0.3),
            'A': np.log(0.35 if organism_type == 'prokaryote' else 0.3),
            'T': np.log(0.2),
            'G': np.log(0.1 if organism_type == 'prokaryote' else 0.2)
        },
        'BG': {
            'A': np.log(0.25),
            'T': np.log(0.25),
            'G': np.log(0.25),
            'C': np.log(0.25)
        }
    }

    # Equal initial probabilities
    initial_probabilities = {
        'TATA': np.log(0.25),
        'GC': np.log(0.25),
        'CAAT': np.log(0.25),
        'BG': np.log(0.25)
    }

    # Run forward-backward algorithm
    F, like_f, B, like_b, R = forward_backward(sequence, 
                                              transition_probabilities,
                                              emission_probabilities,
                                              initial_probabilities)
    
    # Find elements and their positions
    elements_found = []
    for state in ['TATA', 'GC', 'CAAT']:
        # Look for stretches where state probability > 0.5
        prob_threshold = 0.5
        in_element = False
        start_pos = 0
        
        for pos in range(len(sequence)):
            if R[state][pos] > prob_threshold and not in_element:
                start_pos = pos
                in_element = True
            elif (R[state][pos] <= prob_threshold or pos == len(sequence)-1) and in_element:
                elements_found.append({
                    "type": state,
                    "position": (start_pos, pos),
                    "confidence": float(np.mean([R[state][i] for i in range(start_pos, pos)]))
                })
                in_element = False

    # Calculate overall promoter confidence
    # Average of highest probabilities for each promoter element
    promoter_confidence = np.mean([max(R[state]) for state in ['TATA', 'GC', 'CAAT']])

    return {
        "confidence_score": float(promoter_confidence),
        "elements_found": elements_found
    }

def main():
    parser = argparse.ArgumentParser(
        description='Score potential promoter sequences.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True,
                      help='File containing the sequence (fasta file)')
    parser.add_argument('--organism', type=str, choices=['prokaryote', 'eukaryote'],
                      default='eukaryote',
                      help='Organism type (default: eukaryote)')
    
    args = parser.parse_args()
    sequence = read_fasta(args.f)
    
    result = score_promoter(sequence, organism_type=args.organism)
    
    print(f"\nAnalyzing {args.organism} sequence")
    print(f"Promoter Confidence Score: {result['confidence_score']:.2f}")
    print("\nElements Found:")
    for element in result["elements_found"]:
        print(f"- {element['type']} at positions {element['position']} "
              f"(confidence: {element['confidence']:.2f})")

if __name__ == "__main__":
    main()