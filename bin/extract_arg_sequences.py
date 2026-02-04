#!/usr/bin/env python3
"""
Extract ARG sequences from annotation files based on AMRFinderPlus results.

This script takes AMRFinderPlus TSV output and extracts the corresponding
protein sequences from annotation files (e.g., Pyrodigal FAA output).
"""

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, Set, TextIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_amrfinder_results(amrfinder_file: Path) -> Dict[str, dict]:
    """
    Parse AMRFinderPlus TSV output and extract protein IDs with metadata.
    
    Args:
        amrfinder_file: Path to AMRFinderPlus TSV output
        
    Returns:
        Dictionary mapping protein IDs to their ARG metadata
    """
    arg_info = {}
    
    opener = gzip.open if amrfinder_file.suffix == '.gz' else open
    
    with opener(amrfinder_file, 'rt') as f:
        # Skip header line
        header = f.readline().strip().split('\t')
        
        # Find relevant column indices
        try:
            protein_id_idx = header.index('Protein identifier')
            gene_symbol_idx = header.index('Gene symbol')
            seq_name_idx = header.index('Sequence name')
            scope_idx = header.index('Scope')
            element_type_idx = header.index('Element type')
            element_subtype_idx = header.index('Element subtype')
            class_idx = header.index('Class')
            subclass_idx = header.index('Subclass')
        except ValueError as e:
            print(f"Error: Could not find expected column in AMRFinderPlus output: {e}", file=sys.stderr)
            sys.exit(1)
        
        # Parse data lines
        for line in f:
            if not line.strip():
                continue
                
            fields = line.strip().split('\t')
            
            if len(fields) <= max(protein_id_idx, gene_symbol_idx):
                continue
                
            protein_id = fields[protein_id_idx]
            
            # Store metadata for this ARG
            arg_info[protein_id] = {
                'gene_symbol': fields[gene_symbol_idx] if gene_symbol_idx < len(fields) else 'Unknown',
                'sequence_name': fields[seq_name_idx] if seq_name_idx < len(fields) else '',
                'scope': fields[scope_idx] if scope_idx < len(fields) else '',
                'element_type': fields[element_type_idx] if element_type_idx < len(fields) else '',
                'element_subtype': fields[element_subtype_idx] if element_subtype_idx < len(fields) else '',
                'class': fields[class_idx] if class_idx < len(fields) else '',
                'subclass': fields[subclass_idx] if subclass_idx < len(fields) else ''
            }
    
    return arg_info


def extract_sequences(annotation_file: Path, arg_ids: Set[str], arg_info: Dict[str, dict]) -> list:
    """
    Extract sequences from annotation file for given protein IDs.
    
    Args:
        annotation_file: Path to annotation FAA file (can be gzipped)
        arg_ids: Set of protein IDs to extract
        arg_info: Dictionary with metadata for each ARG
        
    Returns:
        List of SeqRecord objects for extracted ARGs
    """
    extracted_sequences = []
    opener = gzip.open if annotation_file.suffix == '.gz' else open
    
    with opener(annotation_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Check if this protein ID matches any ARG
            # Handle different ID formats from different annotation tools
            record_id = record.id.split()[0]  # Get the main ID without description
            
            if record_id in arg_ids:
                # Add metadata to the description
                metadata = arg_info[record_id]
                gene_symbol = metadata['gene_symbol']
                arg_class = metadata['class']
                
                # Update the record description
                record.description = f"{record.description} | ARG:{gene_symbol} | Class:{arg_class}"
                extracted_sequences.append(record)
    
    return extracted_sequences


def write_summary(output_file: Path, arg_info: Dict[str, dict], found_ids: Set[str]):
    """
    Write a summary TSV of extracted ARGs.
    
    Args:
        output_file: Path to output summary TSV
        arg_info: Dictionary with metadata for each ARG
        found_ids: Set of protein IDs that were found in annotations
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join([
            'Protein_ID',
            'Gene_Symbol',
            'Sequence_Name',
            'Scope',
            'Element_Type',
            'Element_Subtype',
            'Class',
            'Subclass',
            'Found_In_Annotation'
        ]) + '\n')
        
        # Write data
        for protein_id, metadata in arg_info.items():
            found = 'Yes' if protein_id in found_ids else 'No'
            f.write('\t'.join([
                protein_id,
                metadata['gene_symbol'],
                metadata['sequence_name'],
                metadata['scope'],
                metadata['element_type'],
                metadata['element_subtype'],
                metadata['class'],
                metadata['subclass'],
                found
            ]) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract ARG sequences from annotation files based on AMRFinderPlus results'
    )
    parser.add_argument(
        '--amrfinder',
        type=Path,
        required=True,
        help='AMRFinderPlus TSV output file'
    )
    parser.add_argument(
        '--annotation',
        type=Path,
        required=True,
        help='Annotation FAA file (can be gzipped)'
    )
    parser.add_argument(
        '--output_fasta',
        type=Path,
        required=True,
        help='Output FASTA file for ARG sequences'
    )
    parser.add_argument(
        '--output_summary',
        type=Path,
        required=True,
        help='Output TSV summary file'
    )
    
    args = parser.parse_args()
    
    # Check input files exist
    if not args.amrfinder.exists():
        print(f"Error: AMRFinderPlus file not found: {args.amrfinder}", file=sys.stderr)
        sys.exit(1)
    
    if not args.annotation.exists():
        print(f"Error: Annotation file not found: {args.annotation}", file=sys.stderr)
        sys.exit(1)
    
    # Parse AMRFinderPlus results
    print(f"Parsing AMRFinderPlus results from {args.amrfinder}...")
    arg_info = parse_amrfinder_results(args.amrfinder)
    print(f"Found {len(arg_info)} ARG hits in AMRFinderPlus output")
    
    if not arg_info:
        print("Warning: No ARGs found in AMRFinderPlus output", file=sys.stderr)
        # Create empty output files
        args.output_fasta.touch()
        with open(args.output_summary, 'w') as f:
            f.write('Protein_ID\tGene_Symbol\tSequence_Name\tScope\tElement_Type\tElement_Subtype\tClass\tSubclass\tFound_In_Annotation\n')
        return
    
    # Extract sequences
    print(f"Extracting sequences from {args.annotation}...")
    arg_ids = set(arg_info.keys())
    extracted = extract_sequences(args.annotation, arg_ids, arg_info)
    
    print(f"Extracted {len(extracted)} ARG sequences")
    
    # Write output FASTA
    with open(args.output_fasta, 'w') as f:
        SeqIO.write(extracted, f, 'fasta')
    
    # Write summary
    found_ids = {record.id for record in extracted}
    write_summary(args.output_summary, arg_info, found_ids)
    
    print(f"Output written to:")
    print(f"  FASTA: {args.output_fasta}")
    print(f"  Summary: {args.output_summary}")
    
    # Report any missing sequences
    missing = arg_ids - found_ids
    if missing:
        print(f"\nWarning: {len(missing)} ARG(s) from AMRFinderPlus not found in annotation file:")
        for missing_id in sorted(missing):
            print(f"  - {missing_id} ({arg_info[missing_id]['gene_symbol']})")


if __name__ == '__main__':
    main()
