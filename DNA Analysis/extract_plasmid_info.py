import os
import sys
import json
import argparse
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from snapgene_reader import snapgene_file_to_seqrecord
from Bio.Seq import Seq

# Import the find_orfs function from the existing script
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from snapgene_renderer import find_orfs

def get_best_orf(record, orfs):
    """Find the ORF that contains one of the tags, or fallback to the longest ORF."""
    seq_len = len(record.seq)
    doubled = str(record.seq).upper() * 2
    
    for (start, end, strand, n_aa) in orfs:
        if start < 0:
            dna_seq = doubled[start + seq_len : end + seq_len]
        else:
            dna_seq = doubled[start : end]
        if strand == -1:
            dna_seq = str(Seq(dna_seq).reverse_complement())
            
        protein_seq = str(Seq(dna_seq).translate()).replace('*', '')
        
        # Check if any tag is in the sequence (His, Myc, FLAG)
        if 'HHHHHH' in protein_seq or 'EQKLISEEDL' in protein_seq or 'DYKDDDDK' in protein_seq:
            return (start, end, strand, n_aa)
            
    # Fallback to longest ORF
    if orfs:
        return max(orfs, key=lambda x: x[3])
    return None

def extract_metadata(dna_path):
    record = snapgene_file_to_seqrecord(dna_path)
    seq_len = len(record.seq)
    
    # 1. Tags & Antibiotic Resistance & Promoters
    tags = []
    abx = []
    promoters = []
    
    for f in record.features:
        label = f.qualifiers.get('label', f.qualifiers.get('name', f.qualifiers.get('note', [''])))
        name = label[0] if isinstance(label, list) else label
        name_lower = name.lower()
        
        if 'promoter' in f.type.lower() or 'promoter' in name_lower:
            if name not in promoters:
                promoters.append(name)
        
        if f.type == 'CDS' and name.endswith('R'):
            if name not in abx:
                abx.append(name)
                
        # Common tags
        if any(tag in name_lower for tag in ['his', 'myc', 'flag', 'ha']):
            if name not in tags:
                tags.append(name)

    # 2. Protein sequence
    orfs = find_orfs(record.seq, min_aa=100)
    best_orf = get_best_orf(record, orfs)
    
    protein_seq = ""
    dna_seq = ""
    mw_kda = 0.0
    pi = 0.0
    
    if best_orf:
        start, end, strand, n_aa = best_orf
        
        # Extract sequence
        doubled = str(record.seq).upper() * 2
        if start < 0:
            dna_seq = doubled[start + seq_len : end + seq_len]
        else:
            dna_seq = doubled[start : end]
            
        if strand == -1:
            dna_seq = str(Seq(dna_seq).reverse_complement())
            
        protein_seq = str(Seq(dna_seq).translate())
        
        # Remove all stop codons (*) for analysis
        protein_seq_clean = protein_seq.replace('*', '')
            
        # Calculate MW and pI
        if protein_seq_clean:
            analysis = ProteinAnalysis(protein_seq_clean)
            mw_kda = analysis.molecular_weight() / 1000.0
            pi = analysis.isoelectric_point()
        
    name = os.path.splitext(os.path.basename(dna_path))[0]
    
    return {
        'name': name,
        'size_kda': round(mw_kda, 2),
        'pi': round(pi, 2),
        'tags': tags,
        'abx': abx,
        'promoters': promoters,
        'dna_seq': dna_seq,
        'protein_seq': protein_seq
    }

def main():
    parser = argparse.ArgumentParser(description="Extract metadata and ORFs from SnapGene (.dna) files.")
    parser.add_argument('dna_files', nargs='+', help="Paths to one or more .dna files to process")
    parser.add_argument('--indent', type=int, default=2, help="Indentation level for JSON output (default: 2)")
    
    args = parser.parse_args()
    
    results = []
    for path in args.dna_files:
        if not os.path.exists(path):
            print(f"Warning: File not found: {path}", file=sys.stderr)
            continue
        try:
            res = extract_metadata(path)
            results.append(res)
        except Exception as e:
            print(f"Error processing {path}: {e}", file=sys.stderr)
            
    # Output generalized format as JSON
    print(json.dumps(results, indent=args.indent))

if __name__ == '__main__':
    main()
