import os
import gzip
import json
import pickle
import re
from collections import defaultdict
from vcf import Reader

AA_MAP = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': 'X', '*': 'X'  # Handle stop codons
}

def parse_protein_change(change):
    """Parse HGVS.p notation with enhanced handling of variants"""
    # Expanded regex to handle stop codons and non-standard notations
    match = re.match(r"p\.([A-Za-z]{3})(\d+)([A-Za-z*]{1,3})", change, re.IGNORECASE)
    if not match:
        return None, None
    
    # Normalize amino acid codes
    ref_3 = match.group(1).capitalize()
    alt_3 = match.group(3).capitalize().replace('*', 'Ter')  # Standardize stop codons
    
    # Map to single-letter codes with fallback
    ref_aa = AA_MAP.get(ref_3, 'X')
    alt_aa = AA_MAP.get(alt_3, 'X')
    
    # Skip synonymous and invalid variants
    if ref_aa == alt_aa or 'X' in (ref_aa, alt_aa):
        return None, None
    
    return f"{ref_aa}{match.group(2)}{alt_aa}", match.group(2)

def process_vcf(vcf_path):
    """Process VCF files with robust field handling"""
    sample_id = os.path.basename(vcf_path).split('.')[0]
    results = defaultdict(lambda: {
        'total_variations': 0,
        'variation_positions': defaultdict(int),
        'symbolic_mutations': defaultdict(int)
    })

    try:
        opener = gzip.open if vcf_path.endswith('.gz') else open
        with opener(vcf_path, 'rt') as f:
            vcf_reader = Reader(f)
            
            for record in vcf_reader:
                # Enhanced DP handling
                dp = int(record.INFO.get('DP', [0])[0]) if isinstance(record.INFO.get('DP'), list) else int(record.INFO.get('DP', 0))
                
                # Quality control thresholds
                if record.QUAL < 100 or dp < 8:
                    continue
                
                # Allele frequency calculation
                af = 0.0
                if 'AF' in record.INFO:
                    af_values = record.INFO['AF']
                    af = float(af_values[0] if isinstance(af_values, list) else af_values)
                elif 'AC' in record.INFO and 'AN' in record.INFO:
                    ac = int(record.INFO['AC'][0])
                    an = int(record.INFO['AN'])
                    af = ac / an if an > 0 else 0
                elif 'DP4' in record.INFO:
                    dp4 = list(map(int, record.INFO['DP4']))
                    if len(dp4) >= 4:
                        alt_count = dp4[2] + dp4[3]
                        total = sum(dp4[:4])
                        af = alt_count / total if total > 0 else 0
                
                if af < 0.75:
                    continue
                
                # Genotype quality handling
                gq_passes = True
                if 'GQ' in record.FORMAT:
                    try:
                        gq_values = [sample.data.GQ for sample in record.samples]
                        gq_passes = all(gq >= 20 for gq in gq_values if gq is not None)
                    except AttributeError:
                        gq_passes = False
                
                if gq_passes:
                    process_annotations(record, results)

    except Exception as e:
        print(f"Error processing {vcf_path}: {str(e)}")
    
    return sample_id, dict(results)

def process_annotations(record, results):
    """Process ANN field annotations for a record"""
    for ann in record.INFO.get('ANN', []):
        parts = ann.split('|')
        if len(parts) < 14:
            continue
        
        # Extract relevant fields from ANN
        gene = parts[3]  # Gene name (e.g., 'gyrA')
        protein_change = parts[10]  # Protein change (e.g., 'p.Leu232Pro')
        
        # Parse protein change
        mut_symbol, position = parse_protein_change(protein_change)
        if not mut_symbol or not position:
            continue
        
        # Update results
        results[gene]['total_variations'] += 1
        results[gene]['variation_positions'][position] += 1
        results[gene]['symbolic_mutations'][mut_symbol] += 1

def main(vcf_dir):
    """Process all VCF files in a directory"""
    final_result = {}
    
    for filename in os.listdir(vcf_dir):
        if filename.endswith(('.vcf', '.vcf.gz')):
            vcf_path = os.path.join(vcf_dir, filename)
            sample_id, data = process_vcf(vcf_path)
            if data:
                final_result[sample_id] = data
    
    # Convert defaultdicts to regular dicts
    for sample in final_result:
        for gene in final_result[sample]:
            final_result[sample][gene]['variation_positions'] = dict(
                final_result[sample][gene]['variation_positions'])
            final_result[sample][gene]['symbolic_mutations'] = dict(
                final_result[sample][gene]['symbolic_mutations'])
    
    return final_result

if __name__ == '__main__':
    # Set paths
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    INPUT_DIR = os.path.join(SCRIPT_DIR, "VCFs")
    OUTPUT_DIR = os.path.join(SCRIPT_DIR, "Variations_dictionaries")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Process files
    result = main(INPUT_DIR)
    
    # Save results
    txt_path = os.path.join(OUTPUT_DIR, "snps_results.txt")
    pkl_path = os.path.join(OUTPUT_DIR, "snps_results.pkl")
    
    with open(txt_path, 'w') as f:
        json.dump(result, f, indent=2)
    
    with open(pkl_path, 'wb') as f:
        pickle.dump(result, f)
    
    print(f"Results saved to:\n{txt_path}\n{pkl_path}")