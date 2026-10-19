import pyBigWig
import pandas as pd
import argparse
import os

def calculate_pause_index(bw_paths, bed_path, output_file, sample_names=None):
    # Load gene annotations
    genes = pd.read_csv(bed_path, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'name', 'score', 'strand',
                               'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                               'blockSizes', 'blockStarts'])
    
    # Generate default sample names from filenames
    if sample_names is None:
        sample_names = [os.path.splitext(os.path.basename(path))[0] for path in bw_paths]
    elif len(sample_names) != len(bw_paths):
        raise ValueError("Number of sample names must match number of bigWig files")
    
    # Initialize gene dictionary to collect all data
    genes_dict = {}
    for _, gene in genes.iterrows():
        genes_dict[gene['name']] = {
            'gene': gene['name'],
            'chr': gene['chr'],
            'start': gene['start'],
            'end': gene['end'],
            'strand': gene['strand'],
            'gene_length': gene['end'] - gene['start']
        }
    
    # Process each bigWig file
    for bw_path, sample_name in zip(bw_paths, sample_names):
        bw = pyBigWig.open(bw_path)
        for gene_name in genes_dict:
            gene_data = genes_dict[gene_name]
            chrom = gene_data['chr']
            strand = gene_data['strand']
            
            # Determine TSS and TES based on strand
            if strand == '+':
                tss = gene_data['start']
                tes = gene_data['end']
            else:
                tss = gene_data['end']
                tes = gene_data['start']
            
            # Define regions
            if strand == '+':
                tss_region = (tss, tss + 250)
                gene_body_start = tss + 250
                gene_body_end = tes
            else:
                tss_region = (max(0, tss - 250), tss)
                gene_body_start = max(0, tss - 250)
                gene_body_end = tes + 250
            
            try:
                # Get TSS signal
                tss_signal = bw.stats(chrom, tss_region[0], tss_region[1], exact=True)[0] or 0
                
                # Get gene body signal
                if gene_body_end > gene_body_start:
                    gene_body_signal = bw.stats(chrom, gene_body_start, gene_body_end, exact=True)[0] or 0
                else:
                    gene_body_signal = 0
                
                # Calculate pause index
                pause_index = (tss_signal + 1) / (gene_body_signal + 1)
                
            except RuntimeError as e:
                print(f"Error processing {gene_name} in {sample_name}: {str(e)}")
                tss_signal = float('nan')
                gene_body_signal = float('nan')
                pause_index = float('nan')
            
            # Store results with sample-specific column names
            gene_data[f"{sample_name}_tss"] = tss_signal
            gene_data[f"{sample_name}_body"] = gene_body_signal
            gene_data[f"{sample_name}_pause"] = pause_index
        
        bw.close()
    
    # Convert to DataFrame and save
    result_df = pd.DataFrame(genes_dict.values())
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Calculate pause indices for multiple samples')
    parser.add_argument('bed_file', help='Gene annotations BED file')
    parser.add_argument('output_file', help='Output TSV file')
    parser.add_argument('bw_files', nargs='+', help='bigWig files')
    parser.add_argument('--sample_names', nargs='+', help='Optional custom sample names')
    
    args = parser.parse_args()
    
    if args.sample_names and len(args.sample_names) != len(args.bw_files):
        raise ValueError("Mismatch between sample names and bigWig files")
    
    calculate_pause_index(
        bw_paths=args.bw_files,
        bed_path=args.bed_file,
        output_file=args.output_file,
        sample_names=args.sample_names
    )

if __name__ == "__main__":
    main()