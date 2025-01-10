import sys
from Bio import SeqIO

def parse_attributes(attributes_str):
    """Parse the attributes column of a GFF3 file and return a dictionary."""
    attributes = {}
    for attribute in attributes_str.split(";"):
        if "=" in attribute:
            key, value = attribute.split("=", 1)
            attributes[key.strip()] = value.strip()
    return attributes

def check_gff3_with_fasta(gff3_file, fasta_file):
    try:
        # Read sequence IDs and lengths from the FASTA file
        fasta_seqs = {record.id: len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
        
        has_gene_feature = False  # To check if there is at least one "gene" feature
        
        with open(gff3_file, 'r') as f:
            for line_number, line in enumerate(f, start=1):  # Add line number tracking
                line = line.strip()  # Remove leading/trailing whitespace
                if not line or line.startswith("#"):
                    continue  # Skip comment lines
                
                fields = line.strip().split('\t')
                if len(fields) < 9:  # Ensure there are at least 9 columns in GFF3 format
                    raise ValueError(f"Invalid GFF3 format in line {line_number}: {line.strip()}")
                
                seqid, source, feature_type, start, end, score, strand, phase, attributes_str = fields
                
                # Convert start and end to integers
                try:
                    start = int(start)
                    end = int(end)
                except ValueError:
                    raise ValueError(f"Start and end positions must be integers in line {line_number}: {line.strip()}")
                
                # Check seqid exists in FASTA
                if seqid not in fasta_seqs:
                    raise ValueError(f"Error: GFF3 seqid '{seqid}' not found in FASTA file at line {line_number}.")
                
                # Check start is less than end
                if start > end:
                    raise ValueError(
                        f"Error: Start position ({start}) must be less than end position ({end}) "
                        f"in line {line_number}: {line.strip()}"
                    )
                
                # Check strand is valid
                if strand not in ["+", "-"]:
                    raise ValueError(
                        f"Error: Strand must be '+' or '-' but found '{strand}' in line {line_number}: {line.strip()}"
                    )
                
                # Check annotation range within FASTA sequence length
                if start < 1 or end > fasta_seqs[seqid]:
                    raise ValueError(
                        f"Error: GFF3 annotation for seqid '{seqid}' has an invalid range "
                        f"({start}-{end}), exceeds FASTA sequence length ({fasta_seqs[seqid]}) in line {line_number}: {line.strip()}"
                    )
                
                # Parse attributes
                attr_dict = parse_attributes(attributes_str)
                
                # Check feature-specific attributes
                if feature_type == "gene":
                    has_gene_feature = True
                    if "ID" not in attr_dict:
                        raise ValueError(f"Error: 'gene' feature must have an 'ID' attribute in line {line_number}: {line.strip()}")
                
                elif feature_type == "mRNA":
                    if "ID" not in attr_dict or "Parent" not in attr_dict:
                        raise ValueError(
                            f"Error: 'mRNA' feature must have both 'ID' and 'Parent' attributes in line {line_number}: {line.strip()}"
                        )
                
                else:  # Other features
                    if "Parent" not in attr_dict:
                        raise ValueError(
                            f"Error: Feature type '{feature_type}' must have a 'Parent' attribute in line {line_number}: {line.strip()}"
                        )
        
        if not has_gene_feature:
            raise ValueError("Error: No 'gene' feature found in the GFF3 file.")
        
        print("Validation successful: GFF3 annotations match the FASTA sequences.")
    
    except Exception as e:
        print(f"Validation failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python validate_genome.py <gff3_file> <fasta_file>")
        sys.exit(1)
    
    gff3_file = sys.argv[1]
    fasta_file = sys.argv[2]
    check_gff3_with_fasta(gff3_file, fasta_file)
