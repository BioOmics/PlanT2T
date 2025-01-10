import argparse
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

def main(input_file, output_file):
    with open(output_file, 'w') as outfile:
        outfile.write("SequenceName\tMolecularWeight\tAromaticity\tInstabilityIndex\tIsoelectricPoint\tHelixFraction\tTurnFraction\tSheetFraction\tReducedCysteinesExtinctionCoefficient\tOxidizedCysteinesExtinctionCoefficient\tGRAVY\tAverageFlexibility\tChargeAtpH7.0\n")
        
        for seq in SeqIO.parse(input_file, 'fasta'):
            protein_sequence = seq.seq
            # Check if the sequence is empty to avoid ZeroDivisionError
            if len(protein_sequence) == 0:
                print(f"Warning: Empty sequence for {seq.id}, skipping.")
                continue  # Skip empty sequences
                
            prot_param = ProteinAnalysis(protein_sequence)

            molecular_weight = prot_param.molecular_weight()
            aromaticity = prot_param.aromaticity()
            instability_index = prot_param.instability_index()
            isoelectric_point = prot_param.isoelectric_point()
            sec_struc = prot_param.secondary_structure_fraction()  # [helix, turn, sheet]
            epsilon_prot = prot_param.molar_extinction_coefficient()  # [reduced, oxidized]
            gravy = prot_param.gravy()
            flexibility = prot_param.flexibility()
            if flexibility:
                average_flexibility = sum(flexibility) / len(flexibility)
            else:
                average_flexibility = 0.0  # Default value if flexibility cannot be computed
            charge_at_pH_7 = prot_param.charge_at_pH(7.0)

            outfile.write(f"{seq.id}\t")
            outfile.write(f"{molecular_weight:0.2f}\t")
            outfile.write(f"{aromaticity:0.2f}\t")
            outfile.write(f"{instability_index:0.2f}\t")
            outfile.write(f"{isoelectric_point:0.2f}\t")
            outfile.write(f"{sec_struc[0]:0.2f}\t")  # Helix fraction
            outfile.write(f"{sec_struc[1]:0.2f}\t")  # Turn fraction
            outfile.write(f"{sec_struc[2]:0.2f}\t")  # Sheet fraction
            outfile.write(f"{epsilon_prot[0]}\t")  # Reduced cysteines extinction coefficient
            outfile.write(f"{epsilon_prot[1]}\t")  # Oxidized cysteines extinction coefficient
            outfile.write(f"{gravy:0.2f}\t")  # GRAVY
            outfile.write(f"{average_flexibility:0.2f}\t")  # Average Flexibility
            outfile.write(f"{charge_at_pH_7:0.2f}\n")  # Charge at pH 7.0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze protein parameters from a FASTA file.")
    parser.add_argument('input_file', type=str, help="Path to the input protein FASTA file.")
    parser.add_argument('output_file', type=str, help="Path to the output file where results will be written.")

    args = parser.parse_args()

    main(args.input_file, args.output_file)
