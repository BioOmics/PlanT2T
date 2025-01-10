import argparse
import os
# output: [species_name] [total] [one_dot_end] [one_dot_not_end_or_multiple] [no_dot] 

def count_proteins(file_path):
    total_count = 0
    one_dot_end_count = 0
    one_dot_not_end_or_multiple_count = 0
    no_dot_count = 0
    
    with open(file_path, 'r') as file:
        in_sequence = False
        current_sequence = ""
        
        for line in file:
            if line.startswith('>'):
                if in_sequence:
                    total_count += 1
                    if current_sequence.count('.') == 1:
                        if current_sequence.endswith('.'):
                            one_dot_end_count += 1
                        else:
                            one_dot_not_end_or_multiple_count += 1
                    elif current_sequence.count('.') > 1:
                        one_dot_not_end_or_multiple_count += 1
                    else:
                        no_dot_count += 1
                current_sequence = ""
                in_sequence = True
            elif in_sequence:
                current_sequence += line.strip()
        
        # Check last sequence
        if in_sequence:
            total_count += 1
            if current_sequence.count('.') == 1:
                if current_sequence.endswith('.'):
                    one_dot_end_count += 1
                else:
                    one_dot_not_end_or_multiple_count += 1
            elif current_sequence.count('.') > 1:
                one_dot_not_end_or_multiple_count += 1
            else:
                no_dot_count += 1

    return total_count, one_dot_end_count, one_dot_not_end_or_multiple_count, no_dot_count

def main():
    parser = argparse.ArgumentParser(description='Count protein sequences with and without "."')
    parser.add_argument('file', type=str, help='Path to the protein sequence file')
    args = parser.parse_args()
    
    # Extract species name from the file name
    species_name = os.path.splitext(os.path.basename(args.file))[0]
    
    total, one_dot_end, one_dot_not_end_or_multiple, no_dot = count_proteins(args.file)
    print(f'{species_name}\t{total}\t{one_dot_end}\t{one_dot_not_end_or_multiple}\t{no_dot}')

if __name__ == '__main__':
    main()

