import argparse
from Bio import SeqIO

def calculate_gc_content(sequence):
    """
    计算序列的GC含量百分比，排除'N'碱基
    :param sequence: 输入的DNA序列
    :return: GC含量的百分比
    """
    # 计算G和C的总数
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    
    # 计算A、T、C、G的总数，排除'N'
    total_bases = sequence.count('A') + sequence.count('T') + g_count + c_count
    
    if total_bases == 0:
        return 0.0  # 防止除零错误
    
    # 计算GC含量的百分比
    gc_content = (g_count + c_count) / total_bases * 100
    return gc_content

def calculate_genome_gc(fasta_file):
    """
    计算全基因组的GC含量
    :param fasta_file: 输入的FASTA文件路径
    :return: 全基因组GC含量百分比
    """
    total_sequence = ""
    
    # 合并所有染色体的序列
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_sequence += str(record.seq)
    
    # 计算全基因组GC含量
    return calculate_gc_content(total_sequence)

def calculate_chromosome_gc(fasta_file):
    """
    计算FASTA文件中每个染色体的GC含量，排除'N'碱基
    :param fasta_file: 输入的FASTA文件路径
    :return: 字典，键为染色体名，值为GC含量的百分比
    """
    gc_content_per_chromosome = {}
    
    # 读取FASTA文件中的每条序列
    for record in SeqIO.parse(fasta_file, "fasta"):
        chromosome = record.id
        sequence = str(record.seq)
        
        # 计算该染色体序列的GC含量
        gc_content = calculate_gc_content(sequence)
        gc_content_per_chromosome[chromosome] = gc_content
    
    return gc_content_per_chromosome

def main():
    # 使用argparse解析命令行参数
    parser = argparse.ArgumentParser(description="Calculate GC content for each chromosome in a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the FASTA file containing the genome sequences")
    args = parser.parse_args()
    
    # 获取FASTA文件路径
    fasta_file = args.fasta_file
    
    # 计算全基因组GC含量
    genome_gc = calculate_genome_gc(fasta_file)
    print(f"Genome\t{genome_gc:.2f}%")
    
    # 计算每个染色体的GC含量
    gc_content = calculate_chromosome_gc(fasta_file)
    
    # 输出每个染色体的GC含量
    for chromosome, gc in gc_content.items():
        print(f"{chromosome}\t{gc:.2f}%")

if __name__ == "__main__":
    main()
