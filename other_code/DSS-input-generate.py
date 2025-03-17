import gzip
import csv
import argparse

def process_cx_report(input_file, output_file):
    """
    将每个CpG位点的两条链数据合并成一行
    
    参数:
    input_file: 输入的压缩CX报告文件路径
    output_file: 输出文件路径
    """
    with gzip.open(input_file, 'rt') as f_in, open(output_file, 'w', newline='') as f_out:
        # 创建CSV读写器
        reader = csv.reader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        
        # 写入表头
        writer.writerow(['chr', 'pos', 'N', 'X', 'proportion'])
        
        # 逐对处理行
        odd_line = None
        for line in reader:
            if odd_line is None:
                # 存储奇数行
                odd_line = line
            else:
                # 处理偶数行并与前一个奇数行合并
                even_line = line
                
                # 确保这是同一个CpG位点的两条链
                if odd_line[0] == even_line[0] and int(odd_line[1]) == int(even_line[1]) - 1:
                    chr_name = odd_line[0]
                    # zero-based
                    pos_odd = int(odd_line[1]) - 1
                    # one-based
                    pos = pos_odd + 1
                    
                    # 合并甲基化和非甲基化计数
                    methylated_count = int(odd_line[3]) + int(even_line[3])
                    unmethylated_count = int(odd_line[4]) + int(even_line[4])
                    
                    # 计算甲基化比例
                    total_count = methylated_count + unmethylated_count
                    proportion = methylated_count / total_count if total_count >= 3 else -1
                    
                    # 写入合并后的行
                    if proportion >= 0:
                        writer.writerow([chr_name, pos, total_count, methylated_count, proportion])
                
                # 重置奇数行
                odd_line = None

def main():
    # 设置输入和输出文件路径
    
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='处理CX报告文件并生成DSS输入')
    parser.add_argument('input_file', help='输入的压缩CX报告文件路径')
    parser.add_argument('output_file', help='输出文件路径')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    
    print(f"开始处理文件: {input_file}")
    process_cx_report(input_file, output_file)
    print(f"处理完成，结果已保存至: {output_file}")

if __name__ == "__main__":
    main()