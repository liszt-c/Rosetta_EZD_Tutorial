import csv

def parse_sc_file(sc_file_path, csv_file_path):
    # 打开输入和输出文件
    with open(sc_file_path, 'r') as infile, open(csv_file_path, 'w', newline='') as outfile:
        # 初始化CSV写入器
        writer = csv.writer(outfile)
        
        headers_written = False
        
        for line in infile:
            if line.startswith("SEQUENCE:") or line.startswith("SCORE:"):
                # 去除开头的'SCORE:'或'SEQUENCE:'
                clean_line = line.replace("SCORE:", "").replace("SEQUENCE:", "").strip()
                # 分割行成字段列表
                fields = clean_line.split()

                # 如果是第一行，则写入表头
                if not headers_written and "total_score" in clean_line:
                    header = [field.strip() for field in fields]
                    writer.writerow(header)
                    headers_written = True
                else:
                    # 写入数据行，尝试将数值转换为浮点数
                    writer.writerow([float(x) if x.replace('.','',1).isdigit() or (x.startswith('-') and x[1:].replace('.','',1).isdigit()) else x for x in fields])


# 使用方法
sc_file_path = "./design_results.sc"
csv_file_path = "./design_results.csv"

parse_sc_file(sc_file_path, csv_file_path)