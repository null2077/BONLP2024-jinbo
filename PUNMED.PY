#获取pmid：
import requests
import xml.etree.ElementTree as ET
# 发送请求
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
params = {
    'db': 'pubmed',
    'term': 'cotton fiber [MeSH Terms] and gene',
    'mindate': '2000/01/01',
    'maxdate': '2024/10/15',
    'retmode': 'xml',
    'retmax': '100000'# 返回 XML 格式
}
response = requests.get(url, params=params)
response.raise_for_status()  # 检查请求是否成功
# 解析 XML 数据
root = ET.fromstring(response.content)
# 提取 Id
ids = [id_elem.text for id_elem in root.findall('.//Id')]
# 将 Id 写入文件
with open('curl_id.txt', 'w') as f:
    for id in ids:
        f.write(f"{id}\n")
# 打印 Id 列表
print(ids)

文献信息提取
import requests
import time
# 定义
F_OUT = "result_Pubtator.txt"
F_LIST = "curl_id.txt"
API_KEY = "c4cddb499c834604eabdec17750233b94608"
# 打印
with open(F_OUT, 'w') as f_out:
    f_out.write("In\n") 
# 读取
with open(F_LIST, 'r') as f_list:
    pmids = f_list.readlines()
# 循环
for i, line in enumerate(pmids, start=1):
   pmid = line.strip() 
   url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={pmid}"
     #头部
headers = {'Authorization': f'Bearer {API_KEY}'     }
    try:
    response = requests.get(url)
         response.raise_for_status()  
# 检查
        with open(F_OUT, 'a') as f_out: 
# 追加
            f_out.write(response.text)  
# 写入
        print(f"{i} -th result out of {len(pmids)} is processing..")
    except requests.exceptions.RequestException as e:
        print(f"Error processing PMID {pmid}: {e}")
    time.sleep(1)

统计gene类型个数及云图实现
#统计gene类型个数：
from collections import defaultdict
# 初始化一个字典来存储基因类型及其出现次数
gene_counts = defaultdict(int)
# 读取文件
with open('result_Pubtator.txt', 'r', encoding='utf-8') as file:
    # 读取所有行
    lines = file.readlines()
# 遍历
for line in lines:
    # 假设基因类型和"Gene"之间可能有空格，并假设基因类型后面直接跟着"Gene"
    # 使用正则表达式来匹配基因类型
    import re
    match = re.search(r'([A-Z]+\w+)\s*Gene', line)
    if match:
        # 获取匹配到的基因类型
        gene_type = match.group(1)
        # 更新基因类型计数
        gene_counts[gene_type] += 1
# 打印基因类型及其出现次数
for gene_type, count in gene_counts.items():
    print(f"Gene Type: {gene_type}, Count: {count}")
# 如果需要，可以打印到文件中
with open('gene_counts.txt', 'w', encoding='utf-8') as output_file:
    for gene_type, count in gene_counts.items():
        output_file.write(f"Gene Type: {gene_type}, Count: {count}\n")

生成词云图：
import re
from wordcloud import WordCloud
import matplotlib.pyplot as plt
# 读取基因类型和计数数据
with open('gene_counts.txt', 'r', encoding='utf-8') as file:
    gene_data= file.read()
# 使用正则表达式提取基因名称和它们的出现次数
gene_counts = {}
for line in gene_data.strip().split('\n'):
    match = re.search(r'Gene Type: (\w+), Count: (\d+)', line)
    if match:
        gene_name, count = match.groups()
        gene_counts[gene_name] = int(count)
# 创建一个以词频为字典的字符串，用于生成词云
text = ' '.join([gene + ' ' * count for gene, count in gene_counts.items()])
# 创建词云对象
wordcloud = WordCloud(width=800, height=400, background_color='white').generate(text)
# 显示词云图像
plt.figure(figsize=(10, 5))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis('off')
plt.show()

import sys
import time
import requests
from urllib.parse import quote
import re
import os
from bs4 import BeautifulSoup

# 读取文件并提取计数大于等于4的基因
def extract_genes(file_path):
    genes = []
    with open(file_path, 'r') as file:
        for line in file:
            if "Count:" in line:
                parts = line.split(', ')
                gene_type = parts[0].split(': ')[1].strip()
                count = int(parts[1].split(': ')[1].strip())
                if count >= 4:
                    genes.append(gene_type)
    return genes

def get_pathway_info(gene_list):
    kegg_ids = []
    for gene in gene_list:
        try:
            # 对基因名称进行url编码
            encoded_gene = quote(gene)
            query_url = f"http://rest.kegg.jp/find/genes/{encoded_gene}"
            response = requests.get(query_url, timeout=10)
            if response.status_code == 200:
                # 获取网页内容
                content = response.text
    
                # 使用正则表达式查找所有以 "ath:" 开头的KEGG ID
                found_kegg_ids = re.findall(r'ath:AT[0-9G]+', content)
                kegg_ids.extend(found_kegg_ids)
    
                print(f'KEGG IDs for {gene} have been found.')
            else:
                print(f'Failed to retrieve data for {gene}, HTTP status code: {response.status_code}')
        except requests.exceptions.RequestException as e:
            print(f'Request failed for gene: {gene}, error: {e}')
        time.sleep(1)  # 避免过于频繁的请求
    return kegg_ids


# 函数用于提取基因家族名
def extract_gene_family_from_kegg(kegg_id):
    url = f"http://www.kegg.jp/entry/{kegg_id}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            # 寻找包含"Symbol"的<th>标签
            symbol_label = soup.find('th', class_='th11 deft tal vtop', text='Symbol')
            if symbol_label:
                # 获取"Symbol"对应的<td>标签
                gene_family_cell = symbol_label.find_next_sibling('td', class_='td11 defd')
                if gene_family_cell:
                    # 提取基因家族名
                    gene_family_name = gene_family_cell.div.get_text(strip=True)
                    return gene_family_name
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for KEGG ID {kegg_id}: {e}")
    return "Not found"

# 主函数
def main():
    file_path = 'gene_counts.txt'
    gene_list = extract_genes(file_path)
    kegg_ids = get_pathway_info(gene_list)

    # KEGG ID文件的路径
    kegg_ids_file_path = 'kegg_ids.txt'

    # 将KEGG ID写入文件
    with open(kegg_ids_file_path, 'w') as file:
        for kegg_id in kegg_ids:
            file.write(kegg_id + '\n')
    print(f'All KEGG IDs have been written to {kegg_ids_file_path}')

    # 读取文件中的KEGG ID并提取基因家族名
    with open(kegg_ids_file_path, 'r') as file:
        kegg_ids = file.read().splitlines()

    # 存储KEGG ID和对应的基因家族名
    gene_families = {}

    for kegg_id in kegg_ids:
        gene_family = extract_gene_family_from_kegg(kegg_id)
        gene_families[kegg_id] = gene_family
        print(f"KEGG ID: {kegg_id}, Gene family: {gene_family}")

    # 如果需要，可以将结果写入到文件
    with open('gene_families.txt', 'w') as file:
        for kegg_id, gene_family in gene_families.items():
            file.write(f"{kegg_id}\t{gene_family}\n")

# 运行主函数
if __name__ == "__main__":
main()
