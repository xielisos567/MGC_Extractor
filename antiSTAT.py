import argparse
import os
import glob
from bs4 import BeautifulSoup
import pandas as pd

def main(input_dir, output_file):
	#dir_path = os.listdir(input_dir)
	# print(dir_path)
	# dir_path = ["SRR10156263_bin_1", "SRR10156267_bin_2"]
	all_data = pd.DataFrame(None, columns=["MAGs", "Genome", "Region", "Type", "From", "To"])
	print (all_data)
	for infile in glob.glob(f'{input_dir}/*/index.html'):
		print (infile)
		with open(infile, "r", encoding="utf-8") as file:
			html_content = file.read()
			soup = BeautifulSoup(html_content, 'html.parser')
			table = soup.find("div", id="single-record-tables")
			# 重置列表
			mags_list = []
			genome_list = []
			region_list = []
			type_list = []
			from_list = []
			to_list = []
			try:
				tables = table.find_all("tr", class_="linked-row odd")
			except:
				continue
			# 为每行数据提取相关信息
			for row in table.find_all("tr", class_="linked-row odd"):
				# MAGs 值
				mags = os.path.abspath(infile).split('/')[-2]
				print (mags)
				# 从 "record-overview-header" 提取 Genome 值
				header_div = row.find_previous("div", class_="record-overview-header")
				genome = header_div.text.split("(original name was: ")[-1].replace(")", "").strip() if header_div else None

				# 提取 Region, Type, From, 和 To 值
				region = row.find("td", class_="regbutton").a.text.replace("&nbsp", "").strip()
				type_ = row.find("a", class_="external-link").text
				from_ = int(row.find("td", class_="digits").text.replace(",",""))
				to = int(row.find("td", class_="digits table-split-left").text.replace(",", ""))

				# 将提取的值添加到相应的列表中
				mags_list.append(mags)
				genome_list.append(genome)
				region_list.append(region)
				type_list.append(type_)
				from_list.append(from_)
				to_list.append(to)

			# 使用pandas将数据整合到一个DataFrame中
			final_df = pd.DataFrame({
				"MAGs": mags_list,
				"Genome": genome_list,
				"Region": region_list,
				"Type": type_list,
				"From": from_list,
				"To": to_list
			})
			if len(all_data) == 0:
				all_data = final_df
			else:
				all_data = pd.concat([all_data, final_df], ignore_index=True)
			print(f"============加载数据:{len(final_df)}条============")
	all_data.to_excel(output_file, index=False)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='数据采集')
	parser.add_argument('-i', '--input', required=True, help='输入文件夹路径')
	parser.add_argument('-o', '--output', required=True, help='输出文件路径')
	args = parser.parse_args()

	main(args.input, args.output)
