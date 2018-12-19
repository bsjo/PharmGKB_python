#!/usr/bin/env python

import requests, StringIO, zipfile, os, datetime, re, string



class PharmGKB_parsing():
	def __init__(self):
		print("START\n")
		temp_today = datetime.datetime.today()
		self.today = str(temp_today).split(" ")[0]
		print(temp_today)

	def PharmGKB_download(self, url_CLI, url_ANNO, input_path):
		data = requests.get(url_CLI)
		f = StringIO.StringIO()
		f.write(data.content)
		input_zip = zipfile.ZipFile(f)
		r_data = [input_zip.read(i) for i in input_zip.namelist() if i == 'clinical_ann_metadata.tsv'][0]
		self.r_list = r_data.splitlines()
		#PharmGKB making
		if not os.path.exists(input_path): os.makedirs(input_path)
		PharmGKB_down = '{0}/PharmGKB_ClinicalVariants_{1}'.format(input_path, self.today)
		with open(PharmGKB_down, 'w') as w:
			w.write(r_data)

		#PharmGKB variant annotation
		data_a = requests.get(url_ANNO)
		f_a = StringIO.StringIO()
		f_a.write(data_a.content)
		input_zip_a = zipfile.ZipFile(f_a)
		r_data_a = [input_zip_a.read(i) for i in input_zip_a.namelist() if i == 'variants.tsv'][0]
		self.variant_anno = r_data_a.splitlines()
		PharmGKB_down = '{0}/PharmGKB_Variants_{1}'.format(input_path, self.today)
		with open(PharmGKB_down, 'w') as w:
			w.write(r_data_a)

	def PharmGKB_processing(self, output_path):
		if not os.path.exists(output_path): os.makedirs(output_path)
		output_path_n = os.path.join(output_path, "Parsed_PharmGKB_{0}".format(self.today))

		anno_hash = {}
		for line in self.variant_anno:
			line_split = line.split('\t')
			Variant_key = line_split[1] ; Location = line_split[5]
			if not re.search('GRCh38|hg38', Variant_key):
				if Location:
					Location = re.sub("\[GRCh37\]chr\d+\:", "", Location)
					anno_hash[Variant_key] = Location

		with open(output_path_n, "w") as w:
			header = '\t'.join(['Chromosome', 'Position', 'Variant', 'Gene', 'Disease', 'Drug', 'Drug_Type', 'Evidence', 'Evidence_Count']) + '\n'
			w.write(header)
			for line in self.r_list[1:]:
				line_split = line.split('\t')
				Chrom=line_split[14]; Variant=line_split[1]; Gene=line_split[2];
				Evidence=line_split[3]; Evidence_Count=line_split[10]; Drug_Type=line_split[4];
				Drug=line_split[11]; Disease=line_split[12]
				Anno_location = ""
				if Variant in anno_hash:
					Anno_location = anno_hash[Variant]
				out_result = '\t'.join([Chrom, Anno_location, Variant, Gene, Disease, Drug, Drug_Type, Evidence, Evidence_Count]) + '\n'
				w.write(out_result)

if __name__ == '__main__':

	url_CLI = "https://s3.pgkb.org/data/annotations.zip"
	url_ANNO = "https://s3.pgkb.org/data/variants.zip"
	input_path = 'Input'
	output_path = 'Output'

	#PharmGKB_Class
	go = PharmGKB_parsing()
	#PharmGKB_CancerBiomarkers_Download_lastest
	go.PharmGKB_download(url_CLI, url_ANNO, input_path)
	#PharmGKB_Parsing
	go.PharmGKB_processing(output_path)

	print("\nEND")