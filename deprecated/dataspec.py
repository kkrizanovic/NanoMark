#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys;
import subprocess;

import basicdefines

if __name__ == "__main__":
	pass;

def execute_command(command):
	sys.stderr.write('[DataSpec] Running command: "%s".\n' % (command));
	subprocess.call(command, shell=True);

### Hack to mimic Enum.
class DataType:
	[none, nanopore, pacbio, illumina, hybridnanopore, hybridpacbio] = range(6);
	all_types_string = ['none', 'nanopore', 'pacbio', 'illumina', 'hybridnanopore', 'hybridpacbio'];
	string_type_to_enum = {'none': none, 'nanopore': nanopore, 'pacbio': pacbio, 'illumina': illumina, 'hybridnanopore': hybridnanopore, 'hybridpacbio': hybridpacbio};

class DataSpec:
	### Parameter description:
	###		- name - the name of the dataset; will be used to name folders for the data and the results.
	###		- data_folder - folder where raw data will be downloaded, unpacked and converted if necessary.
	###		- data_type - specifies whether the dataset should be taken as a nanopore, PacBio, Illumina, or a hybrid combination. Use one of the types defined in DataType class.
	###		- nanopore_files - filenames specifying reads obtained by nanopore sequencing.
	###		- pacbio_files - filenames specifying reads obtained by PacBio sequencing.
	###		- illumina_files - filenames specifying reads obtained by Illumina sequencing.
	###		- data_urls - URLs to raw data files that will be downloaded to path ($data_folder/$name).
	def __init__(self, name, data_folder, data_type, nanopore_files, pacbio_files, illumina_files, data_urls):
		self.name = name;
		self.data_type = data_type;	### nanopore, pacbio, illumina, hybridnanopore, hybridpacbio
		self.data_folder = data_folder;
		self.nanopore_files = nanopore_files;
		self.pacbio_files = pacbio_files;
		self.illumina_files = illumina_files;
		self.data_urls = data_urls;

	def clear(self):
		self.name = '';
		self.data_type = DataType.none;
		self.data_folder = '';
		self.nanopore_files = [];
		self.pacbio_files = [];
		self.illumina_files = [];
		self.data_urls = [];

	def download(self):
		for data_url in self.data_urls:
			execute_command('cd %s; wget %s' % (self.data_folder, self.data_url));

	def save_spec(self, out_file):
		fp = None;
		try:
			fp = open(out_file, 'w');
		except Exception, e:
			sys.stderr.write('ERROR: Could not open file "%s" for writing!' % (out_file));
			return;

		fp.write('name: %s\n' % self.name);
		fp.write('data_type: %s\n' % (DataType.all_types_string[self.data_type]));
		fp.write('data_folder: %s\n' % (self.data_folder));
		fp.write('nanopore_files: %s\n' % (', '.join(self.nanopore_files)));
		fp.write('pacbio_files: %s\n' % (', '.join(self.pacbio_files)));
		fp.write('illumina_files: %s\n' % (', '.join(self.illumina_files)));
		fp.write('data_urls: %s\n' % (', '.join(self.data_urls)));

		fp.close();

	def load_spec(self, in_file):
		fp = None;
		try:
			fp = open(in_file, 'r');
		except Exception, e:
			sys.stderr.write('ERROR: Could not open file "%s" for reading!' % (in_file));
			return;

		for line in fp:
			line = line.strip();
			if (len(line) == 0):
				continue;
			split_line = line.split(': ');
			if (split_line[0] == 'name'):
				self.name = ': '.join(split_line[1:]);	### The ': '.join() is here because such pattern could have occured in the string as well.
			elif (split_line[0] == 'data_type'):
				self.data_type = DataType.string_type_to_enum[': '.join(split_line[1:])];
			elif (split_line[0] == 'data_folder'):
				self.data_folder = ': '.join(split_line[1:]);
			elif (split_line[0] == 'nanopore_files'):
				self.nanopore_files = (': '.join(split_line[1:])).split(', ');
			elif (split_line[0] == 'pacbio_files'):
				self.pacbio_files = (': '.join(split_line[1:])).split(', ');
			elif (split_line[0] == 'illumina_files'):
				self.illumina_files = (': '.join(split_line[1:])).split(', ');
			elif (split_line[0] == 'data_urls'):
				self.data_urls = (': '.join(split_line[1:])).split(', ');

		fp.close();



dataset1 = DataSpec(name='dataset1_30x2d',
					data_folder='%s/' % (basicdefines.DATASETS_PATH),
					data_type=DataType.nanopore,
					nanopore_files=['ERX708228.fastq', 'ERX708229.fastq', 'ERX708230.fastq', 'ERX708231.fastq'],
					pacbio_files=[],
					illumina_files=[],
					data_urls=[	'ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar',
								'ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar',
								'ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar',
								'ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar']);

dataset1.save_spec('../temp/test1.spec');
# dataset1.clear();
# dataset1.save_spec('../temp/test2.spec');
# dataset1.load_spec('../temp/test1.spec');
# dataset1.save_spec('../temp/test3.spec');
