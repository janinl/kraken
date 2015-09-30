#ifndef CONFIG_PARAMS_HPP
#define CONFIG_PARAMS_HPP

#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<map>

using namespace std;

namespace configuration {
	
	struct config_parameters {
	
		string kraken_output_dir;
		string taxid_to_genus_filename;
		int kmer_length;
		string sequences_type;
	
	};
	
	void readConfig(config_parameters& config) {
		ifstream in("../config/config.txt");
		string line;
		cout << "here" << endl;
		while (getline(in, line)) {
			cout << line << endl;
			istringstream s(line.substr(line.find("=")+1));
			if (line.find("kraken_output_dir") != -1){
				s >> config.kraken_output_dir;
			}
			else if (line.find("taxid_to_genus_filename") != -1) {
				s >> config.taxid_to_genus_filename;
			}
			else if (line.find("kmer_length") != -1) {
				s >> config.kmer_length;
			}
			else if (line.find("sequences_type") != -1) {
				s >> config.sequences_type;
			}
		}
		in.close();
	}
	
}
#endif