#include<dirent.h>
#include<iostream>
#include<stdio.h>
#include<cstdlib>
#include<iostream>
#include<string.h>
#include<fstream>
#include<string>
#include<sstream>
#include<unordered_map>
#include<inttypes.h>
#include<string>
#include<fstream>
#include<streambuf>
#include<vector>
#include<iterator>
#include<cerrno>
#include<unordered_set>
#include<math.h>
#include<algorithm>
#include<omp.h>
#include <sys/stat.h>
#include "config_constants.hpp"

using namespace std;
using namespace kraken_constants;	

struct genus_divergence {
	string genus;
	double divergence;
};

genus_divergence get_KL_divergence_for_genus(const char * genus_dir_path) {

	string genus_dir_path_str(genus_dir_path);
	string genus_dir_path_substr = genus_dir_path_str.substr(0, genus_dir_path_str.size()-1);
	size_t delim = genus_dir_path_substr.find_last_of("/");
	string genus = genus_dir_path_substr.substr(delim+1);
	cout << "Processing: " << genus << endl;
	
	DIR *dir;
	struct dirent *ent;
	int numThreadsReadingSameFile = 10;
	omp_set_num_threads(numThreadsReadingSameFile);	
	vector<unordered_map<uint64_t, long long int>> kmer_counts_vect( numThreadsReadingSameFile, unordered_map<uint64_t, long long int>() );
	
	if ((dir = opendir (genus_dir_path)) != NULL) {
		while ((ent = readdir (dir)) != NULL) {
			string fileName = ent->d_name;
			if (fileName != "." && fileName != ".." ){
				
				string filePath = genus_dir_path + fileName;
				ifstream in(filePath);
				string contents((istreambuf_iterator<char>(in)),istreambuf_iterator<char>());
				in.close();
				
				//several threads read the same file content:
				#pragma omp parallel for 
				for (int thNum = 0; thNum < numThreadsReadingSameFile; ++thNum) {
					stringstream datastream(contents);
					uint64_t kmer;
					while(datastream >> kmer) {
						//each thread works on the kmers ending with a particular digit
						int key = kmer % numThreadsReadingSameFile;
						if (key == thNum) {
							kmer_counts_vect[thNum][kmer]++;

						}
					}			
				}
			}
		}
	}
	
	else {
		perror ("");
	}
	closedir(dir);
		
	long long int big_sum = 0;
	long long int big_kmer_count = 0;
	for (int i = 0; i < numThreadsReadingSameFile; i++) {
		unordered_map<uint64_t, long long int> kmer_counts = kmer_counts_vect[i];
		long long int sum = accumulate(begin(kmer_counts), end(kmer_counts), 0, [](long long int prev, const pair<uint64_t,long long int>& p) { return prev+p.second; });
		big_sum += sum;
		big_kmer_count += kmer_counts.size();
	}
		
	double expectedRelFrequency = 1.0/big_kmer_count;
	double KLdivergence = 0.0;
		
	for (int i = 0; i < numThreadsReadingSameFile; i++) {
		unordered_map<uint64_t, long long int> &kmer_counts = kmer_counts_vect[i];
		for (auto &kmerCount : kmer_counts) {
			double kmerRelFrequency = (double)kmerCount.second/big_sum;
			KLdivergence += kmerRelFrequency * log(kmerRelFrequency/expectedRelFrequency);
		}
	}	
		
	genus_divergence genus_div;
	genus_div.genus = genus;
	genus_div.divergence = KLdivergence;
	
	return genus_div;
}


void get_divergence_file(string kmer_dir_path_str) {
	
	DIR *dir = opendir(kmer_dir_path_str.c_str());
	struct dirent *entry = readdir(dir);
	vector<string> genus_dir_paths;
    while (entry != NULL){
        if ((entry->d_type == DT_DIR)){
			printf("%s\n", entry->d_name);
			string dir_name(entry->d_name);
			if ((dir_name != ".") && (dir_name != "..")) {
				string genus_dir_name(entry->d_name);
				string genus_dir_path_str = kmer_dir_path_str +  genus_dir_name + "/";
				genus_dir_paths.push_back(genus_dir_path_str);
			}
		}	
        entry = readdir(dir);
    }
    for (int i = 0; i < genus_dir_paths.size(); i++) {
		cout << genus_dir_paths[i] << endl;
    }
    
    //find path to the results directory:
	string kmer_dir_path_substr = kmer_dir_path_str.substr(0, kmer_dir_path_str.size()-1);
	size_t delim = kmer_dir_path_substr.find_last_of("/");
	string fasta_dir_path = kmer_dir_path_substr.substr(0,delim);
	//cout << "fasta dir : " << fasta_dir_path << endl;
    
	string results_dir_path_str = fasta_dir_path + "/ResultsFolder"; 
	struct stat sb;
	
	//check if the ResultsFolder exists and if so write output to divergence file:
	if (stat(results_dir_path_str.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
		ofstream divergence_file(results_dir_path_str + "/KLdivergence.div");
		divergence_file << "Genus,KLdivergence\n";	
		for (int i = 0; i < genus_dir_paths.size(); i++) {
			genus_divergence genus_div = get_KL_divergence_for_genus(genus_dir_paths[i].c_str());
			divergence_file << genus_div.genus << "," << genus_div.divergence<<endl;
		}
		divergence_file.close();
		closedir(dir);
	}
	
	else {
		cout << results_dir_path_str << " does not exist!" << endl;
	}
}

	
int main() {
	
	DIR *dir = opendir(kraken_output_dir.c_str());
	vector<string> kmer_dir_paths;
	struct dirent *entry = readdir(dir);
	struct stat sb;
	while (entry != NULL) {
		if (entry -> d_type == DT_DIR) {
			string fasta_dir_str(entry -> d_name);
			//cout << fasta_dir_str << endl;
			if ((fasta_dir_str != ".") && (fasta_dir_str != "..")) {
				string kmer_dir_path_str = kraken_output_dir + fasta_dir_str + "/KmersFolder/";
				//cout << kmer_dir_path_str << endl;
				if (stat(kmer_dir_path_str.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
					kmer_dir_paths.push_back(kmer_dir_path_str);
				}
			}
		}
		entry = readdir(dir);
	}
	closedir(dir);
	
	for (int i = 0; i < kmer_dir_paths.size(); i++) {
		//cout << kmer_dir_paths[i] << endl;
		get_divergence_file(kmer_dir_paths[i]);
	}
	
	return 0;
}