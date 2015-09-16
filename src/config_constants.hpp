#ifndef MYLIB_CONSTANTS_HPP
#define MYLIB_CONSTANTS_HPP

#include<string>
using namespace std;

namespace kraken_constants {
	
	string kraken_output_dir = "/illumina/scratch/tmp/users/avoicu/krakenCLL/";
	string taxid_to_genus_filename = "/illumina/scratch/tmp/users/avoicu/taxIdToGenus.txt";
	int kmer_length = 31;
}

#endif
