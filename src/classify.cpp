/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"
#include<iterator>
#include<unordered_map>
#include<numeric>
#include<iostream>
#include<string>
#include<utility>
#include<math.h>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<stdio.h>
#include<string.h>
#include<algorithm>
#include <stdio.h>      
#include <stdlib.h>

//#include<EXTERN.h> 
//#include<perl.h> 


using namespace std;
using namespace kraken;

/*
EXTERN_C void boot_DynaLoader (pTHX_ CV* cv);

EXTERN_C void xs_init(pTHX) {
	const char *file = __FILE__;
	newXS("DynaLoader::boot_DynaLoader", boot_DynaLoader, file);
}

*/
//const size_t DEF_WORK_UNIT_SIZE = 500000;
const size_t DEF_WORK_UNIT_SIZE = 100;
/*
static PerlInterpreter *my_perl;
STRLEN strLength;


void getNCBIdatabase() {
	
	const char *importNCBI = "use Bio::LITE::Taxonomy::NCBI;";
	eval_pv(importNCBI, TRUE);
	
	const char *taxNCBI =
		"$taxNCBI = Bio::LITE::Taxonomy::NCBI->new (" 
			"db=>'/home/avoicu/localperl/lib/site_perl/5.22.0/Bio/LITE/Taxonomy/NCBI.pm'," 
			"names=> '/home/avoicu/Bio-LITE-Taxonomy-NCBI-0.07/t/data/names.dmp'," 
			"nodes=>'/home/avoicu/Bio-LITE-Taxonomy-NCBI-0.07/t/data/nodes.dmp'" 
			");";
	eval_pv(taxNCBI, TRUE);
}	

void initializePerlInterpreter(int argc, char **argv, char **env) {
	
	char *embedding[] = { "", "-e", "0" }; 
	PERL_SYS_INIT3(&argc, &argv, &env); 
	my_perl = perl_alloc();
	perl_construct(my_perl); 
	perl_parse(my_perl, xs_init, 3, embedding, NULL); 
	PL_exit_flags |= PERL_EXIT_DESTRUCT_END;
	perl_run(my_perl);

	ENTER;
	SAVETMPS;
	
}

void exitPerlInterpreter() {
	
	FREETMPS;	
	LEAVE;
	perl_destruct(my_perl);
	perl_free(my_perl);
	PERL_SYS_TERM(); 
	
}

string getGenusName(uint64_t predictedTaxID) {
	char getGenusForPredictedTaxID[150];
	strcpy(getGenusForPredictedTaxID, "$predictedGenus = $taxNCBI->get_term_at_level(");
	string taxID;
	ostringstream o;
	o << predictedTaxID;
	taxID += o.str();
	strcat(getGenusForPredictedTaxID, taxID.c_str());
	strcat(getGenusForPredictedTaxID, ",\"genus\");");
	//cout << getGenusForPredictedTaxID << endl;
	eval_pv(getGenusForPredictedTaxID, TRUE);
	const char *genusName = SvPV(get_sv("predictedGenus", FALSE), strLength);
	//eval_pv("$predictedGenus = undef;", TRUE);
	return genusName;
}

*/

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename);
unordered_map<string, unordered_map<uint64_t, int>> classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss);
string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);
set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);

int Num_threads = 1;
string DB_filename, Index_filename, Nodes_filename;
bool Quick_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
uint32_t Minimum_hit_count = 1;
map<uint32_t, uint32_t> Parent_map;
KrakenDB Database;
string Classified_output_file, Unclassified_output_file, Kraken_output_file;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Kraken_output;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;

uint64_t total_classified = 0;
uint64_t total_sequences = 0;
uint64_t total_bases = 0;

unordered_map<uint32_t, string> taxIDtoGenus;


int main(int argc, char **argv, char **env) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif
  
//initializePerlInterpreter(argc, argv, env);
//getNCBIdatabase();

cout << "classify main" << endl;  
parse_command_line(argc, argv);
  if (! Nodes_filename.empty())
    Parent_map = build_parent_map(Nodes_filename);

  if (Populate_memory)
    cerr << "Loading database... ";

  QuickFile db_file;
  db_file.open_file(DB_filename);
  if (Populate_memory)
    db_file.load_file();
  Database = KrakenDB(db_file.ptr());
  KmerScanner::set_k(Database.get_k());

  QuickFile idx_file;
  idx_file.open_file(Index_filename);
  if (Populate_memory)
    idx_file.load_file();
  KrakenDBIndex db_index(idx_file.ptr());
  Database.set_index(&db_index);

  if (Populate_memory)
    cerr << "complete." << endl;

  if (Print_classified) {
    if (Classified_output_file == "-")
      Classified_output = &cout;
    else
      Classified_output = new ofstream(Classified_output_file.c_str());
  }

  if (Print_unclassified) {
    if (Unclassified_output_file == "-")
      Unclassified_output = &cout;
    else
      Unclassified_output = new ofstream(Unclassified_output_file.c_str());
  }

  if (! Kraken_output_file.empty()) {
    if (Kraken_output_file == "-")
      Print_kraken = false;
    else
      Kraken_output = new ofstream(Kraken_output_file.c_str());
  }
  else
    Kraken_output = &cout;

  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (int i = optind; i < argc; i++)
    process_file(argv[i]);
  gettimeofday(&tv2, NULL);

  report_stats(tv1, tv2);
  //exitPerlInterpreter();
  
  return 0;
}

void report_stats(struct timeval time1, struct timeval time2) {
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;

  cerr << "\r";
  fprintf(stderr, 
          "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
          (unsigned long long) total_sequences, total_bases / 1.0e6, seconds,
          total_sequences / 1.0e3 / (seconds / 60),
          total_bases / 1.0e6 / (seconds / 60) );
  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) (total_sequences - total_classified),
          (total_sequences - total_classified) * 100.0 / total_sequences);
}

void process_file(char *filename) {
	string file_str(filename);
	DNASequenceReader *reader;
	DNASequence dna;
	unordered_map<string, unordered_map<uint64_t, int>> genusToKmerCounts;
	const char * genus_divergence_file = "divergence_file.div";
	ofstream out;
	out.open(genus_divergence_file);
	out << "Genus,KLdivergence\n";
  
	if (Fastq_input)
		reader = new FastqReader(file_str);
	else
		reader = new FastaReader(file_str);
		
	#pragma omp parallel
	{
		vector<DNASequence> work_unit;
		ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;
		
		while (reader->is_valid()) {
			work_unit.clear();
			size_t total_nt = 0;
			#pragma omp critical(get_input)
			{
				while (total_nt < Work_unit_size) {
					dna = reader->next_sequence();
					if (! reader->is_valid())
						break;
					work_unit.push_back(dna);
					total_nt += dna.seq.size();
				}
			}
			if (total_nt == 0)
				break;
		
			kraken_output_ss.str("");
			classified_output_ss.str("");
			unclassified_output_ss.str("");
		
			for (size_t j = 0; j < work_unit.size(); j++) {
				//cout << "Got to sequence " << j << endl;
				unordered_map<string, unordered_map<uint64_t, int>> currentGenusToKmerCounts = 
								classify_sequence( work_unit[j], kraken_output_ss,
								classified_output_ss, unclassified_output_ss);
				//better: return outcome of classify_sequence as a struct of uint32_t and unordered_map<uint64_t, int>
				if (!currentGenusToKmerCounts.empty()) {
					for (auto &currentGenusKmerCount : currentGenusToKmerCounts) {
						string currentGenus = currentGenusKmerCount.first;
						unordered_map<uint64_t, int> currentKmerCounts = currentGenusKmerCount.second;
						//if current tax id not in the overall taxIds map:
						if (genusToKmerCounts.find(currentGenus) == genusToKmerCounts.end()) {
							genusToKmerCounts[currentGenus] = currentKmerCounts;
						}
						//if it's already there: 
						else {
							unordered_map<uint64_t, int> storedKmerCounts = genusToKmerCounts[currentGenus];
							for (auto &kmerCounts: currentKmerCounts) {
								uint64_t kmerIndex = kmerCounts.first;
								int kmerCount = kmerCounts.second;
								//if kmer not stored put it in the map:
								if (storedKmerCounts.find(kmerIndex) == storedKmerCounts.end()) {
									storedKmerCounts[kmerIndex] = kmerCount;
								}
								//else increase the count with the current count:
								else {
									storedKmerCounts[kmerIndex] += kmerCount;
								}
							}
							genusToKmerCounts[currentGenus] = storedKmerCounts;
						}
					}
				}
			}
			#pragma omp critical(write_output)
			{
				if (Print_kraken)
					(*Kraken_output) << kraken_output_ss.str();
				if (Print_classified)
					(*Classified_output) << classified_output_ss.str();
				if (Print_unclassified)
					(*Unclassified_output) << unclassified_output_ss.str();
					total_sequences += work_unit.size();
					total_bases += total_nt;
					cerr << "\rProcessed " << total_sequences << " sequences (" << total_bases << " bp) ...";
				}
			}
		}  // end parallel section
	//cout << endl;
	delete reader;
	
	//go through genusToKmerCounts and compute the KL divergence for each taxID/genus:
	if (!genusToKmerCounts.empty()) {
		for (auto &genusKmerCountsPair : genusToKmerCounts) {
			//cout << genusKmerCountsPair.first << endl;
			unordered_map<uint64_t, int> kmerCounts = genusToKmerCounts[genusKmerCountsPair.first];
			int sum = accumulate(begin(kmerCounts), end(kmerCounts), 0, [](int prev, const pair<uint64_t,int>& p) { return prev+p.second; });
			//cout << "Size: " << kmerCounts.size() << endl;	
			double expectedRelFrequency = 1.0/kmerCounts.size();
			//cout << "Exp rel freq: " << expectedRelFrequency << endl;
			double KLdivergence = 0.0;
			for (auto &kmerCount : kmerCounts) {
				double kmerRelFrequency = (double)kmerCount.second/sum;
				//cout << kmerCount.first << " : " << kmerRelFrequency << endl;
				KLdivergence += kmerRelFrequency * log(kmerRelFrequency/expectedRelFrequency);
			}
			//cout << "KL divergence for " << genusKmerCountsPair.first << " : " << KLdivergence << endl;
			out << genusKmerCountsPair.first; 
			out << ","; 
			out << KLdivergence;
			out << "\n";
		}
	}
	out.close();
}

unordered_map<string, unordered_map<uint64_t, int>> classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss) {
	vector<uint32_t> taxa;
	vector<uint8_t> ambig_list;
	map<uint32_t, uint32_t> hit_counts;
	uint64_t *kmer_ptr;
	uint32_t taxon = 0;
	uint32_t hits = 0;  // only maintained if in quick mode
	
	uint64_t current_bin_key;
	int64_t current_min_pos = 1;
	int64_t current_max_pos = 0;
	
	//vector<uint64_t> kmers;
	//string dnaHeader = "";
	
	unordered_map<string, unordered_map<uint64_t, int>> genusToKmerCount;
	unordered_map<uint64_t, int> kmerIndexCount;
	
	int count = 0;
	//number of kmers in a sequence: seqlength - kmer size + 1 
	//def kmer size for kraken = 31
	int totalKmers = dna.seq.size() - 31 + 1;
	
	if (dna.seq.size() >= Database.get_k()) {
		//dnaHeader = dna.header_line;
		KmerScanner scanner(dna.seq);
		
		while ((kmer_ptr = scanner.next_kmer()) != NULL) {
			taxon = 0;
			if (scanner.ambig_kmer()) {
				ambig_list.push_back(1);
			}
			else {
				ambig_list.push_back(0);
				uint32_t *val_ptr = Database.kmer_query(
									Database.canonical_representation(*kmer_ptr),
									&current_bin_key,
									&current_min_pos, &current_max_pos
									);
				//koss << Database.canonical_representation(*kmer_ptr) << " ";
				taxon = val_ptr ? *val_ptr : 0;
				if (taxon != 0 && taxon != 1 && taxon != 2) {
					//cout << "Assigned kmer to taxon " << taxon << endl;
					count += 1;
					uint64_t kmerIndex = Database.canonical_representation(*kmer_ptr);
					//kmers.push_back(Database.canonical_representation(*kmer_ptr));
					if (kmerIndexCount.find(kmerIndex) == kmerIndexCount.end()) {
						kmerIndexCount.insert(pair<uint64_t, int>(kmerIndex, 1));
					}
					else {
						kmerIndexCount[kmerIndex] += 1;
					}
					hit_counts[taxon]++;
					if (Quick_mode && ++hits >= Minimum_hit_count)
						break;
				}
			}
			taxa.push_back(taxon);
		}
	}
		
	uint32_t call = 0;
	if (Quick_mode)
		call = hits >= Minimum_hit_count ? taxon : 0;
	else {
		call = resolve_tree(hit_counts, Parent_map);
		//
	}
	if (call != 0 && call != 1 && call != 2 && count >= totalKmers/2){
		#pragma omp atomic
		total_classified++;
		cout << "Sequence classified to tax ID " << call << endl;
		if (total_classified % 100 == 0) {
			cout << "Classified " << total_classified << endl;
		}
	}	
	if (Print_unclassified || Print_classified) {
		ostringstream *oss_ptr = call ? &coss : &uoss;
		bool print = call ? Print_classified : Print_unclassified;
		if (print) {
			if (Fastq_input) {
				(*oss_ptr) << "@" << dna.header_line << endl
					<< dna.seq << endl
					<< "+" << endl
					<< dna.quals << endl;
			}
			else {
				(*oss_ptr) << ">" << dna.header_line << endl
					<< dna.seq << endl;
			}
		}
	}

	if (! Print_kraken)
		exit(1);

	if (call != 0 && call != 1 && call != 2 && count >= totalKmers/2) {
		koss << endl;
		koss << "C\t";
	}
	
	else {
		if (Only_classified_kraken_output)
			exit(1); //return;
		koss << endl;
		koss << "U\t";
	}
	koss << dna.id << "\t" << dna.seq.size() << "\t";

	if (Quick_mode) {
		koss << "Q:" << hits;
	}
 
	else {
		if (taxa.empty())
			koss << "0:0";
		else
			koss << hitlist_string(taxa, ambig_list);
	}
	
	if (call != 0 && call != 1 && call != 2 && count >= totalKmers/2) {
		//koss << endl << dnaHeader << endl;
		//convert call (taxid) to genus prior to insertion into the map:
		//string genus = getGenusName(call);
		
		
		//move inside if!!!!
		string taxID;
		ostringstream o;
		o << call;
		taxID += o.str();
		
		
		if (taxIDtoGenus.find(call) == taxIDtoGenus.end()) {
			int random = rand();
			string randomStr;
			ostringstream convert;
			convert << random; 
			randomStr = convert.str();
			string genus_file = "/home/avoicu/genusFiles/genus_" + randomStr + ".txt";
			
			char command[200];
			strcpy(command, "/home/avoicu/localperl/bin/perl /home/avoicu/kraken-orig/getGenus.pl ");
			strcat(command, taxID.c_str());
			strcat(command, " ");
			strcat(command, genus_file.c_str());
			
			//cout << endl;
			//return int
			system(command);
			//string genus;
			ifstream infile(genus_file);
			string genus;
			getline(infile, genus);
			cout << "Genus: " << genus << endl;
			if (genus.compare("undef") != 0 && genus.compare("") != 0) {
				genusToKmerCount[genus] = kmerIndexCount;
				taxIDtoGenus.insert(pair<uint32_t, string>(call, genus));
				//cout << "Identified genus: " << genus << endl;
			}	
			
		}
		else {
			string genus = taxIDtoGenus[call];
			genusToKmerCount[genus] = kmerIndexCount;
		}
		/*infile.seekg(0, infile.end);
		int length = infile.tellg();
		cout << "Length " << length << endl;
		if (length != 0) {
			string line;
				getline(infile, line);
				cout << line << endl;
				string genus = line;
			}
			infile.close();
		//note: close file then delete it using system(rm..);
			
		//}
		
		/*	
		if (!kmers.empty()) {
			copy(kmers.begin(), kmers.end()-1, ostream_iterator<uint64_t>(koss, " "));
			//koss << "HERE: " << kmers.back() << endl;
			//cout << kmers.back() << endl;
		}
		*/
	}
	
	if (!genusToKmerCount.empty()) {
		/*for (auto &genusKmerCountPair : genusToKmerCount) {
			cout << genusKmerCountPair.first << endl;
			unordered_map<uint64_t, int> kmerCounts = genusToKmerCount[genusKmerCountPair.first];
			for (auto &kmerCount : kmerCounts) {
				cout << kmerCount.first << " : " << kmerCount.second << endl;
				 
			}
		}*/
		return genusToKmerCount;	
	}
	else {
		return unordered_map<string, unordered_map<uint64_t, int>>();
	}
}


string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig) {
	int64_t last_code;
	int code_count = 1;
	ostringstream hitlist;

	if (ambig[0])   { last_code = -1; }
	else            { last_code = taxa[0]; }

  for (size_t i = 1; i < taxa.size(); i++) {
    int64_t code;
    if (ambig[i]) { code = -1; }
    else          { code = taxa[i]; }

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code >= 0) {
    hitlist << last_code << ":" << code_count;
  }
  else {
    hitlist << "A:" << code_count;
  }
  return hitlist.str();
}

set<uint32_t> get_ancestry(uint32_t taxon) {
  set<uint32_t> path;

  while (taxon > 0) {
    path.insert(taxon);
    taxon = Parent_map[taxon];
  }
  return path;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  int sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qfcC:U:M")) != -1) {
    switch (opt) {
      case 'd' :
        DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 't' :
        sig = atoi(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive thread count");
        #ifdef _OPENMP
        Num_threads = sig;
        omp_set_num_threads(Num_threads);
        #endif
        break;
      case 'n' :
        Nodes_filename = optarg;
        break;
      case 'q' :
        Quick_mode = true;
        break;
      case 'm' :
        sig = atoi(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive minimum hit count");
        Minimum_hit_count = sig;
        break;
      case 'f' :
        Fastq_input = true;
        break;
      case 'c' :
        Only_classified_kraken_output = true;
        break;
      case 'C' :
        Print_classified = true;
        Classified_output_file = optarg;
        break;
      case 'U' :
        Print_unclassified = true;
        Unclassified_output_file = optarg;
        break;
      case 'o' :
        Kraken_output_file = optarg;
        break;
      case 'u' :
        sig = atoi(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive work unit size");
        Work_unit_size = sig;
        break;
      case 'M' :
        Populate_memory = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filename.empty()) {
    cerr << "Missing mandatory option -d" << endl;
    usage();
  }
  if (Index_filename.empty()) {
    cerr << "Missing mandatory option -i" << endl;
    usage();
  }
  if (Nodes_filename.empty() && ! Quick_mode) {
    cerr << "Must specify one of -q or -n" << endl;
    usage();
  }
  if (optind == argc) {
    cerr << "No sequence data files specified" << endl;
  }
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -n filename      NCBI Taxonomy nodes file" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -f               Input is in FASTQ format" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -h               Print this message" << endl
       << endl
       << "At least one FASTA or FASTQ file must be specified." << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}
