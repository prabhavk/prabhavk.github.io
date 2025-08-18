#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>  
#include <numeric> 
#include <random>
#include <sstream>
#include <stdio.h>
#include <string>
#include <tuple>
#include <vector>


class MST_vertex;
class MST;
class SEM_vertex;
class clique;
class cliqueTree;
class SEM;

using namespace std;

struct mt_error : runtime_error {
    using runtime_error::runtime_error;
};


class EMManager
{
private:
	default_random_engine generator;
	vector <string> sequenceNames;
	map <string,unsigned char> mapDNAtoInteger;		
	ofstream emt_logFile;
	int numberOfLargeEdgesThreshold;
	int numberOfHiddenVertices = 0;
	int edgeWeightThreshold;	
	chrono::steady_clock::time_point start_time;
	chrono::steady_clock::time_point current_time;
	chrono::steady_clock::time_point t_start_time;
	chrono::steady_clock::time_point t_end_time;
	chrono::steady_clock::time_point m_start_time;
	chrono::steady_clock::time_point m_end_time;
	chrono::duration<double> timeTakenToComputeEdgeAndVertexLogLikelihoods;
	chrono::duration<double> timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	chrono::duration<double> timeTakenToComputeSubtree;
	chrono::duration<double> timeTakenToComputeSupertree;
	chrono::duration<double> timeTakenToRootViaEdgeLoglikelihoods;
	chrono::duration<double> timeTakenToRootViaRestrictedSEM;
	string fastaFileName;
	string phylipFileName;
	string topologyFileName;
	string prefix_for_output_files;
	string ancestralSequencesString;
	string loglikelihood_node_rep_file_name;	
	string probabilityFileName_pars;
	string probabilityFileName_diri;
	string probabilityFileName_pars_root;
	string probabilityFileName_diri_root;
	string probabilityFileName_best;
	double max_ll_pars;
	double max_ll_diri;
	string MSTFileName;
	string GMMparametersFileName;
	string distance_measure_for_NJ = "Hamming";
	bool apply_patch = false;
	bool grow_tree_incrementally = false;
	bool flag_topology = false;
    bool flag_set_gmm_parameters = false;
	int ComputeHammingDistance(string seq1, string seq2);
	int ComputeHammingDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2);
	int GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices);
	MST * M;
	SEM * P;
	SEM * p;
	void WriteOutputFiles();
	bool debug;
	bool verbose;
	bool localPhyloOnly;	
	bool useChowLiu;
	bool modelSelection; 	
    int max_iter;		
	string supertree_method;
	int numberOfVerticesInSubtree;
	string GetSequenceListToWriteToFile(map <string, vector <unsigned char>> compressedSeqMap, vector <vector <int> > sitePatternRepetitions);
	vector <string> must_have;
	vector <string> may_have;
    int num_repetitions;
	int max_EM_iter;
	double conv_thresh;
    public:
    EMManager(string sequenceFileNameToSet,
              string input_format,
              string topologyFileNameToSet,
              string prefix_for_output_files,
              int num_repetitions,
              int max_iter,
              double conv_threshold);
	~EMManager();
	double max_log_lik;
	double max_log_lik_pars;
	double max_log_lik_diri;
	double max_log_lik_ssh;
	void SetDNAMap();
	void SetThresholds();
	void EMTRackboneWithOneExternalVertex();
	void EMTRackbone_k2020_preprint();
	void EMgivenInputTopology();
	void RootSuperTree();		
	void EMTRackboneWithRootSEMAndMultipleExternalVertices();
	void EMTRackboneOverlappingSets();
	void EMTRackboneOnlyLocalPhylo();
	void EMforCompleteData();
	void ReadAllSequences(string complete_sequence_file_name);
	void main(string init_criterion, bool root_search);
	void EM_main();
    void EMparsimony();
	void EMdirichlet();
	void SetprobFileforSSH();
    void EMssh();
	string EncodeAsDNA(vector<unsigned char> sequence);
	vector<unsigned char> DecompressSequence(vector<unsigned char>* compressedSequence, vector<vector<int>>* sitePatternRepeats);	
};
