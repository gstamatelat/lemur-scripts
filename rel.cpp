#include <iostream>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>
#include <omp.h>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include "indri/QueryEnvironment.hpp"
#include "indri/KrovetzStemmer.hpp"

using namespace std;
using namespace lemur::api;
using namespace indri::api;
using namespace indri::parse;
using namespace indri::utility;
using namespace indri::index;

/* Parameters */
bool verbose = false;
string index0;
string query;
int prec = 4;
int print = -1;
double beta = 1;
double k = 1.0;
int model = 0;
unsigned int display = 0;

/* Globals */
vector<DOCID_T>::size_type docsSize = -1;
long documentsCount = -1;

/* Class Term */
class Term {
public:
	TERMID_T id;
	UINT64 cross;
	UINT64 results;
	vector<DOCID_T> docs;
public:
	double Precision()
	{
		double prec = (double)cross/(double)results;
		return prec;
	}
	double Fbeta()
	{
		double rec = (double)cross/(double)docsSize;
		double prec = (double)cross/(double)results;
		return (1+pow(beta,2))*(prec*rec)/((pow(beta,2))*prec+rec);
	}
	double MI()
	{
		double M = (double)documentsCount;
		double N11 = (double)cross;
		double N10 = (double)results - (double)cross; // can be zero
		double N01 = (double)docsSize - (double)cross; // can be zero
		double N00 = M - (double)docsSize - N10;
		double N1_ = N10 + N11;
		double N_1 = N01 + N11;
		double N0_ = N00 + N01;
		double N_0 = N00 + N10;
		double pc1 = (double)docsSize / M;
		double pc0 = 1.0 - pc1;
		double Hc = - ( pc1*(log(pc1)/log(2)) + pc0*(log(pc0)/log(2)) );
		double MI = (N11/M) * (log((M * N11)/(N1_ * N_1))/log(2));
		if(N01 > 0) {
			MI += (N01/M) * (log((M * N01)/(N0_ * N_1))/log(2));
		}
		if(N10 > 0) {
			MI += (N10/M) * (log((M * N10)/(N1_ * N_0))/log(2));
		}
		MI += (N00/M) * (log((M * N00)/(N0_ * N_0))/log(2));
		MI /= Hc;
		//if (MI != MI) { MI = 0; }
		return MI;
	}
public:
	Term()
	{
		id = 0;
		cross = 0;
		results = 0;
	}
	~Term() { }
};

/* Helper functions */
int num_digits(int x)
{
	return x > 0 ? (int) log10 ((double) x) + 1 : 1;
}

string getdir(string dir)
{
	string out("");
	DIR *dp;
	struct dirent *dirp;
	dp  = opendir(dir.c_str());
	while (out == "" || out == ".." || out == ".")
	{
		dirp = readdir(dp);
		out = string(dirp->d_name);
	}
	closedir(dp);
	return out;
}

bool TermIdCmp(TERMID_T lhs, TERMID_T rhs)
{
	return lhs < rhs;
}

bool TermFbetaCmp(Term a, Term b)
{
	//return (a.Fbeta() > b.Fbeta());
	return (a.MI() > b.MI());
}

void getDocuments(vector<DOCID_T> &docids)
{
	QueryEnvironment qe;
	qe.addIndex(index0);

	documentsCount = qe.documentCount();

        vector<ScoredExtentResult> results;

	if (model == 1)
	{
		// Preprocessing
		boost::algorithm::trim(query);
		while (std::string::npos != query.find("  ")) { boost::replace_all(query, "  ", " "); }
		std::vector<std::string> strs;
		boost::split(strs, query, boost::is_any_of(" "));

		// Get minDF
		INT64 minDF = qe.documentCount(strs.at(0));
		for (size_t i = 0; i < strs.size(); i++)
		{
			if (qe.documentCount(strs.at(i)) < minDF) {
				minDF = qe.documentCount(strs.at(i));
			}
		}
		// Major bug here
		//if (minDF == 0) { minDF = 1; }

		// Get andDF
		vector<ScoredExtentResult> and_results;
		string andquery = "#uw(" + query + ")";
		and_results = qe.runQuery(andquery, documentsCount);
		vector<ScoredExtentResult>::size_type andDF = and_results.size();
		if (andDF == 0) { andDF = 1; }

		// Get gDF
		double gDF = sqrt((double)andDF*(double)minDF);

		// Run query
		//string orquery = "#or(" + query + ")";
		string orquery = query;
		results = qe.runQuery(orquery, (int)ceil((double)gDF));
	}
	else
	{
		results = qe.runQuery(query, documentsCount);
	}

        for (vector<ScoredExtentResult>::size_type i = 0; i < results.size(); i++)
        {
                docids.push_back(results.at(i).document);
        }

	qe.close();
}

int main(int argc, char * argv[])
{
	/* Print help */
	if (argc == 1)
	{
		//cout << "\033[0;31m";
		cout << "Usage: rel [--display=display] [--verbose] [--print=count] [--prec=digits] [--beta=beta] [--model=model] --index=/path/to/index/ --query=query" << endl;
		cout << "Notes: --index:   location MUST include the trailing slash" << endl;
		cout << "       --query:   can be a valid indri language construct like --query=\"#1(human rights)\", see --model also" << endl;
		cout << "       --verbose: prints aditional information, defaults to no verbose" << endl;
		cout << "       --print:   prints top count terms, negative means all, 0 means no output, defaults to -1" << endl;
		cout << "       --prec:    sets the console precision for fp numbers, defaults to 4" << endl;
		cout << "       --k:       sets the k-anonymity value, defaults to 1.0" << endl;
		cout << "       --beta:    sets the beta parameter for Fbeta sorting, defaults to 1" << endl;
		cout << "       --model:   sets the multi-term model for the initial query, available values: 0 for indri model and 1 for gDF, defaults to 0, if set to gDF, query may only be terms seperated by space like --query=\"human rights\"" << endl;
		cout << "       --display: what to display: binary flags as integer, cross docs, defaults to 0 (don't show)" << endl;
		//cout << "\033[0m";
		return 0;
	}

	/* Read parameters */
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]) == "--verbose") {
			verbose = true;
		} else if (string(argv[i]).find("--index=") == 0) {
			index0 = string(argv[i]).substr(string(argv[i]).find("=")+1);
		} else if (string(argv[i]).find("--query=") == 0) {
			query = string(argv[i]).substr(string(argv[i]).find("=")+1);
		} else if (string(argv[i]).find("--print=") == 0) {
			print = atoi(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
		} else if (string(argv[i]).find("--prec=") == 0) {
			prec = atoi(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
		} else if (string(argv[i]).find("--k=") == 0) {
			k = atof(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
		} else if (string(argv[i]).find("--beta=") == 0) {
			beta = atof(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
		} else if (string(argv[i]).find("--model=") == 0) {
                        model = atoi(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
		} else if (string(argv[i]).find("--display=") == 0) {
			display = atoi(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
		} else {
			cout << "Unrecognised parameter " << string(argv[i]) << ". Execute ./rel for help" << endl;
			return 0;
		}
	}

	/* Open index so that it's ready when we need it */
	string indexLocation = index0 + "index/";
	indexLocation = indexLocation + getdir(indexLocation) + "/";
	DiskIndex ind; ind.open(indexLocation, "");

	/* Get OpenMP */
	if (verbose)
	{
		int nthreads; //TODO: How can number of threads by signed? ;0 Check return type to avoid casting
		#pragma omp parallel for
		for (int i = 0; i < 24; i++)
		{
			nthreads = omp_get_num_threads();
		}
		cout << "Parallel algorithms are using " << nthreads << " threads" << endl;
	}

	/* Get documents */
	vector<DOCID_T> docids;
	getDocuments(docids);
	docsSize = docids.size();

	/* Get intersections (for recall) */
	if (verbose) { cout << "Documents count: " << docsSize << " -- 0%"; }
	bool(*fn_pt)(TERMID_T,TERMID_T) = TermIdCmp;
	map<TERMID_T,Term,bool(*)(TERMID_T,TERMID_T)> terms(fn_pt);
	for (vector<DOCID_T>::size_type i = 0; i < docsSize; i++)
	{
		if (verbose)
		{
			//TODO: embarassing way. better is to convert to double to avoid size_type overflow
			for (int sp = 0; sp < num_digits(100*(i)/docsSize) + 1; sp++) { cout << "\b"; }
			cout << 100*(i+1)/docsSize << "%";
		}
		const TermList * tmpTermList = ind.termList(docids.at(i));
		greedy_vector<TERMID_T> tt = tmpTermList->terms();
		sort(tt.begin(), tt. end());
		tt.erase(unique(tt.begin(), tt.end()), tt.end());
		//Embarassing way of removing duplicates? Which is faster?
		//greedy_vector<TERMID_T>::iterator it = unique(tt.begin(), tt.end()); tt.resize(it - tt.begin());
		for (size_t j = 0; j < tt.size(); j++)
		{
			TERMID_T tempId = tt.at(j);
			Term tmp; tmp.id = tempId;
			terms.insert( pair<TERMID_T,Term>(tempId,tmp) );
			terms[tempId].cross++;
			if (display == 1) { terms[tempId].docs.push_back(docids.at(i)); }
		}
		delete tmpTermList;
	}
	if (verbose) { cout << endl; }

	/* Get results (for precision) */
	// It's probably much better here to use indri::index::VocabularyIterator* iter = ind.vocabularyIterator() rather than ind.term((*it).first)
	map<TERMID_T,Term,bool(*)(TERMID_T,TERMID_T)>::size_type termsSize = terms.size();
	map<TERMID_T,Term,bool(*)(TERMID_T,TERMID_T)>::iterator it;
	map<TERMID_T,Term,bool(*)(TERMID_T,TERMID_T)>::size_type iter = 1;
	if (verbose) { cout << "Terms count: " << termsSize << " -- 0%"; }
	for (it = terms.begin(); it != terms.end(); it++)
	{
		if (verbose)
                {
			//TODO: embarassing way. better is to convert to double to avoid size_type overflow
			for (int sp = 0; sp < num_digits(100*(iter-1)/termsSize) + 1; sp++) { cout << "\b"; }
			cout << 100*(iter)/termsSize << "%";
			iter++;
                }
		(*it).second.results = ind.documentCount(ind.term((*it).first));
	}
	if (verbose) { cout << endl; }

	/* Remove terms with precision higher than 1/k */
	if (k > 1.0) {
		double anonymity = 1.0 / k;
		it = terms.begin();
		while (it != terms.end()) {
			if ((*it).second.Precision() > anonymity) {
				terms.erase(it++);
			} else {
				it++;
			}
		}
	}

	/* Convert map to vector and sort */
	vector<Term> terms_v;
	for (it = terms.begin(); it != terms.end(); it++)
	{
		terms_v.push_back((*it).second);
	}
	sort(terms_v.begin(), terms_v.end(), TermFbetaCmp);

	/* Print results */
	for (vector<Term>::size_type i = 0; ((vector<Term>::size_type)print < 0 || i < (vector<Term>::size_type)print) && i < terms_v.size(); i++)
	{
		// TODO: Possible bug in extra space after sentence, remove it
		cout << std::fixed << std::setprecision(prec) << terms_v.at(i).MI() << " "
			<< ind.term(terms_v.at(i).id) << " "
			<< terms_v.at(i).cross << " "
			<< terms_v.at(i).results << " ";
		for (vector<DOCID_T>::size_type j = 0; j < terms_v.at(i).docs.size(); j++)
		{
			cout << terms_v.at(i).docs.at(j) << " ";
		}
		cout << endl;
	}

	/* Ok, close the index */
	ind.close();

	return 0;
}
