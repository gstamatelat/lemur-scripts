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
int model = 0;
unsigned int display = 0;

/* Globals */
QueryEnvironment qe;
vector<DOCID_T>::size_type docsSize = -1;

/* Class Term */
class Term {
public:
	TERMID_T id;
	vector<DOCID_T> docs;
public:
	Term()
	{
		id = 0;
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
	return (a.docs.size() > b.docs.size());
}

void getDocuments(string query0, vector<DOCID_T> &docids, int n)
{
	vector<ScoredExtentResult> results = qe.runQuery(query0, n);

        for (vector<ScoredExtentResult>::size_type i = 0; i < results.size(); i++)
        {
                docids.push_back(results.at(i).document);
        }
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
		cout << "       --beta:    sets the beta parameter for Fbeta sorting, defaults to 1" << endl;
		cout << "       --model:   sets the multi-term model for the initial query, available values: 0 for indri model and 1 for gDF, defaults to 0, if set to gDF, query may only be terms seperated by space like --query=\"human rights\"" << endl;
		cout << "       --display: what to display: binary flags as integer, cross docs, defaults to 0 (don't show)" << endl;
		//cout << "\033[0m";
		return 0;
	}

	/* Read parameters */
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]).find("--index=") == 0) {
			index0 = string(argv[i]).substr(string(argv[i]).find("=")+1);
			qe.addIndex(index0);
		} else if (string(argv[i]).find("--query=") == 0) {
			query = string(argv[i]).substr(string(argv[i]).find("=")+1);
		} else if (string(argv[i]).find("--print=") == 0) {
			print = atoi(string(argv[i]).substr(string(argv[i]).find("=")+1).c_str());
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
	getDocuments(query, docids, 50);
	docsSize = docids.size();

	/* Get intersections (for recall) */
	if (verbose) { cout << "Documents count: " << docsSize << " -- 0%"; }
	bool(*fn_pt)(TERMID_T,TERMID_T) = TermIdCmp;
	map<TERMID_T,Term,bool(*)(TERMID_T,TERMID_T)> terms(fn_pt);
	for (vector<DOCID_T>::size_type i = 0; i < docsSize; i++)
	{
		const TermList * tmpTermList = ind.termList(docids.at(i));
		greedy_vector<TERMID_T> tt = tmpTermList->terms();
		sort(tt.begin(), tt. end());
		tt.erase(unique(tt.begin(), tt.end()), tt.end());
		for (size_t j = 0; j < tt.size(); j++)
		{
			TERMID_T tempId = tt.at(j);
			Term tmp; tmp.id = tempId;
			terms.insert( pair<TERMID_T,Term>(tempId,tmp) );
		}
		delete tmpTermList;
	}

	/* Convert map to vector*/
	vector<Term> terms_v;
	map<TERMID_T,Term,bool(*)(TERMID_T,TERMID_T)>::iterator it;
	for (it = terms.begin(); it != terms.end(); it++)
	{
		terms_v.push_back((*it).second);
	}

	/* Remove the empty term and the initial */
	for (vector<Term>::size_type i = 0; i < terms_v.size(); i++)
        {
                if (ind.term(terms_v.at(i).id).empty() || ind.term(terms_v.at(i).id) == query) {
			terms_v.erase(terms_v.begin()+i);
			i--;
			cout << "removed" << endl;
		}
        }

	//cout << "Terms size: " << terms_v.size() << endl;

	/* Get 1000 first documents of each term */
	for (vector<Term>::size_type i = 0; i < terms_v.size(); i++)
	{
		//cout << ind.term(terms_v.at(i).id) << " "; cout.flush();
		//cout << i << " "; cout.flush();
		getDocuments(ind.term(terms_v.at(i).id), terms_v.at(i).docs, 1000);
	}
	//cout << endl;

	/* Remove unrelevant documents */
	for (vector<Term>::size_type j = 0; j < terms_v.size(); j++)
	{
		vector<DOCID_T> tmpVector;
		for (vector<DOCID_T>::size_type i = 0; i < docids.size(); i++)
		{
			if (find(terms_v.at(j).docs.begin(), terms_v.at(j).docs.end(), docids.at(i)) != terms_v.at(j).docs.end())
			{
				tmpVector.push_back(docids.at(i));
			}
		}
		terms_v.at(j).docs = tmpVector;
	}

	/* Sort */
	sort(terms_v.begin(), terms_v.end(), TermFbetaCmp);

	/* Check if full coverage with the first 50 terms */
	set<DOCID_T> docset;
	for (vector<Term>::size_type j = 0; j < 10; j++)
	{
		for (vector<DOCID_T>::size_type i = 0; i < terms_v.at(j).docs.size(); i++)
		{
			docset.insert(terms_v.at(j).docs.at(i));
		}
	}

	/* Print coverage */
	cout << docset.size() << "/" << docsSize << endl;

	/* Sort and print top results */
	/*sort(terms_v.begin(), terms_v.end(), TermFbetaCmp);
	for (vector<Term>::size_type i = 0; ((vector<Term>::size_type)print < 0 || i < (vector<Term>::size_type)print) && i < terms_v.size(); i++)
	{
		cout << std::fixed << terms_v.at(i).docs.size() << " " << ind.term(terms_v.at(i).id) << endl;
	}*/

	/* Ok, close the index */
	ind.close();
	qe.close();

	return 0;
}
