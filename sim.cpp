#include <iostream>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>
#include <omp.h>
#include <math.h>
#include "indri/QueryEnvironment.hpp"
#include "indri/KrovetzStemmer.hpp"

using namespace std;
using namespace lemur::api;
using namespace indri::api;
using namespace indri::parse;
using namespace indri::utility;
using namespace indri::index;

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

/* Parameters */
string index0;
string term;

class Term {
public:
        string name;
	double distance;
        double results;
public:
        Term()
        {
        }
        ~Term() { }
};

/* Helper functions */
int termCmp(Term a, Term b)
{
	// Maybe we should take into account results too ...
	return a.distance > b.distance;
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

string getCommonCharacters(string string1, string string2, int allowedDistance)
{
	int str1_len = string1.length();
	int str2_len = string2.length();
	string temp_string2 = string2;

	string commonCharacters = "";
	for (int i = 0; i < str1_len; i++)
	{
		bool noMatch = true;

		// compare if char does match inside given allowedDistance
		// and if it does add it to commonCharacters
		for (int j = max(0, i - allowedDistance); noMatch && j < min(i + allowedDistance + 1, str2_len ); j++)
		{
			if(temp_string2[j] == string1[i])
			{
				noMatch = false;
				commonCharacters.append(1,string1[i]);
				temp_string2[j] = '?';
				//temp_string2.erase(j, 1);
			}
		}
	}
	return commonCharacters;
}

double Jaro(string string1, string string2)
{
	int str1_len = string1.length();
	int str2_len = string2.length();

	// theoretical distance
	int distance = floor(min( str1_len, str2_len ) / 2.0);

	// get common characters
	string commons1 = getCommonCharacters(string1, string2, distance);
	string commons2 = getCommonCharacters(string2, string1, distance );

	//cout << commons1 << endl << commons2 << endl;

	int commons1_len; int commons2_len;
	if( (commons1_len = commons1.length()) == 0) return 0;
	if( (commons2_len = commons2.length()) == 0) return 0;

	// calculate transpositions
	double transpositions = 0;
	int upperBound = min( commons1_len, commons2_len );
	for(int i = 0; i < upperBound; i++)
	{
		if( commons1[i] != commons2[i] )
			transpositions++;
	}
	transpositions /= 2.0;

	// return the Jaro distance
	return (commons1_len/(double)str1_len + commons2_len/(double)str2_len + (commons1_len - transpositions)/(double)commons1_len) / 3.0;
}

int getPrefixLength(string string1, string string2, int MINPREFIXLENGTH = 4 )
{
	int n = min( MINPREFIXLENGTH, min (string1.length(), string2.length()) );

	for(int i = 0; i < n; i++)
	{
		if( string1[i] != string2[i] )
		{
			// return index of first occurrence of different characters
			return i;
		}
	}

	// first n characters are the same
	return n;
}

double JaroWinkler(string string1, string string2, double PREFIXSCALE = 0.1 ){

  double JaroDistance = Jaro( string1, string2 );

  int prefixLength = getPrefixLength( string1, string2 );

  return JaroDistance + prefixLength * PREFIXSCALE * (1.0 - JaroDistance);
}

int main(int argc, char * argv[])
{
	/* Print help */
	if (argc == 1)
	{
		//cout << "\033[0;31m";
		cout << "Usage: sim --index=/path/to/index/ --term=term" << endl;
		cout << "Notes: --index: location MUST include the trailing slash" << endl;
		cout << "       --term:  term" << endl;
		//cout << "\033[0m";
		return 0;
	}

	/* Read parameters */
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]).find("--index=") == 0) {
			index0 = string(argv[i]).substr(string(argv[i]).find("=")+1);
		} else if (string(argv[i]).find("--term=") == 0) {
			term = string(argv[i]).substr(string(argv[i]).find("=")+1);
		} else {
			cout << "Unrecognised parameter " << string(argv[i]) << ". Execute ./sim for help" << endl;
			return 0;
		}
	}

	/* Open index so that it's ready when we need it */
	string indexLocation = index0 + "index/";
	indexLocation = indexLocation + getdir(indexLocation) + "/";
	DiskIndex ind; ind.open(indexLocation, "");

	/* Gather similar terms */
	// Maybe use frequentVocabularyIterator with a non stemmed index
	vector<Term> terms;
	indri::index::VocabularyIterator* iter = ind.vocabularyIterator();
	iter->startIteration();
	while(!iter->finished())
	{
		indri::index::DiskTermData* entry = iter->currentEntry();
		indri::index::TermData* termData = entry->termData;
		Term tmpTerm;
		tmpTerm.name = termData->term;
		tmpTerm.distance = JaroWinkler(termData->term, term);
		tmpTerm.results = (double)termData->corpus.documentCount / (double)ind.documentCount();
		terms.push_back(tmpTerm);
		iter->nextEntry();
	}
	delete iter;
	ind.close();

	/* Sort and print */
	sort(terms.begin(), terms.end(), termCmp);
	for (int i = 0; i < 10; i++)
	{
		cout << terms.at(i).distance << " " << terms.at(i).name << " " << terms.at(i).results << endl;
	}

	return 0;
}
