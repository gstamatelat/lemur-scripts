#include <iostream>
#include "indri/KrovetzStemmer.hpp"

#define INSERT_NEW_LINE true

int main(int argc, char * argv[])
{
	indri::parse::KrovetzStemmer ks;
	for (int i = 1; i < argc; i++)
	{
		if (i > 1) std::cout << " ";
		std::cout << ks.kstem_stemmer(argv[i]);
	}
	if (INSERT_NEW_LINE) std::cout << std::endl;

	return 0;
}
