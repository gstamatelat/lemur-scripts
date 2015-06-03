#include <iostream>
#include "indri/UnparsedDocument.hpp"
#include "indri/TextTokenizer.hpp"
#include "indri/TokenizedDocument.hpp"

int main(int argc, char * argv[])
{
        std::string s("some test i'm james $word$");

        indri::parse::UnparsedDocument *doc = new indri::parse::UnparsedDocument();
        doc->content = s.c_str();
        doc->contentLength = s.length();

        indri::parse::TextTokenizer tok(false, true);
        indri::parse::TokenizedDocument *tdoc = tok.tokenize(doc);

        indri::utility::greedy_vector<char *> terms = tdoc->terms;
        for (int i = 0; i < terms.size(); i++) {
                std::cout << terms[i] << " ";
        }

        return 0;
}
