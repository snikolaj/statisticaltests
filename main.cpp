#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace std;
#define OCCURRENCES_SIZE 256

struct chiSquareResult {
    double criticalValue;
    unsigned int df;
};

vector<unsigned char> openFileAndReturnBuffer(string fileName = "actual_results.BIN") {
    ifstream input(fileName, std::ios::binary);
    vector<unsigned char> buffer(std::istreambuf_iterator<char>(input), {});
    return buffer;
}

double logBase(double a, double b) {
    return log10(a) / log10(b);
}

chiSquareResult chiSquareGOF(vector<unsigned char>& buffer) {
    // returns critical value and degrees of freedom
    chiSquareResult result;
    double occurrences[OCCURRENCES_SIZE];
    double expected[OCCURRENCES_SIZE];

    // initialize variables (yeah you shouldn't memset doubles, but with IEE-754 it's fine)
    result.criticalValue = 0.0;
    result.df = OCCURRENCES_SIZE - 1;
    memset(occurrences, 0, OCCURRENCES_SIZE * sizeof(double));
    for (unsigned int i = 0; i < OCCURRENCES_SIZE; i++) {
        expected[i] = buffer.size() / OCCURRENCES_SIZE;
    }

    for (unsigned int i = 0; i < buffer.size(); i++) {
        occurrences[buffer[i]] += 1.0;
    }

    for (unsigned int i = 0; i < OCCURRENCES_SIZE; i++) {
        double chi = occurrences[i] - expected[i];
        result.criticalValue += ((chi * chi) / (double)expected[i]);
    }
    
    return result;
    // x^2 GOF test https://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm
    // https://www.codeproject.com/Articles/432194/How-to-Calculate-the-Chi-Squared-P-Value
}

double shannonEntropy(vector<unsigned char>& buffer) {
    double finalEntropy = 0.0;
    
    double occurrences[OCCURRENCES_SIZE];
    memset(occurrences, 0, OCCURRENCES_SIZE * sizeof(double));
    for (unsigned int i = 0; i < buffer.size(); i++) {
        occurrences[buffer[i]] += 1.0;
    }

    /*for (unsigned int i = 0; i < OCCURRENCES_SIZE; i++) {
        cout << i << ": " << occurrences[i] << endl;
    }*/

    for (unsigned int i = 0; i < OCCURRENCES_SIZE; i++) {
        double intermediate = (occurrences[i] / buffer.size()) * log2(occurrences[i] / (double)buffer.size());
        finalEntropy -= (isnan(intermediate)) ? 0.0 : intermediate;
    }
    return finalEntropy;

    // inspired by http://rosettacode.org/wiki/Entropy
}

double mean(vector<unsigned char>& buffer) {
    double meanValue = 0.0;
    for (unsigned int i = 0; i < buffer.size(); i++) {
        meanValue += (double)buffer[i];
    }
    return meanValue / (double)buffer.size();
}

void testSuite(vector<unsigned char>& buffer, string name) {
    chiSquareResult temp = chiSquareGOF(buffer);
    cout << "Name:                " << name << endl;
    cout << "Size of data:        " << buffer.size() << endl;
    cout << "Entropy:             " << shannonEntropy(buffer) << " (expected 8)" << endl;
    cout << "Mean:                " << mean(buffer) << " (expected 128)" << endl;
    cout << "X^2 critical value:  " << temp.criticalValue << " df=" << temp.df << endl;
    cout << endl;
}

int main()
{
    cout << "===== Random Number Test Suite - Stefan =====\n";
    cout << "=====          it's alright!            =====\n\n";
    string name = "Random.org Dataset";
    vector<unsigned char> buffer = openFileAndReturnBuffer("RandomNumbers");
    testSuite(buffer, name);

    name = "Unconditioned Test Data";
    buffer = openFileAndReturnBuffer();
    testSuite(buffer, name);
}
