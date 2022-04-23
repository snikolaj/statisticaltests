#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <bitset>
#include <intrin.h>

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

int referencePopcount(uint8_t num) {
    return _mm_popcnt_u32(num);
}

unsigned int popcount(int num) {
    int retVal = 0;
    for (uint8_t i = 0; i < 8; i++) {
        retVal += num / (2 << i);
    }
    return num - retVal;
}

chiSquareResult chiSquareGOFbyte(vector<unsigned char>& buffer) {
    // returns critical value and degrees of freedom
    chiSquareResult result;
    double occurrences[OCCURRENCES_SIZE];
    double expected[OCCURRENCES_SIZE];

    // initialize variables 
    result.criticalValue = 0.0;
    result.df = OCCURRENCES_SIZE - 1;
    for (unsigned int i = 0; i < OCCURRENCES_SIZE; i++) {
        occurrences[i] = 0;
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

chiSquareResult chiSquareGOFbit(vector<unsigned char>& buffer) {
    chiSquareResult result;
    // 9 possible results for 8 bit population count - 0 to 8
    double occurrences[9];
    double expected[9];

    // df = 9 - 1
    result.criticalValue = 0.0;
    result.df = 8;
    
    // zero out arrays, memset just didn't want to work
    for (unsigned int i = 0; i < 9; i++) {
        expected[i] = 0;
        occurrences[i] = 0;
    }

    // get frequency of population counts given a uniform distribution of 8 bit integers
    for (unsigned int i = 0; i < 256; i++) {
        expected[referencePopcount(i)] += 1.0;
    }

    // count occurrence of all possible population counts in sample
    for (unsigned int i = 0; i < buffer.size(); i++) {
        occurrences[referencePopcount(buffer[i])] += 1.0;
    }

    // normalize (?) the counts
    for (int i = 0; i < 9; i++) {
        occurrences[i] /= (buffer.size() / 256);
    }

    // chi square
    for (unsigned int i = 0; i < 9; i++) {
        double chi = occurrences[i] - expected[i];
        result.criticalValue += ((chi * chi) / expected[i]);
    }

    return result;
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
    chiSquareResult chisquarebyte = chiSquareGOFbyte(buffer);
    chiSquareResult chisquarebit = chiSquareGOFbit(buffer);

    cout << "Name:                       " << name << endl;
    cout << "Size of data:               " << buffer.size() << endl;
    cout << "Entropy:                    " << shannonEntropy(buffer) << " (expected 8)" << endl;
    cout << "Mean:                       " << mean(buffer) << " (expected 127.5)" << endl;
    cout << "X^2 critical value (byte):  " << chisquarebyte.criticalValue << " df=" << chisquarebyte.df << endl;
    cout << "X^2 critical value (bit):   " << chisquarebit.criticalValue << " df=" << chisquarebit.df << endl;
    cout << endl;
}

int main()
{
    const unsigned int bufsize = 4096;
    cout << "===== Random Number Test Suite - Stefan =====\n\n";
    string name = "Random.org Dataset";
    vector<unsigned char> fixedbuf(bufsize);
    vector<unsigned char> buffer = openFileAndReturnBuffer("RandomNumbers");
    for (unsigned int i = 0; i < fixedbuf.size(); i++) {
        fixedbuf[i] = buffer[i];
    }
    testSuite(fixedbuf, name);

    name = "Unconditioned Test Data";
    buffer = openFileAndReturnBuffer();
    for (unsigned int i = 0; i < fixedbuf.size(); i++) {
        fixedbuf[i] = buffer[i];
    }
    testSuite(fixedbuf, name);

}
