#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>

#include "FFT.h"

using namespace std;


vector<double> multiplyFFT(const vector<double>& a, const vector<double>& b) {

	vector<complex<double>> fa(a.begin(), a.end());
	vector<complex<double>> fb(b.begin(), b.end());

	int n = std::pow(2, std::ceil(std::log2(a.size() + b.size())));

	// resize to the nearest power of 2
	fa.resize(n);
	fb.resize(n);

	FFT(fa);
	FFT(fb);

	// perform element-wise multiplication in the frequency domain
	for (int i = 0; i < n; ++i) {
		fa[i] *= fb[i];
	}

	inverseFFT(fa);
	vector<double> ans(n);
	for (int i = 0; i < n; ++i) {
		ans[i] = fa[i].real();
	}

	return ans;
}

vector<double> multiplyNaive(const vector<double>& a, const vector<double>& b) {

	size_t n = std::max(a.size(), b.size());

	vector<double> ans((n - 1) + n, 0.0);

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			ans[i + j] += a[i] * b[j];
		}
	}

	return ans;
}


int main() {

	vector<double> a = { 3, 2, 1, 5 };
	vector<double> b = { 0, 0, 5, 2 };


	vector<double> resultNaive = multiplyNaive(a, b);
	vector<double> resultFFT = multiplyFFT(a, b);

	cout << "Naive: \n";
	for (size_t i = 0; i < resultNaive.size(); ++i) {
		cout << resultNaive[i] << ", ";
	}
	cout << "\n\nFFT: ";

	for (size_t i = 0; i < resultFFT.size(); ++i) {
		cout << resultFFT[i] << ", ";
	}
	cout << "\n";

	return 0;
}
