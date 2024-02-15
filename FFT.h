#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

#define PI 3.1415926535

// https://cp-algorithms.com/algebra/fft.html

// FFT uses divide-and-conquer to efficiently compute the Discrete Fourier Transform (DFT)
// We can use FFT to multiply two polynomial in O(n * log(n)) time complexity. The idea is
// that in the frequency domain, polynomial multiplication can be done with element-wise
// multiplication, wich is O(n), so we aply FFT in both polynomials, perform element-wise
// multiplication and than go back to the time domain with the inverse FFT. 

// we use divide-and-conquer this way:

// A(x) = a_0 * x^0 + a_1 * x^1 + ... + a_n-1 * x^n-1 , n is a power of 2, n > 1

// we can decompose A(x) in two other polynomials, one with even powers and other with odd:

// A_0(x) = a_0 * x^0 + a_2 * x^1 + ... + a_n-2 * x^(n/2 - 1)
// A_1(x) = a_1 * x^0 + a_3 * x^1 + ... + a_n-1 * x^(n/2 - 1)

// this way, A(x) = A_0(x^2) + x * A_1(x^2), so if we have FFT(A_0) and FFT(A_1), we can
// "combine" them into FFT(A) in O(n) time, getting total time complexity O(n * log(n))

// all the following is just a naive implementation of the FFT algorithm

void FFT(std::vector<std::complex<double>>& a) {
	size_t n = a.size();
	size_t n_2 = n / 2;
	if (n == 1) {
		return;
	}

	// separate the odd and even terms
	std::vector<std::complex<double>> a_0(n_2);
	std::vector<std::complex<double>> a_1(n_2);
	for (size_t i = 0; i < n_2; ++i) {
		a_0[i] = a[i * 2];
		a_1[i] = a[i * 2 + 1];
	}

	FFT(a_0);
	FFT(a_1);

	std::complex<double> w(std::cos(-2 * PI / n), std::sin(-2 * PI / n)); // e^ix = cos(x) + i * sin(x)
	std::complex<double> wPowerJ(1);

	for (size_t j = 0; j < n_2; ++j) {
		a[j] = a_0[j] + wPowerJ * a_1[j];
		a[j + n_2] = a_0[j] - wPowerJ * a_1[j];

		wPowerJ *= w;
	}
}

void inverseFFT(std::vector<std::complex<double>>& a) {
	size_t n = a.size();
	size_t n_2 = n / 2;
	if (n == 1) {
		return;
	}

	// separate the odd and even terms
	std::vector<std::complex<double>> a_0(n_2);
	std::vector<std::complex<double>> a_1(n_2);
	for (size_t i = 0; i < n_2; ++i) {
		a_0[i] = a[i * 2];
		a_1[i] = a[i * 2 + 1];
	}

	inverseFFT(a_0);
	inverseFFT(a_1);

	std::complex<double> w(std::cos(2 * PI / n), std::sin(2 * PI / n)); // e^ix = cos(x) + i * sin(x)
	std::complex<double> wPowerJ(1);

	for (size_t j = 0; j < n_2; ++j) {
		a[j] = (a_0[j] + wPowerJ * a_1[j]) * 0.5;
		a[j + n_2] = (a_0[j] - wPowerJ * a_1[j]) * 0.5;

		wPowerJ *= w;
	}
}


#endif
