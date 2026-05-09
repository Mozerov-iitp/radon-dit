
//inverse_fourier_transform
#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iterator>
#include <iostream>

#ifndef M_PI
#define M_PI 3.1415926535897932045
#endif
using complex_dbl = std::complex<double>;
complex_dbl* ext_clm(complex_dbl* in, int n, int x) {
	complex_dbl* out = new complex_dbl[n];
	int i = n; while (i--) {

		out[i] = in[x + n * i];

	}
	return out;
}
void put_clm(complex_dbl* out, complex_dbl* put, int n, int x) {

	int i = n; while (i--) {

		out[i * n + x] = put[i];

	}

}

unsigned int bitReverse(unsigned int x, int log2n)

{
	int n = 0;
	int mask = 0x1;
	for (int i = 0; i < log2n; i++)

	{
		n <<= 1;
		n |= (x & 1);
		x >>= 1;
	}
	return n;
}
int get_log2n(int n) {
	int log2n = 0; while (n - 1) { n /= 2; log2n++; }

	return log2n;
}
void fft_X(float* a, complex_dbl* b, int X)
{
	int log2n = get_log2n(X);
	const complex_dbl J(0, 1);
	int n = 1 << log2n;
	for (unsigned int i = 0; i < n; ++i)
	{
		b[bitReverse(i, log2n)] = { (double)a[i],0 };
	}

	for (int s = 1; s <= log2n; ++s)

	{
		int m = 1 << s;
		int m2 = m >> 1;
		complex_dbl w(1, 0);
		complex_dbl wm = exp(-J * (M_PI / m2));
		for (int j = 0; j < m2; ++j)

		{
			for (int k = j; k < n; k += m)

			{
				complex_dbl t = w * b[k + m2];
				complex_dbl u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
			}
			w *= wm;
		}
	}
}
void fft_X(double* a, complex_dbl* b, int X)
{
	int log2n = get_log2n(X);
	const complex_dbl J(0, 1);
	int n = 1 << log2n;
	for (unsigned int i = 0; i < n; ++i)
	{
		b[bitReverse(i, log2n)] = { a[i],0 };
	}

	for (int s = 1; s <= log2n; ++s)

	{
		int m = 1 << s;
		int m2 = m >> 1;
		complex_dbl w(1, 0);
		complex_dbl wm = exp(-J * (M_PI / m2));
		for (int j = 0; j < m2; ++j)

		{
			for (int k = j; k < n; k += m)

			{
				complex_dbl t = w * b[k + m2];
				complex_dbl u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
			}
			w *= wm;
		}
	}
}
void fft_inv_X(complex_dbl* a, double* bb, int X)
{
	complex_dbl* b = new complex_dbl[X];
	int log2n = get_log2n(X);
	const complex_dbl J(0, -1);
	int n = 1 << log2n;
	for (unsigned int i = 0; i < n; ++i)
	{
		b[bitReverse(i, log2n)] = a[i];
	}

	for (int s = 1; s <= log2n; ++s)

	{
		int m = 1 << s;
		int m2 = m >> 1;
		complex_dbl w(1, 0);
		complex_dbl wm = exp(-J * (M_PI / m2));
		for (int j = 0; j < m2; ++j)

		{
			for (int k = j; k < n; k += m)

			{
				complex_dbl t = w * b[k + m2];
				complex_dbl u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
			}
			w *= wm;
		}
	}
	for (unsigned int i = 0; i < n; ++i)bb[i] = b[i].real() / n;
	delete[] b;
}
void fft_inv(complex_dbl* a, complex_dbl* b, int log2n)
{

	const complex_dbl J(0, -1);
	int n = 1 << log2n;
	for (unsigned int i = 0; i < n; ++i)
	{
		b[bitReverse(i, log2n)] = a[i];
	}

	for (int s = 1; s <= log2n; ++s)

	{
		int m = 1 << s;
		int m2 = m >> 1;
		complex_dbl w(1, 0);
		complex_dbl wm = exp(-J * (M_PI / m2));
		for (int j = 0; j < m2; ++j)

		{
			for (int k = j; k < n; k += m)

			{
				complex_dbl t = w * b[k + m2];
				complex_dbl u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
			}
			w *= wm;
		}
	}
	for (unsigned int i = 0; i < n; ++i)b[i] /= n;
}
void fft_2D_inv(complex_dbl* a, complex_dbl* b, int n)
{
	int log2n = get_log2n(n);
	int i = n; while (i--)fft_inv(&a[i * n], &b[i * n], log2n);


	{
		complex_dbl* bb = new complex_dbl[n];
		complex_dbl* bc = new complex_dbl[n * n];

		i = n; while (i--) {
			complex_dbl* pt = ext_clm(b, n, i);
			fft_inv(pt, bb, log2n);
			put_clm(bc, bb, n, i);


			delete[] pt;
		}
		memcpy(b, bc, sizeof(complex_dbl) * n * n);
		delete[] bc;

		delete[] bb;
	}
}
double* fft_2D_inv(complex_dbl* a, int n)
{
	double* ret = new double[n * n];
	complex_dbl* b = new complex_dbl[n * n];
	fft_2D_inv(a, b, n);
	double mean_out = 0;
	int p = n * n; while (p--) {//semiperiod shift 
		int x = p % n - n / 2, y = p / n - n / 2;
		int x_hT = x < 0 ? x + n : x;
		int y_hT = y < 0 ? y + n : y;
		int p_hT = x_hT + y_hT * n;
		mean_out += ret[p_hT] =  b[p].real();
		
	}
	delete[] b;

	return ret;
}

double get_mean_from_sino(float* sino, int T, int A) {
	int p = A * T; double ret = 0;
	while (p--)ret += sino[p];
	return ret / A;


}
void gauss_flt(double* in, double sigx, double sigy, int X, int Y, int w)
{
	double* out = new double[X * Y];
	double* out0 = new double[X * Y];

	// Calculate window sizes based on sigma values
	// The window covers 3 sigma in each direction, so total width = 6*sigma + 1
	int wx = (int)(3.0 * sigx) * 2 + 1;
	int wy = (int)(3.0 * sigy) * 2 + 1;

	// If user provided a window size limit (w), use the smaller of calculated or w
	if (w > 0) {
		wx = std::min(wx, w);
		wy = std::min(wy, w);
	}

	// Ensure window sizes are odd (so we have a clear center pixel)
	wx = wx | 1;  // Force to odd number using bitwise OR
	wy = wy | 1;

	// Half window sizes (distance from center to edge)
	int hwx = wx / 2;
	int hwy = wy / 2;

	// Create 2D Gaussian kernel
	double* kernel = new double[wx * wy];
	double sum = 0.0;  // For normalization

	// Fill the 2D Gaussian kernel
	for (int ky = 0; ky < wy; ky++) {
		// Vertical distance from kernel center
		int dy = ky - hwy;
		// Gaussian function value for y-direction
		double gy = exp(-(dy * dy) / (2 * sigy * sigy));

		for (int kx = 0; kx < wx; kx++) {
			// Horizontal distance from kernel center
			int dx = kx - hwx;
			// Gaussian function value for x-direction
			double gx = exp(-(dx * dx) / (2 * sigx * sigx));

			// 2D Gaussian = product of 1D Gaussians (separable property)
			kernel[ky * wx + kx] = gx * gy;
			sum += kernel[ky * wx + kx];
		}
	}

	// Normalize the kernel so that all coefficients sum to 1
	// This ensures the filtered image maintains the same overall brightness
	for (int i = 0; i < wx * wy; i++) {
		kernel[i] /= sum;
	}

	// Apply convolution to each pixel in the image
	for (int y = 0; y < Y; y++) {
		for (int x = 0; x < X; x++) {
			double pixel = 0.0;  // Accumulator for filtered pixel value

			// Convolution: sum of (neighbor pixel * kernel coefficient)
			for (int ky = 0; ky < wy; ky++) {
				// Calculate neighbor's y-coordinate
				int ny = y + ky - hwy;

				// Handle boundary conditions using mirroring
				// This prevents artifacts at image edges
				if (ny < 0) ny = -ny;                    // Mirror above top edge
				if (ny >= Y) ny = 2 * Y - ny - 1;        // Mirror below bottom edge

				for (int kx = 0; kx < wx; kx++) {
					// Calculate neighbor's x-coordinate
					int nx = x + kx - hwx;

					// Handle boundary conditions for x-direction
					if (nx < 0) nx = -nx;                    // Mirror left of left edge
					if (nx >= X) nx = 2 * X - nx - 1;        // Mirror right of right edge

					// Add weighted contribution from this neighbor
					pixel += in[ny * X + nx] * kernel[ky * wx + kx];
				}
			}

			// Store the filtered pixel value
			out[y * X + x] = pixel;
		}
	}

	// Copy filtered image back to input array
	int p = X * Y; while (p--)out0[p] = out[p];
	memcpy(in, out0, sizeof(double) * X * Y);

	// Clean up allocated memory
	delete[] out;
	delete[] out0;
	delete[] kernel;
}
void gauss_flt(float * in, double sigx, double sigy, int X, int Y, int w)
{
	double* out = new double[X * Y];
	float* out0 = new float[X * Y];

	// Calculate window sizes based on sigma values
	// The window covers 3 sigma in each direction, so total width = 6*sigma + 1
	int wx = (int)(3.0 * sigx) * 2 + 1;
	int wy = (int)(3.0 * sigy) * 2 + 1;

	// If user provided a window size limit (w), use the smaller of calculated or w
	if (w > 0) {
		wx = std::min(wx, w);
		wy = std::min(wy, w);
	}

	// Ensure window sizes are odd (so we have a clear center pixel)
	wx = wx | 1;  // Force to odd number using bitwise OR
	wy = wy | 1;

	// Half window sizes (distance from center to edge)
	int hwx = wx / 2;
	int hwy = wy / 2;

	// Create 2D Gaussian kernel
	double* kernel = new double[wx * wy];
	double sum = 0.0;  // For normalization

	// Fill the 2D Gaussian kernel
	for (int ky = 0; ky < wy; ky++) {
		// Vertical distance from kernel center
		int dy = ky - hwy;
		// Gaussian function value for y-direction
		double gy = exp(-(dy * dy) / (2 * sigy * sigy));

		for (int kx = 0; kx < wx; kx++) {
			// Horizontal distance from kernel center
			int dx = kx - hwx;
			// Gaussian function value for x-direction
			double gx = exp(-(dx * dx) / (2 * sigx * sigx));

			// 2D Gaussian = product of 1D Gaussians (separable property)
			kernel[ky * wx + kx] = gx * gy;
			sum += kernel[ky * wx + kx];
		}
	}

	// Normalize the kernel so that all coefficients sum to 1
	// This ensures the filtered image maintains the same overall brightness
	for (int i = 0; i < wx * wy; i++) {
		kernel[i] /= sum;
	}

	// Apply convolution to each pixel in the image
	for (int y = 0; y < Y; y++) {
		for (int x = 0; x < X; x++) {
			double pixel = 0.0;  // Accumulator for filtered pixel value

			// Convolution: sum of (neighbor pixel * kernel coefficient)
			for (int ky = 0; ky < wy; ky++) {
				// Calculate neighbor's y-coordinate
				int ny = y + ky - hwy;

				// Handle boundary conditions using mirroring
				// This prevents artifacts at image edges
				if (ny < 0) ny = -ny;                    // Mirror above top edge
				if (ny >= Y) ny = 2 * Y - ny - 1;        // Mirror below bottom edge

				for (int kx = 0; kx < wx; kx++) {
					// Calculate neighbor's x-coordinate
					int nx = x + kx - hwx;

					// Handle boundary conditions for x-direction
					if (nx < 0) nx = -nx;                    // Mirror left of left edge
					if (nx >= X) nx = 2 * X - nx - 1;        // Mirror right of right edge

					// Add weighted contribution from this neighbor
					pixel += in[ny * X + nx] * kernel[ky * wx + kx];
				}
			}

			// Store the filtered pixel value
			out[y * X + x] = pixel;
		}
	}

	// Copy filtered image back to input array
	int p = X * Y; while (p--)out0[p] = out[p];
	memcpy(in, out0, sizeof(float) * X * Y);

	// Clean up allocated memory
	delete[] out;
	delete[] out0;
	delete[] kernel;
}