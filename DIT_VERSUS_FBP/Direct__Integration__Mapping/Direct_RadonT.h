#pragma once
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <limits>
#include "inverse_fourier_transform.h"
#include "Generate_sino.h"
#include <numeric>

bool all_images = 0;
int method_iml; //0 DIM; 1 FBP; 2 FBP_NORM;
typedef struct {
	double psnr_img;
	double psnr_sng;
	double ssim_img;
	double var;

} TYPE_ERRORS;
TYPE_ERRORS DIM_ERR, FBP_ERR, FBP_NORM, ROT_ERR; 
typedef struct {
	float* sngm;
	unsigned char* img;
	int X;
	int A;

} TST_RADON_IMG;


#define IMAGE_IN 0
#define SINO 1
#define IMAGE_GT 2
#define IMAGE_REC 3
#define IMAGE_ERR 4
#define RESULT_TXT 5
#define IMAGE_REC_NOISE 6
#define IMAGE_ERR_NOISE 7
#define RESULT_TXT_NOISE 8
#define IMAGE_ERR_Ñ 9


char* concat(const char* s1, const char* s2)
{
	char* result = new char[(strlen(s1) + strlen(s2) + 1)];//+1 for the null-terminator

	strcpy(result, s1);
	strcat(result, s2);
	return result;
}
const char* gen_name( int func, const char* name, const char* folder, int A)
{   
	std::string str = std::to_string(A);
	const char* num_A = str.c_str();
	const char* name_png = NULL;
	const char* fld = NULL;
	const char* img_ret = NULL;
	const char* name_pfm = NULL;
	const char* fld0 = NULL;
	const char* fld_1 = NULL;
	const char* fld_fld = NULL;
	const char* img_pfm = NULL;
	const char* name_A = NULL;
	const char* name_rec = NULL;
	const char* separator = "\\";


	switch (func) {
	case IMAGE_IN:
		name_png = concat(name, ".png");
		fld = concat(folder, separator);
		img_ret = concat(fld, name_png);
		return img_ret;
		break;
	case SINO:
		name_pfm = concat(name, ".pfm");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_pfm);
		return img_ret;
		break;
	case IMAGE_GT:	
		name_png = concat(name, ".png");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case IMAGE_REC:
		name_rec = concat(name, "_rec_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".png");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case IMAGE_REC_NOISE:
		name_rec = concat(name, "_rec_noise_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".png");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case IMAGE_ERR:
		name_rec = concat(name, "_err_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".png");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case IMAGE_ERR_Ñ:
		name_rec = concat(name, "_err_ñ_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".png");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case RESULT_TXT:
		name_rec = concat(name, "_err_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".txt");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case RESULT_TXT_NOISE:
		name_rec = concat(name, "_err_noise_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".txt");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	case IMAGE_ERR_NOISE:
		name_rec = concat(name, "_err_");
		name_A = concat(name_rec, num_A);
		name_png = concat(name_A, ".png");
		fld0 = concat(folder, separator);
		fld_1 = concat(fld0, name);
		fld_fld = concat(fld_1, separator);
		img_ret = concat(fld_fld, name_png);
		return img_ret;
		break;
	default:
		return "NULL";
	}


	//std::string str = std::to_string(A);
	//const char* num_A = str.c_str();
	//char* im_out = concat(name, "_C.png");
	//char* name_ = concat(name, "_A");
	//char* name_num = concat(name_, num_A);
	//char* im_in = concat(name, ".png");
	//char* im_recon = concat(name_num, "_BA_REC.png");
	//char* im_err = concat(name, "_C_err.png");
	//char* im_pfm = concat(name, "_cir.pfm");
	
}
unsigned char* test_crcs_byte(int X, int R1, int R2, int vl1, int vl2) {
	unsigned char* out = new unsigned char[X * X];
	int p = X * X; while (p--) {
		int x = p % X - X / 2, y = p / X - X / 2;
		double r = sqrt(x * x + y * y);
		if (r > R1)out[p] = 0;
		else {
			if (r > R2)out[p] = vl1;
			else out[p] = vl2 + vl1;


		}



	}
	return out;

}
int log2_int( int n) {
	if (n == 0) return -1;

	int result = -1;
	while (n) {
		result++;
		n >>= 1;
	}
	return result;
}
int power_int(int base, int exponent) {
	int result = 1;
	for (int i = 0; i < exponent; i++) {
		result *= base;
	}
	return result;
}
unsigned char* test_file_byte(unsigned char* in, int X_in, int Y,  int * X) {
	
	int N = *X = power_int(2, log2_int(X_in)); 
	
	unsigned char* out = new unsigned char[N * N];
	int p = N*N; while (p--) {
		int x = p % N - N/2, y = p /N - N/2;
		int xp = x + N / 2;
		int yp = y + N / 2;
		int yp0 = yp < Y ? yp : Y - 1;
		int p0 = xp + yp0 * X_in;
		double r = sqrt(x * x + y * y);
		bool do_it = r < (N) / 2 - 1;
		
		if (do_it)out[p] = in[p0];
		else out[p] = 0;

	}
	return out;

}


void test_file_byte_sino(TST_RADON_IMG* circle, int A,unsigned char* in, int X_in, int Y, int* X)
{
	// generate sinogramm and byte image (as a ground truth) 
	// for the next DIM process
	circle->img = test_file_byte(in, X_in, Y, X);
	circle->sngm =discreteRadonTransformDoublePtr(
		prepareForRadon(circle->img, *X),
		*X, A, *X);
		//test_sino(*X, A, circle->img);
	circle->A = A; 
	circle->X = *X;
	

}
double* forward_result( int X,  int A, double* in)
{
return discreteRadonTransformDoublePtr(prepareForRadon(in, X),X, A, X);


}



/**
 * Calculates Peak Signal-to-Noise Ratio (PSNR) between two images.
 *
 * @param result    Pointer to the result image pixel array (8-bit unsigned)
 * @param gr_truth  Pointer to the ground truth image pixel array (8-bit unsigned)
 * @param size      Number of pixels to compare (width * height * channels)
 * @return          PSNR value in decibels (dB), or -1.0 in case of error
 *
 * @note Expects 8-bit images (pixel values 0-255)
 * @note Returns INFINITY for identical images (MSE = 0)
 */
double get_PSNR(double* result, unsigned char* gr_truth, int size) {
	// 1. Validate input parameters
	if (result == nullptr || gr_truth == nullptr) {
		std::cerr << "Error: null pointer to image data" << std::endl;
		return -1.0;
	}

	if (size <= 0) {
		std::cerr << "Error: invalid array size" << std::endl;
		return -1.0;
	}

	// 2. Calculate Mean Squared Error (MSE)
	double mse = 0.0;
	long long diff; // Use long long to prevent overflow
	int X = sqrt(size); int cnt = 0;
	for (int i = 0; i < size; i++) {
		int x = i % X - X / 2, y = i / X - X / 2;
		bool do_it = sqrt((double)x * x + y * y) < X / 2;
		if (do_it) {
			diff =
				static_cast<int>(result[i]) - static_cast<int>(gr_truth[i]);
			cnt++;
		}
		else diff = 0;
		mse += static_cast<double>(diff * diff);
	}
	mse /= static_cast<double>(cnt);

	// 3. Handle the case of identical images
	if (mse == 0.0) {
		return INFINITY;
	}

	// 4. Calculate PSNR
	const double MAX_PIXEL_VALUE = 255.0;
	double psnr = 10.0 * log10((MAX_PIXEL_VALUE * MAX_PIXEL_VALUE) / mse);

	return psnr;
}

double get_PSNR(unsigned char* result, unsigned char* gr_truth, int size) {
	// 1. Validate input parameters
	if (result == nullptr || gr_truth == nullptr) {
		std::cerr << "Error: null pointer to image data" << std::endl;
		return -1.0;
	}

	if (size <= 0) {
		std::cerr << "Error: invalid array size" << std::endl;
		return -1.0;
	}

	// 2. Calculate Mean Squared Error (MSE)
	double mse = 0.0;
	long long diff; // Use long long to prevent overflow
	int X = sqrt(size);
	for (int i = 0; i < size; i++) {
		int x = i % X - X / 2, y = i / X - X / 2;
		bool do_it = sqrt((double)x * x + y * y) < X / 2;
		diff = static_cast<int>(result[i]) - static_cast<int>(gr_truth[i]);
		if(do_it)mse += static_cast<double>(diff * diff);
	}
	mse /= static_cast<double>(size);

	// 3. Handle the case of identical images
	if (mse == 0.0) {
		return INFINITY;
	}

	// 4. Calculate PSNR
	const double MAX_PIXEL_VALUE = 255.0;
	double psnr = 10.0 * log10((MAX_PIXEL_VALUE * MAX_PIXEL_VALUE) / mse);

	return psnr;
}
double get_PSNR_(double* result, unsigned char* gr_truth, int size) {
	// 1. Validate input parameters
	int p = size; double max = 0;
	while (p--) {
		if (gr_truth[p] > max)max = gr_truth[p];
	}
	if (result == nullptr || gr_truth == nullptr) {
		std::cerr << "Error: null pointer to image data" << std::endl;
		return -1.0;
	}

	if (size <= 0) {
		std::cerr << "Error: invalid array size" << std::endl;
		return -1.0;
	}

	// 2. Calculate Mean Squared Error (MSE)
	double mse = 0.0;
	double diff; // Use long long to prevent overflow
	int cnt = 0;
	for (int i = 0; i < size; i++) {



		{
			diff =
				(double)(result[i]) - (double)(gr_truth[i]);
			cnt++;
		}
		mse += static_cast<double>(diff * diff);
	}
	mse /= static_cast<double>(cnt);

	// 3. Handle the case of identical images
	if (mse == 0.0) {
		return INFINITY;
	}

	// 4. Calculate PSNR
	const double MAX_PIXEL_VALUE = max;
	double psnr = 10.0 * log10((MAX_PIXEL_VALUE * MAX_PIXEL_VALUE) / mse);

	return psnr;
}

double get_PSNR(double* result, float* gr_truth, int size) {
	// 1. Validate input parameters
	int p = size; double max = 0; 
	while (p--) {
		if (gr_truth[p] > max)max = gr_truth[p];
	}
	if (result == nullptr || gr_truth == nullptr) {
		std::cerr << "Error: null pointer to image data" << std::endl;
		return -1.0;
	}

	if (size <= 0) {
		std::cerr << "Error: invalid array size" << std::endl;
		return -1.0;
	}

	// 2. Calculate Mean Squared Error (MSE)
	double mse = 0.0;
double diff; // Use long long to prevent overflow
 int cnt = 0; 
	for (int i = 0; i < size; i++) {



	    {
			diff =
				(double)(result[i]) - (double)(gr_truth[i]);
			cnt++;
		}
		mse += static_cast<double>(diff * diff);
	}
	mse /= static_cast<double>(cnt);

	// 3. Handle the case of identical images
	if (mse == 0.0) {
		return INFINITY;
	}

	// 4. Calculate PSNR
	const double MAX_PIXEL_VALUE = max;
	double psnr = 10.0 * log10((MAX_PIXEL_VALUE * MAX_PIXEL_VALUE) / mse);

	return psnr;
}
void get_mu_sig(unsigned char* results, int size, double*mu, double* sig) {
	
	int p = size; double mn = 0;
	while (p--) {
		mn += results[p];
	}
	*mu = mn /= size;
	p = size; double sg = 0;
	while (p--) {
		double vl = mn - results[p];
		sg += vl * vl;
	}
	*sig = sqrt(sg / size);

}

void PSNR(unsigned char* result, unsigned char* gr_truth, int size) {

	double psnr_value = get_PSNR(result, gr_truth, size);
	double mu, sig, mu_r, sig_r; 
	get_mu_sig(gr_truth, size, &mu, &sig);
	get_mu_sig(result, size, &mu_r, &sig_r);
	std::cout << "ground truth mean = " << mu << 
		" std deviation = " << sig << std::endl;
	std::cout << "reconstruction result mean = " << mu_r <<
		" std deviation = " << sig_r << std::endl;
	std::cout << " result mean /(gt mean) = " << mu_r / mu <<
		" std deviation = " << sig_r / sig << std::endl;

	if (psnr_value != -1.0) {
		if (std::isinf(psnr_value)) {
			std::cout << "Images are identical (PSNR = infinity)" << std::endl;
		}
		else {
			std::cout << "PSNR = " << psnr_value << " dB" << std::endl;
		}
	}
}
void get_mean_sigma(unsigned char* in, int total_pixels, double* mean, double* sigma) {

	if (total_pixels == 0) {
		*mean = 0.0;
		*sigma = 0.0;
		return;
	}

	// First pass: calculate mean
	double sum = 0.0;
	for (int i = 0; i < total_pixels; ++i) {
		sum += in[i];
	}
	*mean = sum / total_pixels;

	// Second pass: calculate variance
	double sum_sq_diff = 0.0;
	for (int i = 0; i < total_pixels; ++i) {
		double diff = in[i] - *mean;
		sum_sq_diff += diff * diff;
	}

	// Calculate standard deviation
	*sigma = sqrt(sum_sq_diff / total_pixels);
}
void get_mean_sigma(double* in, int total_pixels, double* mean, double* sigma) {

	if (total_pixels == 0) {
		*mean = 0.0;
		*sigma = 0.0;
		return;
	}

	// First pass: calculate mean
	double sum = 0.0;
	for (int i = 0; i < total_pixels; ++i) {
		sum += in[i];
	}
	*mean = sum / total_pixels;

	// Second pass: calculate variance
	double sum_sq_diff = 0.0;
	for (int i = 0; i < total_pixels; ++i) {
		double diff = in[i] - *mean;
		sum_sq_diff += diff * diff;
	}

	// Calculate standard deviation
	*sigma = sqrt(sum_sq_diff / total_pixels);
}
void PSNR(double* result, float* gr_truth, int size) {

	double psnr_value = get_PSNR(result, gr_truth, size);

	if (psnr_value != -1.0) {
		if (std::isinf(psnr_value)) {
			std::cout << "Images are identical (PSNR = infinity)" << std::endl;
		}
		else {
			std::cout << "PSNR of sinogram = " << psnr_value << " dB" << std::endl;
		}
	}
}
double MAE(unsigned char* recover, unsigned char* gr_truth, int size) {
	double mae = 0;
	int X = sqrt(size);
	for (int i = 0; i < size; i++) {
		int x = i % X - X / 2, y = i / X - X / 2;
		bool do_it = sqrt((double)x * x + y * y) < X / 2;
		if(do_it)mae += abs(recover[i] - gr_truth[i]);
	}
	return mae / size / 255 * 100;
}
double MAE(double* recover, unsigned char* gr_truth, int size) {
	double mae = 0;
	int X = sqrt(size); int cnt = 0; 
	for (int i = 0; i < size; i++) {
		int x = i % X - X / 2, y = i / X - X / 2;
		bool do_it = sqrt((double)x * x + y * y) < X / 2;
		if (do_it) { mae += abs(recover[i] - gr_truth[i]); cnt++; }
	}
	return mae / cnt / 255 * 100;
}
double SSIM(double * recover, unsigned char* gr_truth, int size) {
	// Constants for SSIM calculation (as per the original paper)
	const double C1 = 6.5025;  // (0.01 * 255)^2
	const double C2 = 58.5225; // (0.03 * 255)^2
	const double C3 = 29.26125; // C2 / 2

	// Initialize variables for statistics
	double sum_x = 0.0, sum_y = 0.0;
	double sum_x2 = 0.0, sum_y2 = 0.0, sum_xy = 0.0;

	// Calculate mean values and prepare for variance/covariance
	for (int i = 0; i < size; i++) {
		double x = static_cast<double>(recover[i]);
		double y = static_cast<double>(gr_truth[i]);

		sum_x += x;
		sum_y += y;
		sum_x2 += x * x;
		sum_y2 += y * y;
		sum_xy += x * y;
	}

	// Calculate means
	double mean_x = sum_x / size;
	double mean_y = sum_y / size;

	// Calculate variances and covariance
	double var_x = (sum_x2 / size) - (mean_x * mean_x);
	double var_y = (sum_y2 / size) - (mean_y * mean_y);
	double cov_xy = (sum_xy / size) - (mean_x * mean_y);

	// Ensure variances are non-negative (numerical stability)
	var_x = std::max(var_x, 0.0);
	var_y = std::max(var_y, 0.0);

	// Calculate SSIM components
	double luminance = (2 * mean_x * mean_y + C1) / (mean_x * mean_x + mean_y * mean_y + C1);
	double contrast = (2 * std::sqrt(var_x) * std::sqrt(var_y) + C2) / (var_x + var_y + C2);
	double structure = (cov_xy + C3) / (std::sqrt(var_x) * std::sqrt(var_y) + C3);

	// Combined SSIM value
	double ssim = luminance * contrast * structure;

	return ssim;
}

double SSIM(unsigned char* recover, unsigned char* gr_truth, int size) {
	// Constants for SSIM calculation (as per the original paper)
	const double C1 = 6.5025;  // (0.01 * 255)^2
	const double C2 = 58.5225; // (0.03 * 255)^2
	const double C3 = 29.26125; // C2 / 2

	// Initialize variables for statistics
	double sum_x = 0.0, sum_y = 0.0;
	double sum_x2 = 0.0, sum_y2 = 0.0, sum_xy = 0.0;

	// Calculate mean values and prepare for variance/covariance
	for (int i = 0; i < size; i++) {
		double x = static_cast<double>(recover[i]);
		double y = static_cast<double>(gr_truth[i]);

		sum_x += x;
		sum_y += y;
		sum_x2 += x * x;
		sum_y2 += y * y;
		sum_xy += x * y;
	}

	// Calculate means
	double mean_x = sum_x / size;
	double mean_y = sum_y / size;

	// Calculate variances and covariance
	double var_x = (sum_x2 / size) - (mean_x * mean_x);
	double var_y = (sum_y2 / size) - (mean_y * mean_y);
	double cov_xy = (sum_xy / size) - (mean_x * mean_y);

	// Ensure variances are non-negative (numerical stability)
	var_x = std::max(var_x, 0.0);
	var_y = std::max(var_y, 0.0);

	// Calculate SSIM components
	double luminance = (2 * mean_x * mean_y + C1) / (mean_x * mean_x + mean_y * mean_y + C1);
	double contrast = (2 * std::sqrt(var_x) * std::sqrt(var_y) + C2) / (var_x + var_y + C2);
	double structure = (cov_xy + C3) / (std::sqrt(var_x) * std::sqrt(var_y) + C3);

	// Combined SSIM value
	double ssim = luminance * contrast * structure;

	return ssim;
}
double SSIM_windowed(unsigned char* recover, unsigned char* gr_truth, int width, int height) {
	// SSIM constants
	const double K1 = 0.01;
	const double K2 = 0.03;
	const double L = 255.0; // dynamic range for 8-bit images
	const double C1 = (K1 * L) * (K1 * L);
	const double C2 = (K2 * L) * (K2 * L);

	const int window_size = 2;
	const int stride = 2; // Non-overlapping windows as in original implementation

	double total_ssim = 0.0;
	int window_count = 0;

	// Slide window over the image
	for (int y = 0; y <= height - window_size; y += stride) {
		for (int x = 0; x <= width - window_size; x += stride) {
			double sum_x = 0.0, sum_y = 0.0;
			double sum_x2 = 0.0, sum_y2 = 0.0, sum_xy = 0.0;

			// Calculate statistics for current window
			for (int j = 0; j < window_size; j++) {
				for (int i = 0; i < window_size; i++) {
					int idx = (y + j) * width + (x + i);
					double pixel_x = static_cast<double>(recover[idx]);
					double pixel_y = static_cast<double>(gr_truth[idx]);

					sum_x += pixel_x;
					sum_y += pixel_y;
					sum_x2 += pixel_x * pixel_x;
					sum_y2 += pixel_y * pixel_y;
					sum_xy += pixel_x * pixel_y;
				}
			}

			int window_pixels = window_size * window_size;
			double mean_x = sum_x / window_pixels;
			double mean_y = sum_y / window_pixels;

			double var_x = (sum_x2 / window_pixels) - (mean_x * mean_x);
			double var_y = (sum_y2 / window_pixels) - (mean_y * mean_y);
			double cov_xy = (sum_xy / window_pixels) - (mean_x * mean_y);

			// Calculate SSIM for this window
			double numerator = (2 * mean_x * mean_y + C1) * (2 * cov_xy + C2);
			double denominator = (mean_x * mean_x + mean_y * mean_y + C1) * (var_x + var_y + C2);

			double window_ssim = numerator / denominator;
			total_ssim += window_ssim;
			window_count++;
		}
	}

	// Return average SSIM over all windows
	return (window_count > 0) ? total_ssim / window_count : 0.0;
}
void ERRORS(std::ofstream& file,
	unsigned char* result, unsigned char* gr_truth, 
	double* result_f, float* gr_truth_f,
	int size, int size_f) {
	double mean, sigma, mean_gt, sigma_gt;
	get_mean_sigma(result, size, &mean, &sigma);
	get_mean_sigma(gr_truth, size, &mean_gt, &sigma_gt);
	double mae = MAE(result, gr_truth, size);
	double ssim = SSIM(result, gr_truth, size);
	double psnr_value = get_PSNR(result, gr_truth, size);
	double psnr_value_f = get_PSNR(result_f, gr_truth_f, size_f);
	std::cout << "mean GT = " << mean_gt
		<< " - mean result = " << mean  
		<< " - mean / mean GT = " << mean/ mean_gt <<
		std::endl;
	file << "mean GT = " << mean_gt
		<< " - mean result = " << mean 
		<< " - mean / mean GT = " << mean / mean_gt <<
		std::endl;
	std::cout << "sigma GT = " << sigma_gt
		<< " - sigma result = " << sigma 
		<< " - sigma / sigma GT = " << sigma / sigma_gt <<
		std::endl;
	
	file << "sigma GT = " << sigma_gt
		<< " - sigma result = " << sigma  
		<< " - sigma / sigma GT = " << sigma / sigma_gt <<
		std::endl;
	if (psnr_value != -1.0) {
		if (std::isinf(psnr_value)) {
			std::cout << "Images are identical (PSNR = infinity)" << std::endl;
			file << "Images are identical (PSNR = infinity)" << std::endl;
		}
		else {
			std::cout << "PSNR of reconstruction = " << psnr_value << " dB" << std::endl;
			file << "PSNR of reconstruction = " << psnr_value << " dB" << std::endl;
		}
	}
	
	std::cout << "MAE of reconstruction = " << mae << " \%" << std::endl;
	file << "MAE of reconstruction = " << mae << " \%" << std::endl;
	std::cout << "SSIM of reconstruction = " << ssim << " \%" << std::endl;
	file << "SSIM of reconstruction = " << ssim << " \%" << std::endl;
	if (psnr_value_f != -1.0) {
		if (std::isinf(psnr_value_f)) {
			std::cout << "Images are identical (PSNR = infinity)" << std::endl;
			file << "Images are identical (PSNR = infinity)" << std::endl;
		}
		else {
			std::cout << "PSNR of sinogram = " << psnr_value_f << " dB" << std::endl;
			file << "PSNR of sinogram = " << psnr_value_f << " dB" << std::endl;
		}
	}

}
void ERRORS(std::ofstream& file,
	double * result, unsigned char* gr_truth,
	double* result_f, float* gr_truth_f,
	int size, int size_f) {
	double mean, sigma, mean_gt, sigma_gt;
	get_mean_sigma(result, size, &mean, &sigma);
	get_mean_sigma(gr_truth, size, &mean_gt, &sigma_gt);
	double mae = MAE(result, gr_truth, size);
	double ssim = SSIM(result, gr_truth, size);
	double psnr_value = get_PSNR(result, gr_truth, size);
	double psnr_value_f = get_PSNR(result_f, gr_truth_f, size_f);

	std::cout 
		<< " - sigma / sigma GT = " << sigma / sigma_gt <<
		std::endl;

	file 
		<< " - sigma / sigma GT = " << sigma / sigma_gt <<
		std::endl;
	if (psnr_value != -1.0) {
		if (std::isinf(psnr_value)) {
			std::cout << "Images are identical (PSNR = infinity)" << std::endl;
			file << "Images are identical (PSNR = infinity)" << std::endl;
		}
		else {
			std::cout << "PSNR of reconstruction = " << psnr_value << " dB" << std::endl;
			file << "PSNR of reconstruction = " << psnr_value << " dB" << std::endl;
		}
	}

	std::cout << "SSIM of reconstruction = " << ssim << " \%" << std::endl;
	file << "SSIM of reconstruction = " << ssim << " \%" << std::endl;
	if (psnr_value_f != -1.0) {
		if (std::isinf(psnr_value_f)) {
			std::cout << "Images are identical (PSNR = infinity)" << std::endl;
			file << "Images are identical (PSNR = infinity)" << std::endl;
		}
		else {
			std::cout << "PSNR of sinogram = " << psnr_value_f << " dB" << std::endl;
			file << "PSNR of sinogram = " << psnr_value_f << " dB" << std::endl;
		}
	}
	if (all_images) {
		if (method_iml == 0) {
			DIM_ERR.psnr_img += psnr_value;
			DIM_ERR.psnr_sng += psnr_value_f;
			DIM_ERR.ssim_img += ssim;
			DIM_ERR.var += sigma / sigma_gt;

		}
		if (method_iml == 1) {
			FBP_ERR.psnr_img += psnr_value;
			FBP_ERR.psnr_sng += psnr_value_f;
			FBP_ERR.ssim_img += ssim;
			FBP_ERR.var += sigma / sigma_gt;

		}
		if (method_iml == 2) {
			FBP_NORM.psnr_img += psnr_value;
			FBP_NORM.psnr_sng += psnr_value_f;
			FBP_NORM.ssim_img += ssim;
			FBP_NORM.var += sigma / sigma_gt;

		}
	}
}//DIM_ERR, FBP_ERR, FBP_NORM, ROT_ERR; 
void MEAN_ERRORS() {
	std::cout << std::endl;
	if(DIM_ERR.psnr_img){
		std::cout << "=== Direct Integration Method mean Errors ===" <<
			std::endl;
		std::cout << std::endl;
		std::cout << " PSNR of reconstructed image  = " << DIM_ERR.psnr_img / 10 <<
			std::endl;
		std::cout << " SSIM of reconstructed image  = " << DIM_ERR.ssim_img / 10 <<
			std::endl;
		std::cout << " sigma / sigma GT  = " << DIM_ERR.var / 10 <<
			std::endl;
		std::cout << " PSNR  of inverse   sinogram  = " << DIM_ERR.psnr_sng / 10 <<
			std::endl;
		std::cout << "=============================================" <<
			std::endl;
		std::cout << std::endl;
	}
	if (FBP_ERR.psnr_img) {
		std::cout << "===             FBP  mean Errors          ===" <<
			std::endl;
		std::cout << " PSNR of reconstructed image  = " << FBP_ERR.psnr_img / 10 <<
			std::endl;
		std::cout << " SSIM of reconstructed image  = " << FBP_ERR.ssim_img / 10 <<
			std::endl;
		std::cout << " sigma / sigma GT  = " << FBP_ERR.var / 10 <<
			std::endl;
		std::cout << " PSNR  of inverse   sinogram  = " << FBP_ERR.psnr_sng / 10 <<
			std::endl;
		std::cout << "=============================================" <<
			std::endl;
		std::cout << std::endl;
	}
	if (FBP_NORM.psnr_img) {
		std::cout << "=== FBP normalized            mean Errors ===" <<
			std::endl;
		std::cout << " PSNR of reconstructed image  = " << FBP_NORM.psnr_img / 10 <<
			std::endl;
		std::cout << " SSIM of reconstructed image  = " << FBP_NORM.ssim_img / 10 <<
			std::endl;
		std::cout << " sigma / sigma GT  = " << FBP_NORM.var / 10 <<
			std::endl;
		std::cout << " PSNR  of inverse   sinogram  = " << FBP_NORM.psnr_sng / 10 <<
			std::endl;
		std::cout << "=============================================" <<
			std::endl;
		std::cout << std::endl;
	}
	if (ROT_ERR.psnr_img) {
		std::cout << "===     Double Rotation mean Errors       ===" <<
			std::endl;
		std::cout << " PSNR of reconstructed image  = " << ROT_ERR.psnr_img / 10 <<
			std::endl;
		std::cout << " PSNR  of inverse   sinogram  = " << ROT_ERR.psnr_sng / 10 <<
			std::endl;
	}

	

	
}//DIM_ERR, FBP_ERR, FBP_NORM, ROT_ERR; 
void Norm_with_GT(
	double* result, unsigned char* gr_truth, int size) {
	double mean, sigma, mean_gt, sigma_gt;
	get_mean_sigma(result, size, &mean, &sigma);
	get_mean_sigma(gr_truth, size, &mean_gt, &sigma_gt);
	double* out = new double[size];
	int p = size; while (p--)out[p] =
		(result[p] - mean) / sigma * sigma_gt + mean_gt;
	memcpy(result, out, sizeof(double) * size);
	delete[]out;
}
void truncate_double_to_byte(double* in,
	unsigned char* result,  int T) {
	double min = in[T * T - 1], max = min;
	
	int p = T * T; while (p--) {
		if (min > in[p])min = in[p];
		if (max < in[p]) max = in[p];

	}
	double sc = max - min;

	    p = T * T; while (p--) {
		int x = p % T - T / 2, y = p / T - T / 2;
		 {  
			double vl = (in[p]); if (vl < 0)vl = 0;

			
			result[p] = vl/max*255;
			
		}
		

	}

}
void truncate_double_to_byte_ñ(double* in,
	unsigned char* result, unsigned char* gt, unsigned char* err_m, int T) {
	
	int p = T * T; while (p--) {
		int x = p % T - T / 2, y = p / T - T / 2;
		if (sqrt((double)x * x + y * y) < T / 2) {
			int vl = in[p] + 0.5;
			if (vl > 255)vl = 255; if (vl < 0)vl = 0;
			result[p] = vl;
			int err = (gt[p] - vl)*10;
			int errp = err >= 0 ? err : 0; if (errp > 255)errp = 255;
			int errm = err <  0 ? abs(err) : 0; if (errm > 255)errm = 255;
			
			err_m[p] = errm;
			err_m[p + T * T] = errp > errm ? errp : errm;
			err_m[p + 2*T * T] = errp;
		}
		else { err_m[p] = 0; result[p] = 0; }

	}

}

void truncate_double_to_byte(double* in,
	unsigned char* result, unsigned char* gt, unsigned char* err_m, int T) {
	int p = T * T; while (p--) {
		int x = p % T - T / 2, y = p / T - T / 2;
		if(sqrt((double)x*x+y*y)<T/2){
			int vl = in[p] + 0.5;
			if (vl > 255)vl = 255; if (vl < 0)vl = 0;
			result[p] = vl;
			int err = abs(gt[p] - vl) * 10; if (err > 255)err = 255;
			err_m[p] = err;
		}
		else { err_m[p] = 0; result[p] = 0;  }
	
	}

}
double generate_gaussian_noise(double sigma) {
	static std::default_random_engine generator;
	static std::normal_distribution<double> distribution(0.0, 1.0);
	return distribution(generator) * sigma;
}
void add_shift_noise(float* image, int X, int Y, int Sh,  double sh_pr) {

	float* out = new float[X * Y];
	int  sh_r = sh_pr * X;
	
		for (int y = 0; y < Y; y++) {
			int sh = (sh_pr * X < 2) ? 0 : rand() % sh_r - sh_r / 2;
			
			
			for (int x = 0; x < X; x++){
				int xp = x + sh + Sh;
				if (xp < 0)xp = 0;
				if (xp > X - 1)xp = X - 1;
				out[xp + y * X] = image[x + y * X];
			}
				
		}
	memcpy(image, out, sizeof(float) * X * Y);
	delete[] out;
}
void add_poisson_noise(float* image, int size, double photon_count_I0) {
	if (!image || size <= 0 || photon_count_I0 <= 0) return;

	static std::mt19937 gen(std::random_device{}());

	for (int i = 0; i < size; ++i) {
		// 1. Convert sinogram value (attenuation) back to photon count
		// image[i] = -ln(I / I0) => I = I0 * exp(-image[i])
		double expected_photons = photon_count_I0 * std::exp(-static_cast<double>(image[i]));

		// 2. Apply Poisson distribution based on the expected number of photons
		std::poisson_distribution<int> dist(expected_photons);
		double noisy_photons = static_cast<double>(dist(gen));

		// 3. Convert back to attenuation space (sinogram)
		// New value = -ln(noisy_I / I0)
		// Clamp to 1 photon minimum to avoid ln(0)
		if (noisy_photons < 1.0) noisy_photons = 1.0;

		image[i] = static_cast<float>(-std::log(noisy_photons / photon_count_I0));

		// 4. Ensure no negative attenuation (physically impossible)
		if (image[i] < 0.0f) image[i] = 0.0f;
	}
}

void add_gaussian_noise(float* image, int size, double noise_level_percent) {
	if (!image || size <= 0 || noise_level_percent <= 0.0) return;

	// 1. Calculate the mean value to use as a baseline for noise scaling
	// This is more robust than using max_value as it ignores isolated outliers.
	double sum = std::accumulate(image, image + size, 0.0);
	double mean_value = sum / size;

	// 2. Define the standard deviation (sigma) relative to the mean
	double sigma = mean_value * noise_level_percent;

	// 3. Initialize the random engine and distribution
	// 'static' ensures the generator is initialized only once for better performance.
	static std::mt19937 gen(std::random_device{}());
	std::normal_distribution<float> dist(0.0f, static_cast<float>(sigma));

	for (int i = 0; i < size; ++i) {
		float noisy_value = image[i] + dist(gen);

		// 4. Clamping:
		// We only clamp to 0 because negative values are physically impossible 
		// in a sinogram (line integrals of attenuation). 
		// We do NOT clamp the upper bound to preserve the Gaussian distribution's tail.
		image[i] = (noisy_value < 0.0f) ? 0.0f : noisy_value;
	}
}double bi_cube(unsigned char* in, int x, int y, int X, double angle) {
	// Apply inverse rotation to get coordinates in original image
	// Rotation matrix for -angle:
	// [ cos(-a)  -sin(-a) ] = [ cos(a)   sin(a) ]
	// [ sin(-a)   cos(-a) ]   [ -sin(a)  cos(a) ]
	double cos_a = cos(angle);  // cos(-angle) = cos(angle)
	double sin_a = sin(angle);  // sin(-angle) = -sin(angle)

	// Inverse rotation: (x_orig, y_orig) = R^(-1) * (x, y)
	// where R^(-1) is rotation by -angle
	double x_rot = x * cos_a + y * sin_a;    // x*cos(angle) + y*sin(angle)
	double y_rot = -x * sin_a + y * cos_a;   // -x*sin(angle) + y*cos(angle)

	// Convert from center-based to corner-based coordinates
	// Origin was at (X/2, X/2), now move to (0, 0)
	double x_orig = x_rot + X / 2.0;
	double y_orig = y_rot + X / 2.0;

	// If coordinates are outside image boundaries, return 0
	if (x_orig < 0 || x_orig >= X - 1 || y_orig < 0 || y_orig >= X - 1) {
		return 0.0;
	}

	// Get integer and fractional parts for interpolation
	int x_int = (int)floor(x_orig);
	int y_int = (int)floor(y_orig);
	double dx = x_orig - x_int;
	double dy = y_orig - y_int;

	/**
	 * Cubic interpolation kernel (Catmull-Rom spline)
	 * This kernel provides good quality interpolation with minimal artifacts
	 *
	 * @param t: distance from interpolation point (-2 <= t <= 2)
	 * @return weight coefficient
	 */
	auto cubic_kernel = [](double t) -> double {
		double abs_t = fabs(t);

		if (abs_t <= 1.0) {
			// Close to interpolation point: cubic polynomial
			return 1.0 + abs_t * abs_t * (1.5 * abs_t - 2.5);
		}
		else if (abs_t <= 2.0) {
			// Further from interpolation point: different cubic
			return 2.0 + abs_t * (-4.0 + abs_t * (2.5 - 0.5 * abs_t));
		}
		return 0.0;  // Outside [-2, 2] range
		};

	// Perform bi-cubic interpolation over 4x4 neighborhood
	double result = 0.0;

	// Iterate over 4x4 neighborhood around the interpolation point
	for (int m = -1; m <= 2; m++) {
		// Calculate weight for x-direction
		double weight_x = cubic_kernel(dx - m);

		for (int n = -1; n <= 2; n++) {
			// Calculate weight for y-direction
			double weight_y = cubic_kernel(dy - n);

			// Get coordinates of sample pixel
			int sample_x = x_int + m;
			int sample_y = y_int + n;

			// Clamp coordinates to image boundaries (edge handling)
			// This uses "clamp to edge" strategy
			if (sample_x < 0) sample_x = 0;
			if (sample_x >= X) sample_x = X - 1;
			if (sample_y < 0) sample_y = 0;
			if (sample_y >= X) sample_y = X - 1;

			// Get pixel value from input image
			// Image is stored in row-major order: row * width + column
			unsigned char pixel_value = in[sample_y * X + sample_x];

			// Add weighted contribution to result
			result += pixel_value * weight_x * weight_y;
		}
	}

	// Ensure result is in valid range [0, 255]
	// Note: Due to interpolation, values might slightly exceed bounds
	if (result < 0.0) result = 0.0;
	if (result > 255.0) result = 255.0;

	return result;
}

double* rotate(unsigned char* in, int X, double angle) {
	double* out = new double[X * X];
	int p = X * X; while (p--) {
		int x = p % X - X / 2, y = p / X - X / 2;
		bool doit = sqrt((double)x * x + y * y) < X / 2;
		if (doit) {
			out[p] = bi_cube(in, x, y, X, angle);

		}
		else out[p] = 0;
	}
	return out;
}
double bi_cube(double* in, int x, int y, int X, double angle) {
	// Apply inverse rotation to get coordinates in original image
	// Rotation matrix for -angle:
	// [ cos(-a)  -sin(-a) ] = [ cos(a)   sin(a) ]
	// [ sin(-a)   cos(-a) ]   [ -sin(a)  cos(a) ]
	double cos_a = cos(angle);  // cos(-angle) = cos(angle)
	double sin_a = sin(angle);  // sin(-angle) = -sin(angle)

	// Inverse rotation: (x_orig, y_orig) = R^(-1) * (x, y)
	// where R^(-1) is rotation by -angle
	double x_rot = x * cos_a + y * sin_a;    // x*cos(angle) + y*sin(angle)
	double y_rot = -x * sin_a + y * cos_a;   // -x*sin(angle) + y*cos(angle)

	// Convert from center-based to corner-based coordinates
	// Origin was at (X/2, X/2), now move to (0, 0)
	double x_orig = x_rot + X / 2.0;
	double y_orig = y_rot + X / 2.0;

	// If coordinates are outside image boundaries, return 0
	if (x_orig < 0 || x_orig >= X - 1 || y_orig < 0 || y_orig >= X - 1) {
		return 0.0;
	}

	// Get integer and fractional parts for interpolation
	int x_int = (int)floor(x_orig);
	int y_int = (int)floor(y_orig);
	double dx = x_orig - x_int;
	double dy = y_orig - y_int;

	/**
	 * Cubic interpolation kernel (Catmull-Rom spline)
	 * This kernel provides good quality interpolation with minimal artifacts
	 *
	 * @param t: distance from interpolation point (-2 <= t <= 2)
	 * @return weight coefficient
	 */
	auto cubic_kernel = [](double t) -> double {
		double abs_t = fabs(t);

		if (abs_t <= 1.0) {
			// Close to interpolation point: cubic polynomial
			return 1.0 + abs_t * abs_t * (1.5 * abs_t - 2.5);
		}
		else if (abs_t <= 2.0) {
			// Further from interpolation point: different cubic
			return 2.0 + abs_t * (-4.0 + abs_t * (2.5 - 0.5 * abs_t));
		}
		return 0.0;  // Outside [-2, 2] range
		};

	// Perform bi-cubic interpolation over 4x4 neighborhood
	double result = 0.0;

	// Iterate over 4x4 neighborhood around the interpolation point
	for (int m = -1; m <= 2; m++) {
		// Calculate weight for x-direction
		double weight_x = cubic_kernel(dx - m);

		for (int n = -1; n <= 2; n++) {
			// Calculate weight for y-direction
			double weight_y = cubic_kernel(dy - n);

			// Get coordinates of sample pixel
			int sample_x = x_int + m;
			int sample_y = y_int + n;

			// Clamp coordinates to image boundaries (edge handling)
			// This uses "clamp to edge" strategy
			if (sample_x < 0) sample_x = 0;
			if (sample_x >= X) sample_x = X - 1;
			if (sample_y < 0) sample_y = 0;
			if (sample_y >= X) sample_y = X - 1;

			// Get pixel value from input image
			// Image is stored in row-major order: row * width + column
			unsigned char pixel_value = in[sample_y * X + sample_x];

			// Add weighted contribution to result
			result += pixel_value * weight_x * weight_y;
		}
	}

	// Ensure result is in valid range [0, 255]
	// Note: Due to interpolation, values might slightly exceed bounds
	if (result < 0.0) result = 0.0;
	if (result > 255.0) result = 255.0;

	return result;
}

double* rotate(double * in, int X, double angle) {
	double* out = new double[X * X];
	int p = X * X; while (p--) {
		int x = p % X - X / 2, y = p / X - X / 2;
		bool doit = sqrt((double)x * x + y * y) < X / 2;
		if (doit) {
			out[p] = bi_cube(in, x, y, X, angle);

		}
		else out[p] = 0;
	}
	return out;
}
double PSNR_rot(unsigned char* in, int X, double angle) {
	double* rot = rotate(in, X, angle);
	double* rot_bk = rotate(rot, X, -angle);
	double ret = get_PSNR_(rot_bk, in, X * X);
	delete[]rot_bk;
	delete[] rot;
	return ret;
}
void PSNR_rot(double *psnr_img, double *psnr_sino, unsigned char* in, int X, int A, double angle,  float* sino_gt) {
	double* rot = rotate(in, X, angle);
	double* rot_bk = rotate(rot, X, -angle);
	*psnr_img = get_PSNR_(rot_bk, in, X * X);
	double* sino_rot = forward_result(X, A, rot_bk);
	*psnr_sino = get_PSNR(sino_rot, sino_gt, X * A);
	delete[] sino_rot; 
	delete[]rot_bk;
	delete[] rot;

}
// Calculates the median value in a window around pixel (x, y)
double mediana_window(double* in, int x, int y, int X, int Y, int wx, int wy) {
	// Boundary check
	if (x < 0 || y < 0 || x >= X || y >= Y) {
		return 0.0f;
	}

	// Vector to store values in the neighborhood window
	std::vector<float> values;
	values.reserve(wx * wy);

	// Half window sizes for centering
	int half_wx = wx / 2;
	int half_wy = wy / 2;

	// Collect values within the window
	for (int dy = -half_wy; dy <= half_wy; dy++) {
		for (int dx = -half_wx; dx <= half_wx; dx++) {
			int curX = x + dx;
			int curY = y + dy;

			// Check if current coordinates are within image bounds
			if (curX >= 0 && curX < X && curY >= 0 && curY < Y) {
				values.push_back(in[curX + curY * X]);
			}
		}
	}

	// If window is empty, return the original pixel value
	if (values.empty()) {
		return in[x + y * X];
	}

	// Sort values to find the median
	std::sort(values.begin(), values.end());

	// Calculate median value
	size_t size = values.size();
	if (size % 2 == 0) {
		// Even number of elements: average of two middle values
		return (values[size / 2 - 1] + values[size / 2]) / 2.0f;
	}
	else {
		// Odd number of elements: middle value
		return values[size / 2];
	}
}

float mediana_window(float* in, int x, int y, int X, int Y, int wx, int wy) {
	// Boundary check
	if (x < 0 || y < 0 || x >= X || y >= Y) {
		return 0.0f;
	}

	// Vector to store values in the neighborhood window
	std::vector<float> values;
	values.reserve(wx * wy);

	// Half window sizes for centering
	int half_wx = wx / 2;
	int half_wy = wy / 2;

	// Collect values within the window
	for (int dy = -half_wy; dy <= half_wy; dy++) {
		for (int dx = -half_wx; dx <= half_wx; dx++) {
			int curX = x + dx;
			int curY = y + dy;

			// Check if current coordinates are within image bounds
			if (curX >= 0 && curX < X && curY >= 0 && curY < Y) {
				values.push_back(in[curX + curY * X]);
			}
		}
	}

	// If window is empty, return the original pixel value
	if (values.empty()) {
		return in[x + y * X];
	}

	// Sort values to find the median
	std::sort(values.begin(), values.end());

	// Calculate median value
	size_t size = values.size();
	if (size % 2 == 0) {
		// Even number of elements: average of two middle values
		return (values[size / 2 - 1] + values[size / 2]) / 2.0f;
	}
	else {
		// Odd number of elements: middle value
		return values[size / 2];
	}
}

// Applies median filter to the entire image
//Zero - order hold
void mediana(double* in, int X, int Y, int wx, int wy) {
	// Create output buffer
	double* out = new double[X * Y];

	// Process each pixel in the image
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			out[x + y * X] = mediana_window(in, x, y, X, Y, wx, wy);
		}
	}

	// Copy filtered result back to input array
	memcpy(in, out, sizeof(double) * X * Y);

	// Clean up
	delete[] out;
}
void mediana(float* in, int X, int Y, int wx, int wy) {
	// Create output buffer
	float* out = new float[X * Y];

	// Process each pixel in the image
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			out[x + y * X] = mediana_window(in, x, y, X, Y, wx, wy);
		}
	}

	// Copy filtered result back to input array
	memcpy(in, out, sizeof(float) * X * Y);

	// Clean up
	delete[] out;
}
int get_disk() {
	int ret;
	std::cout << "On which drive is the DIT_BIN folder located ?"
		<< std::endl
		<< std::endl;
	std::cout << "C: - 0"
		<< std::endl;
	std::cout << "D: - 1"
		<< std::endl;
	std::cout << "E: - 2"
		<< std::endl;
	std::cout << "F: - 3"
		<< std::endl
		<< std::endl;
	std::cin >> ret;
	std::cout << std::endl;
	return ret;
}
double get_noise_n() {
	int ret;
	std::cout
		<< "Enter the Gaussian noise level (in %) for the original sinogram:"
		<< std::endl
		<< std::endl;
	std::cout << "0.5%     - 0"
		<< std::endl;
	std::cout << "  1%     - 1"
		<< std::endl;
	std::cout << "  2%     - 2"
		<< std::endl;
	std::cout << "  3%     - 3"
		<< std::endl;
	std::cout << "  n%     - n"
		<< std::endl;
	std::cout << std::endl;
	std::cin >> ret;
	std::cout << std::endl;
	if (ret > 20)ret = 20;

	return ret ? ret * 0.01 : 0.005;
}
int get_angle_n() {
	int ret;
	int rt[6] = { 16,40, 80, 200,400, 800 };
	std::cout
		<< "What number of samples do you select for the [0, 180) degree full angle?"
		<< std::endl
		<< std::endl;
	std::cout << "16(2%)    - 0"
		<< std::endl;
	std::cout << "40(5%)    - 1"
		<< std::endl;
	std::cout << "80(10%)   - 2"
		<< std::endl;
	std::cout << "200(25%)  - 3"
		<< std::endl;
	std::cout << "400(50%)  - 4"
		<< std::endl;
	std::cout << "800(100%) - 5"
		<< std::endl;
	std::cout << "any - more than 5"
		<< std::endl
		<< std::endl;
	std::cin >> ret;
	std::cout << std::endl;
	if (ret <= 5)ret = rt[ret];
	
	return ret;
}
int get_method() {
	int ret;
	std::cout << "What reconstruction technique do you wish to apply?"
		<< std::endl
		<< std::endl;
	std::cout << "Direct Integration               - 0"
		<< std::endl;
	std::cout << "Mean-normalized FBP              - 1"
		<< std::endl;
	std::cout << "Mean and variance normalized FBP - 2"
		<< std::endl;
	std::cout << "Run all (sequentially)           - 3"
		<< std::endl;

	std::cout << std::endl;

	std::cin >> ret;
	std::cout << std::endl;
	if (ret > 3)ret = 3;
	return ret;
}
int get_noise_flt() {
	int ret;
	std::cout << "Denoising filter options:"
		<< std::endl
		<< std::endl;
	std::cout << "None                              - 0"
		<< std::endl;
	std::cout << "Pre-processing (Sinogram domain)  - 1"
		<< std::endl;
	std::cout << "Post-processing (Image domain)    - 2"
		<< std::endl;
	std::cout << "Combined                          - 3"
		<< std::endl;

	std::cout << std::endl;

	std::cin >> ret;
	std::cout << std::endl;
	if (ret > 3)ret = 3;
	return ret;
}
int get_interp() {
	int ret;
	std::cout << "What interpolation  technique do you wish to apply?"
		<< std::endl
		<< std::endl;
	std::cout << "Nearest-neighbor  - 0"
		<< std::endl;
	std::cout << "Linear            - 1"
		<< std::endl;
	std::cout << "Cubic Spline      - 2"
		<< std::endl
		<< std::endl;

	std::cin >> ret;
	std::cout << std::endl;
	if (ret > 2)ret = 2;
	return ret;
}
int get_filter() {
	int ret;
	std::cout << "Which filter would you like to use?"
		<< std::endl
		<< std::endl;
	std::cout << "Ramp filter  - 0"
		<< std::endl;

	std::cout << "Shepp-Logan filter  - 1"
		<< std::endl
		<< std::endl;

	std::cin >> ret;
	std::cout << std::endl;
	if (ret > 1)ret = 1;
	return ret;
}
int  get_image_n() {
	int ret;
	const char* Images[10] =
	{ "Arnolfini", "butterfly", "Hunters",
		 "lizard", "seashell", "still-life", "disk","Shepp-Logan", 
		"Dance", "Time" };

	std::cout << "Which test image would you like to select? "
		<< std::endl
		<< std::endl;
	std::cout << Images[0]<<"        -  0"
		<< std::endl;
	std::cout << Images[1] << "        -  1"
		<< std::endl;
	std::cout << Images[2] << "          -  2"
		<< std::endl;
	std::cout << Images[3] << "           -  3"
		<< std::endl;
	std::cout << Images[4] << "         -  4"
		<< std::endl;
	std::cout << Images[5] << "       -  5"
		<< std::endl;
	std::cout << Images[6] << "             -  6"
		<< std::endl;
	std::cout << Images[7] << "      -  7"
		<< std::endl;
	std::cout << Images[8] << "            -  8"
		<< std::endl;
	std::cout << Images[9] << "             -  9"
		<< std::endl;
	std::cout << "All of the above" << " - 10"
		<< std::endl
		<< std::endl;
	
	std::cin >> ret;
	std::cout << std::endl;
	if (ret > 10)ret = 10;
	return ret;
}
const char* gen_folder(int disk, int method) {
	if (disk > 5)disk = 5;
	if (method > 2)method = 2;
	const char* Folders[3] =
	{ "\\DIT_BIN\\DIT_TEST_IMAGES", "\\DIT_BIN\\FBP_TEST_IMAGES",
		"\\DIT_BIN\\FBP_NORM_TEST_IMAGES" };
	const char* Dsk[6] =
	{ "C:", "D:","E:", "G:","F:", "H:" };

	return (const char*)concat(Dsk[disk], Folders[method]);
}
void waitEnter() {
	// Clear buffer if there's leftover input from previous cin >>
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cout << "Press Enter to continue...";
	std::cin.get();
}


