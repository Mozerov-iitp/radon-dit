#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iterator>
#include <iostream>
#include "inverse_fourier_transform.h"

/**
 * Performs cubic interpolation using 4 equidistant sample points.
 * Implements the Catmull-Rom spline formula.
 * @param t - Interpolation parameter in [0, 1] between the middle two points.
 * @param f - Array of 4 function values at x = [-1, 0, 1, 2].
 * @return Interpolated value at position t.
 */
double zero_min_set(int size, double* in) {
    double ret = 0;

    double* out = new double[size];
    int p = size;
    double min = in[p - 1];
    while (p--) {

        min = min > in[p] ? in[p] : min;
    }


    p = size; while (p--) {

        ret += out[p] = (in[p] - min);
    }
    memcpy(in, out, sizeof(double) * size);
    delete[] out;
    return ret / size;
}

double cubic_interpolate(double t, const double* f) {
    // Coefficients for the cubic polynomial: a0 + a1*t + a2*t^2 + a3*t^3
    double a0 = f[1];
    double a1 = 0.5 * (f[2] - f[0]); // Central difference approx. of 1st derivative
    double a2 = f[0] - 2.5 * f[1] + 2.0 * f[2] - 0.5 * f[3];
    double a3 = 0.5 * (f[3] - f[0]) + 1.5 * (f[1] - f[2]);
    return a0 + t * (a1 + t * (a2 + t * a3));
}

/**
 * Performs cubic interpolation on a 1D signal with boundary checking.
 * @param signal - Input signal array.
 * @param index - Floating-point index. interpolation is performed between
 *                indices floor(index) and floor(index)+1.
 * @param size - Total size of the signal array.
 * @return Interpolated value, or 0.0 if out of bounds.
 */
double interpolate_at(const double* signal, double index, int size) {
    int base_idx = static_cast<int>(floor(index));
    double t = index - base_idx; // Fractional part

    // Prepare 4 points for cubic interpolation: [base_idx-1, base_idx, base_idx+1, base_idx+2]
    double values[4];
    for (int i = 0; i < 4; ++i) {
        int fetch_idx = base_idx - 1 + i; // Indices: -1, 0, +1, +2 relative to base
        // Handle boundaries: clamp to nearest valid index
        if (fetch_idx < 0) fetch_idx = 0;
        else if (fetch_idx >= size) fetch_idx = size - 1;
        values[i] = signal[fetch_idx];
    }

    return cubic_interpolate(t, values);
}

/**
 * Reconstructs a single pixel using Filtered Back Projection.
 * @param filtered_sinogram - Ramp-filtered sinogram (size: num_angles * detector_size).
 * @param pixel_idx - Linear index of the pixel in the reconstructed image.
 * @param image_size - Width/height of the square output image.
 * @param detector_size - Number of detector elements.
 * @param num_angles - Number of projection angles.
 * @return Reconstructed pixel value.
 */
double reconstruct_pixel(double* filtered_sinogram, int pixel_idx,
    int image_size, int detector_size, int num_angles) {
    // Convert linear index to 2D image coordinates with origin at center
    int x = pixel_idx % image_size - image_size / 2;
    int y = pixel_idx / image_size - image_size / 2;

    // Check if pixel is within the reconstruction circle (inscribed in square)
    double radius = sqrt(x * x + y * y);


    double angle_step = M_PI / num_angles; // Angular sampling step
    double sum = 0.0;

    // Back-project over all angles
    for (int angle_idx = 0; angle_idx < num_angles; ++angle_idx) {
        double angle = angle_idx * angle_step;

        // Calculate detector coordinate for this pixel and angle
        double detector_coord = x * cos(angle) + y * sin(angle); // Radon transform
        // Convert to array index: shift from [-D/2, D/2] to [0, D-1]
        double array_index = detector_coord + detector_size / 2.0;

        // Perform cubic interpolation in the filtered projection
        int projection_offset = angle_idx * detector_size;
        double interpolated_value = interpolate_at(&filtered_sinogram[projection_offset],
            array_index, detector_size);
        sum += interpolated_value;
    }

    // Normalization is applied AFTER all back-projections (outside this function)
    return sum < 0 ? 0 : sum;
}

/**
 * Main FBP reconstruction algorithm with optional filters.
 * @param sinogram - Input sinogram data (size: num_angles * detector_size).
 * @param detector_size - Number of detector elements.
 * @param num_angles - Number of projection angles (uniformly sampled over 180°).
 * @param filter_type - Type of filter to use: 0 = ramp, 1 = Shepp-Logan
 * @return Reconstructed image (size: detector_size * detector_size).
 *         NOTE: Caller is responsible for freeing the memory.
 */
void cout_end_time(float time) {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "================  FBP End  ================" << std::endl;
    std::cout << std::endl;
    std::cout << "        Execution time: " << time << " seconds" << std::endl;

    std::cout << std::endl;
    std::cout << "===========================================" << std::endl;
}

double* fbp_reconstruct(float* sinogram, int detector_size,
    int num_angles, int filter_type) {

    std::cout << "=======            FBP               ======" << std::endl;
    std::cout << std::endl;
    clock_t start = clock();
    double* filter_weight = new double[detector_size];
    double weight_norm = detector_size / 4;
    int X = detector_size;
    for (int x = 0; x < X; x++) {
        int frequency = (x < X / 2 / 2)
            ? x
            : x - X;
        if (filter_type == 0) {
            filter_weight[x] = abs(frequency);
        }

        else { // Shepp-Logan filter
            if (!frequency) {
                // Avoid division by zero at frequency = 0
                filter_weight[x] = 0;
            }
            else {
                // Shepp-Logan filter: |ω| * sinc(ω/2)
                // sinc(x) = sin(πx)/(πx)
                double omega = 2.0 * M_PI * frequency / X; // Angular frequency

                filter_weight[x] = abs(frequency) *
                    sin(omega * 0.5) / (omega * 0.5);
            }
        }
    }
    for (int x = 0; x < X; x++)filter_weight[x] /= weight_norm;

    int total_sino_samples = num_angles * detector_size;
    double mean_from_sino = 0;
    int p = total_sino_samples;
    while (p--)mean_from_sino += sinogram[p];
    mean_from_sino /= num_angles * detector_size * detector_size;

 

    // 1. Allocate buffers
    complex_dbl* fft_buffer = new complex_dbl[total_sino_samples];
    double* filtered_sinogram = new double[total_sino_samples];

    // 2. 1D FFT of each projection (row)
    for (int angle_idx = 0; angle_idx < num_angles; ++angle_idx) {
        float* projection_in = &sinogram[angle_idx * detector_size];
        complex_dbl* fft_row = &fft_buffer[angle_idx * detector_size];
        fft_X(projection_in, fft_row, detector_size); // Your FFT function
    }

    // 3. Apply filter in frequency domain
    for (int angle_idx = 0; angle_idx < num_angles; ++angle_idx) {
        for (int freq_idx = 0; freq_idx < detector_size; ++freq_idx) {



            complex_dbl& spectrum_val =
                fft_buffer[angle_idx * detector_size + freq_idx];
            spectrum_val *= filter_weight[freq_idx];
        }
    }

    // 4. Inverse FFT to get filtered sinogram
    for (int angle_idx = 0; angle_idx < num_angles; ++angle_idx) {
        std::complex<double>* fft_row = &fft_buffer[angle_idx * detector_size];
        double* filtered_row = &filtered_sinogram[angle_idx * detector_size];
        fft_inv_X(fft_row, filtered_row, detector_size); // Your IFFT function
    }

    // 5. CORRECT NORMALIZATION: Divide by number of angles
    for (int i = 0; i < total_sino_samples; ++i) {
        filtered_sinogram[i] /= num_angles;
    }

    // 6. Allocate output image
    int image_size = detector_size; // Typically image size equals detector size
    int total_pixels = image_size * image_size;
    double* reconstructed_image = new double[total_pixels];

    // 7. Back-projection for each pixel
    int mod_T = total_pixels / 50;
    double mean_from_rec_image = 0;
    for (int pixel_idx = 0; pixel_idx < total_pixels; ++pixel_idx) {

        mean_from_rec_image +=
            reconstructed_image[pixel_idx] = reconstruct_pixel(
                filtered_sinogram, pixel_idx,
                image_size, detector_size, num_angles
            );

        if ((pixel_idx) % mod_T == 0) { std::cout << char(219); }
    }
    mean_from_rec_image /= total_pixels;

    clock_t finish = clock();
    float time = (float)(((double)(finish - start)) / CLOCKS_PER_SEC);
    cout_end_time(time);




    p = total_pixels;
    while (p--)
    reconstructed_image[p] += -mean_from_rec_image + mean_from_sino;

    // Clean up temporary buffers
    delete[] fft_buffer;
    delete[] filtered_sinogram;

    return reconstructed_image;
}

// Overloaded version for backward compatibility (uses ramp filter by default)
//double* fbp_reconstruct(float* sinogram, int detector_size, int num_angles) {
//    return fbp_reconstruct(sinogram, detector_size, num_angles, 0);
//}