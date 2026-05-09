#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>


constexpr double PI = 3.141592653589793238462643383279502884;

// ============================================================================
// BICUBIC INTERPOLATION FUNCTIONS (DOUBLE PRECISION INTERNALLY)
// ============================================================================

/**
 * Cubic Catmull-Rom interpolation for 1D (double precision)
 * @param p Array of 4 control points (double)
 * @param x Interpolation position between p[1] and p[2] (0 to 1)
 * @return Interpolated value (double)
 */
//double bilateral_flt(unsigned char* in, int X, int Y, int x, int y, int W,
//    double sigma_r, double sigma_c) {
//    // Boundary check
//    if (x < 0 || x >= X || y < 0 || y >= Y) return 0.0;
//
//    // Window radius (half of width)
//    int radius = W / 2;
//
//    // Current pixel value (center pixel)
//    double center_value = static_cast<double>(in[y * X + x]);
//
//    // Accumulation variables
//    double sum = 0.0;
//    double weight_sum = 0.0;
//
//    // Precompute constants for efficiency
//    // sigma_r - spatial component, sigma_c - range/intensity component
//    double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;  // Spatial component
//    double two_sigma_c_sq = 2.0 * sigma_c * sigma_c;  // Range component
//
//    // Avoid division by zero issues
//    if (two_sigma_r_sq < 1e-10) two_sigma_r_sq = 1e-10;
//    if (two_sigma_c_sq < 1e-10) two_sigma_c_sq = 1e-10;
//
//    // Iterate over the neighborhood window
//    for (int j = -radius; j <= radius; ++j) {
//        int current_y = y + j;
//
//        // Check Y boundary
//        if (current_y < 0 || current_y >= Y) continue;
//
//        for (int i = -radius; i <= radius; ++i) {
//            int current_x = x + i;
//
//            // Check X boundary
//            if (current_x < 0 || current_x >= X) continue;
//
//            // Get neighbor pixel value
//            double neighbor_value = static_cast<double>(in[current_y * X + current_x]);
//
//            // Compute SPATIAL weight (distance-based Gaussian)
//            // sigma_r controls spatial influence
//            double spatial_dist_sq = static_cast<double>(i * i + j * j);
//            double spatial_weight = exp(-spatial_dist_sq / two_sigma_r_sq);
//
//            // Compute RANGE weight (intensity difference-based Gaussian)
//            // sigma_c controls intensity similarity
//            double range_diff = neighbor_value - center_value;
//            double range_dist_sq = range_diff * range_diff;
//            double range_weight = exp(-range_dist_sq / two_sigma_c_sq);
//
//            // Combined weight = spatial weight × range weight
//            double total_weight = spatial_weight * range_weight;
//
//            // Accumulate weighted sum and total weight
//            sum += neighbor_value * total_weight;
//            weight_sum += total_weight;
//        }
//    }
//
//    // Avoid division by zero (if all weights are zero)
//    if (weight_sum < 1e-10) {
//        return center_value;
//    }
//
//    // Return the normalized weighted average
//    return sum / weight_sum;
//}
//unsigned char* applyBilateralFilter(unsigned char* input, 
//    int width, int height, int window_size,
//    double spatial_sigma, double range_sigma) {
//
//    // Allocate output if needed
//    unsigned char* output = new unsigned char[width * height];
//
//    for (int y = 0; y < height; ++y) {
//        for (int x = 0; x < width; ++x) {
//            // Apply bilateral filter to each pixel
//            double filtered_value = bilateral_flt(input, width, height, x, y,
//                window_size, spatial_sigma, range_sigma);
//
//            // Clamp to valid range [0, 255] and round
//            filtered_value = std::max(0.0, std::min(255.0, filtered_value));
//            output[y * width + x] = static_cast<unsigned char>(filtered_value + 0.5);
//        }
//    }
//    return output;
//}

double cubicInterpolate(double p[4], double x) {
    // Catmull-Rom spline formula (double precision):
    // f(x) = 0.5 * [(-p0 + 3p1 - 3p2 + p3)x³ + (2p0 - 5p1 + 4p2 - p3)x² + (-p0 + p2)x + 2p1]
    double ret =  p[1] + 0.5 * x * (p[2] - p[0] +
        x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
            x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));

    return ret;
}

/**
 * Bicubic interpolation for 2D image with double precision internal calculations
 * @param image Input image as 1D array (float, but converted to double for calculations)
 * @param width Image width
 * @param height Image height
 * @param x X coordinate (double precision)
 * @param y Y coordinate (double precision)
 * @return Interpolated value as float (converted from double)
 */
float bicubicInterpolateDouble(const std::vector<float>& image, int width, int height,
    double x, double y) {
    // Get integer coordinates of the central pixel (using double floor)
    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));

    // Check if central pixel is within bounds
    if (x0 < 0 || y0 < 0 || x0 >= width || y0 >= height) {
        return 0.0f;  // Return 0 for out-of-bounds
    }

    // Collect 4x4 neighborhood for interpolation (stored as double)
    double neighborhood[4][4];

    // Fractional part for interpolation
    double dx = x - x0;
    double dy = y - y0;

    // Ensure dx and dy are in [0, 1) range
    if (dx < 0.0) dx = 0.0;
    if (dx >= 1.0) dx = 0.999999;
    if (dy < 0.0) dy = 0.0;
    if (dy >= 1.0) dy = 0.999999;

    for (int j = -1; j <= 2; j++) {      // -1, 0, 1, 2 (4 rows)
        for (int i = -1; i <= 2; i++) {  // -1, 0, 1, 2 (4 columns)
            int xi = x0 + i;
            int yj = y0 + j;

            double pixelValue = 0.0;

            if (xi >= 0 && xi < width && yj >= 0 && yj < height) {
                // Pixel is inside image bounds - convert float to double
                pixelValue = static_cast<double>(image[yj * width + xi]);
            }
            else {
                // Mirror boundary conditions
                int xc = xi;
                int yc = yj;

                if (xc < 0) xc = -xc;
                else if (xc >= width) xc = 2 * width - xc - 1;

                if (yc < 0) yc = -yc;
                else if (yc >= height) yc = 2 * height - yc - 1;


                pixelValue = static_cast<double>(image[yc * width + xc]);
            }

            neighborhood[j + 1][i + 1] = pixelValue;
        }
    }

    // First, interpolate along X axis for each of the 4 rows (double precision)
    double interpolatedRows[4];
    for (int row = 0; row < 4; row++) {
        double rowValues[4] = {
            neighborhood[row][0],
            neighborhood[row][1],
            neighborhood[row][2],
            neighborhood[row][3]
        };
        interpolatedRows[row] = cubicInterpolate(rowValues, dx);
    }

    // Then interpolate along Y axis using the 4 interpolated row values
    double result = cubicInterpolate(interpolatedRows, dy);

    // Convert final result to float (output precision)
    return static_cast<float>(result);
}
float bicubicInterpolateDouble(const std::vector<double>& image, int width, int height,
    double x, double y) {
    // Get integer coordinates of the central pixel (using double floor)
    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));

    // Check if central pixel is within bounds
    if (x0 < 0 || y0 < 0 || x0 >= width || y0 >= height) {
        return 0.0f;  // Return 0 for out-of-bounds
    }

    // Collect 4x4 neighborhood for interpolation (stored as double)
    double neighborhood[4][4];

    // Fractional part for interpolation
    double dx = x - x0;
    double dy = y - y0;

    // Ensure dx and dy are in [0, 1) range
    if (dx < 0.0) dx = 0.0;
    if (dx >= 1.0) dx = 0.999999;
    if (dy < 0.0) dy = 0.0;
    if (dy >= 1.0) dy = 0.999999;

    for (int j = -1; j <= 2; j++) {      // -1, 0, 1, 2 (4 rows)
        for (int i = -1; i <= 2; i++) {  // -1, 0, 1, 2 (4 columns)
            int xi = x0 + i;
            int yj = y0 + j;

            double pixelValue = 0.0;

            if (xi >= 0 && xi < width && yj >= 0 && yj < height) {
                // Pixel is inside image bounds - convert float to double
                pixelValue = static_cast<double>(image[yj * width + xi]);
            }
            else {
                // Mirror boundary conditions
                int xc = xi;
                int yc = yj;

                if (xc < 0) xc = -xc;
                else if (xc >= width) xc = 2 * width - xc - 1;

                if (yc < 0) yc = -yc;
                else if (yc >= height) yc = 2 * height - yc - 1;


                pixelValue = static_cast<double>(image[yc * width + xc]);
            }

            neighborhood[j + 1][i + 1] = pixelValue;
        }
    }

    // First, interpolate along X axis for each of the 4 rows (double precision)
    double interpolatedRows[4];
    for (int row = 0; row < 4; row++) {
        double rowValues[4] = {
            neighborhood[row][0],
            neighborhood[row][1],
            neighborhood[row][2],
            neighborhood[row][3]
        };
        interpolatedRows[row] = cubicInterpolate(rowValues, dx);
    }

    // Then interpolate along Y axis using the 4 interpolated row values
    double result = cubicInterpolate(interpolatedRows, dy);

    // Convert final result to float (output precision)
    return static_cast<float>(result);
}


// ============================================================================
// DISCRETE RADON TRANSFORM (DOUBLE PRECISION INTERNALLY)
// ============================================================================

/**
 * Discrete Radon Transform with bicubic interpolation
 * INTERNAL COMPUTATIONS: double precision
 * OUTPUT: float precision (for memory efficiency)
 *
 * @param image Input image (square, row-major, float)
 * @param imageSize Size of image (width = height = imageSize)
 * @param numAngles Number of projection angles
 * @param numOffsets Number of offsets (detector positions)
 * @return Sinogram as 1D float array [angle][offset]
 */
std::vector<float> prepareForRadon(const unsigned char* img, int size) {
    int n = size * size;
    std::vector<float> out(n);
    for (int i = 0; i < n; i++) out[i] = (float)img[i];
    return out;
}
std::vector<double> prepareForRadon(double* img, int size) {
    int n = size * size;
    std::vector<double> out(n);
    for (int i = 0; i < n; i++) out[i] = (double)img[i];
    return out;
}
float* discreteRadonTransformDoublePtr(const std::vector<float>& image,
    int imageSize,
    int numAngles,
    int numOffsets) {
    // Validate input
    if (image.size() != static_cast<size_t>(imageSize * imageSize)) {
        throw std::invalid_argument("Image must be square and match imageSize");
    }

    // Allocate output array (raw pointer)
    float* sinogram = new float[numAngles * numOffsets];

    // Initialize with zeros
    std::fill(sinogram, sinogram + numAngles * numOffsets, 0.0f);

    // Image center coordinates (DOUBLE PRECISION)
    double center = static_cast<double>(imageSize) / 2.0;

    // Angle step: from 0 to π (180 degrees)
    double angleStep = PI / static_cast<double>(numAngles);

    // Calculate maximum distance from center
   // double maxDist = std::sqrt(2.0) * center;

    // Offset step between detector positions
   // double offsetStep = (2.0 * maxDist) / static_cast<double>(numOffsets - 1);
    double offsetStep = 1.0;
    double maxDist = offsetStep * (imageSize) / 2.0;
    // 
    std::cout << "------Starting Radon transform computation------" << std::endl;
    std::cout << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    // Precompute trigonometric values for all angles
    std::vector<double> cosAngles(numAngles);
    std::vector<double> sinAngles(numAngles);
    for (int angleIdx = 0; angleIdx < numAngles; angleIdx++) {
        double angle = static_cast<double>(angleIdx) * angleStep;
        cosAngles[angleIdx] = std::cos(angle);
        sinAngles[angleIdx] = std::sin(angle);
    }

    // Precompute offsets
    std::vector<double> offsets(numOffsets);
    for (int offsetIdx = 0; offsetIdx < numOffsets; offsetIdx++) {
        offsets[offsetIdx] = -maxDist + static_cast<double>(offsetIdx) * offsetStep;
    }

    // Main computation loop
//#pragma omp parallel for
    int st_T = numAngles / 40;
    if (st_T < 1)st_T = 1;
    for (int angleIdx = 0; angleIdx < numAngles; angleIdx++) {
        double cosAngle = cosAngles[angleIdx];
        double sinAngle = sinAngles[angleIdx];

        for (int offsetIdx = 0; offsetIdx < numOffsets; offsetIdx++) {
            double offset = offsets[offsetIdx];
            double integral = 0.0;

            // Constant step in world coordinates
            double tStep = 1.;
            double tMin = -static_cast<double>(imageSize) * 1.5;
            double tMax = static_cast<double>(imageSize) * 1.5;

            int numSteps = static_cast<int>((tMax - tMin) / tStep);

            for (int step = 0; step <= numSteps; step++) {
                double t = tMin + static_cast<double>(step) * tStep;

                // Calculate world coordinates
                double worldX = -t * sinAngle + offset * cosAngle;
                double worldY = t * cosAngle + offset * sinAngle;

                // Convert to image coordinates
                double imgX = worldX + center;
                double imgY = worldY + center;

                // Check bounds
                if (imgX >= 0.0 && imgX < static_cast<double>(imageSize - 1) &&
                    imgY >= 0.0 && imgY < static_cast<double>(imageSize - 1)) {

                    // Bicubic interpolation
                    float value = bicubicInterpolateDouble(image, imageSize, imageSize, imgX, imgY);
                    integral += static_cast<double>(value) * tStep;
                }
            }

            // Store result (convert double to float)
            sinogram[angleIdx * numOffsets + offsetIdx] = static_cast<float>(integral);
        }
        if (angleIdx % st_T == 0) { std::cout << char(219); }
    } std::cout << std::endl; std::cout << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::
        duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Radon transform completed in " << duration.count()/1000 << " seconds" << std::endl;
    std::cout << std::endl;
    return sinogram;
}
double* discreteRadonTransformDoublePtr(const std::vector<double>& image,
    int imageSize,
    int numAngles,
    int numOffsets) {
    // Validate input
    if (image.size() != static_cast<size_t>(imageSize * imageSize)) {
        throw std::invalid_argument("Image must be square and match imageSize");
    }

    // Allocate output array (raw pointer)
   double* sinogram = new double[numAngles * numOffsets];

    // Initialize with zeros
    std::fill(sinogram, sinogram + numAngles * numOffsets, 0.0f);

    // Image center coordinates (DOUBLE PRECISION)
    double center = static_cast<double>(imageSize) / 2.0;

    // Angle step: from 0 to π (180 degrees)
    double angleStep = PI / static_cast<double>(numAngles);

    // Calculate maximum distance from center
   // double maxDist = std::sqrt(2.0) * center;

    // Offset step between detector positions
   // double offsetStep = (2.0 * maxDist) / static_cast<double>(numOffsets - 1);
    double offsetStep = 1.0;
    double maxDist = offsetStep * (imageSize) / 2.0;
    // 
    std::cout << "- Starting Radon transform computation of the yet reonstructed image -" << std::endl;
    std::cout << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    // Precompute trigonometric values for all angles
    std::vector<double> cosAngles(numAngles);
    std::vector<double> sinAngles(numAngles);
    for (int angleIdx = 0; angleIdx < numAngles; angleIdx++) {
        double angle = static_cast<double>(angleIdx) * angleStep;
        cosAngles[angleIdx] = std::cos(angle);
        sinAngles[angleIdx] = std::sin(angle);
    }

    // Precompute offsets
    std::vector<double> offsets(numOffsets);
    for (int offsetIdx = 0; offsetIdx < numOffsets; offsetIdx++) {
        offsets[offsetIdx] = -maxDist + static_cast<double>(offsetIdx) * offsetStep;
    }

    // Main computation loop
//#pragma omp parallel for
    int st_T = numAngles / 40;
    if (st_T < 1)st_T = 1;
    for (int angleIdx = 0; angleIdx < numAngles; angleIdx++) {
        double cosAngle = cosAngles[angleIdx];
        double sinAngle = sinAngles[angleIdx];

        for (int offsetIdx = 0; offsetIdx < numOffsets; offsetIdx++) {
            double offset = offsets[offsetIdx];
            double integral = 0.0;

            // Constant step in world coordinates
            double tStep = 1.0;
            double tMin = -static_cast<double>(imageSize) * 1.5;
            double tMax = static_cast<double>(imageSize) * 1.5;

            int numSteps = static_cast<int>((tMax - tMin) / tStep);

            for (int step = 0; step <= numSteps; step++) {
                double t = tMin + static_cast<double>(step) * tStep;

                // Calculate world coordinates
                double worldX = -t * sinAngle + offset * cosAngle;
                double worldY = t * cosAngle + offset * sinAngle;

                // Convert to image coordinates
                double imgX = worldX + center;
                double imgY = worldY + center;

                // Check bounds
                if (imgX >= 0.0 && imgX < static_cast<double>(imageSize - 1) &&
                    imgY >= 0.0 && imgY < static_cast<double>(imageSize - 1)) {

                    // Bicubic interpolation
                    float value = bicubicInterpolateDouble(image, imageSize, imageSize, imgX, imgY);
                    integral += static_cast<double>(value) * tStep;
                }
            }

            // Store result (convert double to float)
            sinogram[angleIdx * numOffsets + offsetIdx] = static_cast<float>(integral);
        }
        if (angleIdx % st_T == 0) { std::cout << char(219); }
    } std::cout << std::endl; std::cout << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::
        duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Radon transform completed in " << duration.count() / 1000 << " seconds" << std::endl;
    std::cout << std::endl;
    return sinogram;
}
std::vector<float> discreteRadonTransformDouble(const std::vector<float>& image,
    int imageSize,
    int numAngles,
    int numOffsets) {
    // Validate input
    if (image.size() != static_cast<size_t>(imageSize * imageSize)) {
        throw std::invalid_argument("Image must be square and match imageSize");
    }

    // Initialize sinogram (output): numAngles × numOffsets (FLOAT OUTPUT)
    std::vector<float> sinogram(numAngles * numOffsets, 0.0f);

    // Image center coordinates (DOUBLE PRECISION)
    double center = static_cast<double>(imageSize - 1) / 2.0;

    // Angle step: from 0 to π (180 degrees) (DOUBLE PRECISION)
    double angleStep = PI / static_cast<double>(numAngles);

    // Calculate maximum distance from center to cover entire rotated image
    // For a square image rotated by 45°, the farthest point is at distance √2 * center
    double maxDist = std::sqrt(2.0) * center;

    // Offset step between detector positions (DOUBLE PRECISION)
    double offsetStep = (2.0 * maxDist) / static_cast<double>(numOffsets - 1);

    std::cout << "Starting Radon transform computation (double precision internally)..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // Precompute trigonometric values for all angles (DOUBLE PRECISION)
    std::vector<double> cosAngles(numAngles);
    std::vector<double> sinAngles(numAngles);
    for (int angleIdx = 0; angleIdx < numAngles; angleIdx++) {
        double angle = static_cast<double>(angleIdx) * angleStep;
        cosAngles[angleIdx] = std::cos(angle);
        sinAngles[angleIdx] = std::sin(angle);
    }

    // Precompute offsets (DOUBLE PRECISION)
    std::vector<double> offsets(numOffsets);
    for (int offsetIdx = 0; offsetIdx < numOffsets; offsetIdx++) {
        offsets[offsetIdx] = -maxDist + static_cast<double>(offsetIdx) * offsetStep;
    }

    // ========================================================================
    // MAIN COMPUTATION LOOP (DOUBLE PRECISION INTERNALLY)
    // ========================================================================
    // Parallelize over angles (each angle is independent)
#pragma omp parallel for
    for (int angleIdx = 0; angleIdx < numAngles; angleIdx++) {
        // Current trigonometric values (DOUBLE PRECISION)
        double cosAngle = cosAngles[angleIdx];
        double sinAngle = sinAngles[angleIdx];

        // For each detector position (offset)
        for (int offsetIdx = 0; offsetIdx < numOffsets; offsetIdx++) {
            // Offset s in line equation: x·cosθ + y·sinθ = s (DOUBLE PRECISION)
            double offset = offsets[offsetIdx];

            // Accumulator for integral (DOUBLE PRECISION)
            double integral = 0.0;

            // ================================================================
            // INTEGRATION ALONG THE LINE (DOUBLE PRECISION)
            // ================================================================
            // Line parameterization:
            // x(t) = -t·sinθ + s·cosθ + center_x
            // y(t) =  t·cosθ + s·sinθ + center_y
            // where t is parameter along the line

            // Constant step in world coordinates (ISOTROPIC SAMPLING)
            // Using double precision for all calculations
            double tStep = 1.0;  // Constant step in world units (DOUBLE)

            // Integration bounds: cover entire diagonal of rotated image
            double tMin = -static_cast<double>(imageSize) * 1.5;  // Slightly larger than necessary
            double tMax = static_cast<double>(imageSize) * 1.5;

            // Calculate number of steps (ensuring we cover the entire range)
            int numSteps = static_cast<int>((tMax - tMin) / tStep);

            // Numerical integration using rectangle rule (DOUBLE PRECISION)
            for (int step = 0; step <= numSteps; step++) {
                double t = tMin + static_cast<double>(step) * tStep;

                // Calculate world coordinates (DOUBLE PRECISION)
                double worldX = -t * sinAngle + offset * cosAngle;
                double worldY = t * cosAngle + offset * sinAngle;

                // Convert to image coordinates (add center offset)
                double imgX = worldX + center;
                double imgY = worldY + center;

                // Check if point is inside or near image bounds
                // We use double precision for bounds checking
                if (imgX >= 0.0 && imgX < static_cast<double>(imageSize - 1) &&
                    imgY >= 0.0 && imgY < static_cast<double>(imageSize - 1)) {

                    // Get interpolated value at (imgX, imgY) - returns float
                    // but internal calculations are double precision
                    float value = bicubicInterpolateDouble(image, imageSize, imageSize, imgX, imgY);

                    // Add to integral (convert value to double for accumulation)
                    // ∫ f(x,y) dt ≈ Σ f(x(tᵢ), y(tᵢ)) · Δt
                    integral += static_cast<double>(value) * tStep;
                }
            }

            // Store result in sinogram (convert double to float for output)
            sinogram[angleIdx * numOffsets + offsetIdx] = static_cast<float>(integral);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Radon transform completed in " << duration.count() << " ms" << std::endl;

    return sinogram;
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * Create Shepp-Logan phantom (standard test image in tomography)
 * Uses double precision for calculations, float for output
 * @param size Image size
 * @return Shepp-Logan phantom image as float
 */
std::vector<float> createSheppLoganPhantomDouble(int size) {
    std::vector<float> image(size * size, 0.0f);
    double center = static_cast<double>(size - 1) / 2.0;

    // Define ellipses: (center_x, center_y, a, b, rotation, value)
    // Parameters from classic Shepp-Logan phantom
    struct Ellipse {
        double cx, cy;  // Center relative to image center (in normalized coordinates -1 to 1)
        double a, b;    // Semi-axes (in normalized coordinates)
        double angle;   // Rotation angle in degrees
        double value;   // Intensity
    };

    // Classic Shepp-Logan parameters (double precision)
    std::vector<Ellipse> ellipses = {
        {0.0, 0.0, 0.69, 0.92, 0.0, 1.0},           // Main ellipse
        {0.0, -0.0184, 0.6624, 0.874, 0.0, -0.8},   // Inner ellipse
        {0.22, 0.0, 0.11, 0.31, -18.0, -0.2},       // Left small ellipse
        {-0.22, 0.0, 0.16, 0.41, 18.0, -0.2},       // Right small ellipse
        {0.0, 0.35, 0.21, 0.25, 0.0, 0.1},          // Top ellipse
        {0.0, 0.1, 0.046, 0.046, 0.0, 0.1},         // Small center ellipse
        {-0.08, -0.605, 0.046, 0.023, 0.0, 0.1},    // Bottom left
        {0.0, -0.605, 0.023, 0.023, 0.0, 0.1},      // Bottom center
        {0.06, -0.605, 0.023, 0.046, 0.0, 0.1}      // Bottom right
    };

    // Convert normalized coordinates to pixel coordinates
    double scale = static_cast<double>(size) / 2.0;

    // Precompute angle conversions
    const double degToRad = PI / 180.0;

    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            // Convert to normalized coordinates [-1, 1] (double precision)
            double nx = (static_cast<double>(x) - center) / scale;
            double ny = (static_cast<double>(y) - center) / scale;

            double pixelValue = 0.0;

            // Sum contributions from all ellipses (double precision calculations)
            for (const auto& ellipse : ellipses) {
                // Translate point to ellipse center
                double tx = nx - ellipse.cx;
                double ty = ny - ellipse.cy;

                // Rotate point by -angle (to align with ellipse axes)
                double radAngle = ellipse.angle * degToRad;
                double cosA = std::cos(radAngle);
                double sinA = std::sin(radAngle);

                double rx = tx * cosA + ty * sinA;
                double ry = -tx * sinA + ty * cosA;

                // Check if point is inside ellipse: (rx/a)² + (ry/b)² ≤ 1
                double term = (rx * rx) / (ellipse.a * ellipse.a) +
                    (ry * ry) / (ellipse.b * ellipse.b);

                if (term <= 1.0) {
                    pixelValue += ellipse.value;
                }
            }

            // Convert to float for output
            image[y * size + x] = static_cast<float>(pixelValue);
        }
    }

    return image;
}

/**
 * Create a simple circle test image with double precision calculations
 * @param size Image size
 * @return Circle image as float
 */
std::vector<float> createCircleTestDouble(int size) {
    std::vector<float> image(size * size, 0.0f);
    double center = static_cast<double>(size - 1) / 2.0;
    double radius = static_cast<double>(size) / 4.0;
    double radiusSq = radius * radius;

    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            double dx = static_cast<double>(x) - center;
            double dy = static_cast<double>(y) - center;
            double distSq = dx * dx + dy * dy;

            // Create a circle with soft edges (double precision calculations)
            if (distSq <= radiusSq) {
                // Simple binary circle
                image[y * size + x] = 1.0f;
            }
        }
    }

    return image;
}

/**
 * Save float data to PGM file
 * @param data Float data
 * @param width Image width
 * @param height Image height
 * @param filename Output filename
 */
void saveFloatToPGM(const std::vector<float>& data, int width, int height,
    const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    // Find min and max for normalization (double precision for accuracy)
    double minVal = static_cast<double>(data[0]);
    double maxVal = static_cast<double>(data[0]);

    for (size_t i = 1; i < data.size(); i++) {
        double val = static_cast<double>(data[i]);
        if (val < minVal) minVal = val;
        if (val > maxVal) maxVal = val;
    }

    double range = maxVal - minVal;
    if (range == 0.0) range = 1.0;

    // Write PGM header
    file << "P5\n" << width << " " << height << "\n255\n";

    // Write pixel data
    for (size_t i = 0; i < data.size(); i++) {
        double normalized = (static_cast<double>(data[i]) - minVal) / range;
        auto pixel = static_cast<unsigned char>(normalized * 255.0);
        file.write(reinterpret_cast<const char*>(&pixel), 1);
    }

    file.close();
    std::cout << "Saved: " << filename << " (range: " << minVal << " to " << maxVal << ")" << std::endl;
}

/**
 * Compare float and double precision results
 * @param sinogramFloat Sinogram computed with float precision
 * @param sinogramDouble Sinogram computed with double precision internally
 * @param tolerance Maximum allowed difference
 */
void comparePrecision(const std::vector<float>& sinogramFloat,
    const std::vector<float>& sinogramDouble,
    float tolerance = 1e-6f) {
    if (sinogramFloat.size() != sinogramDouble.size()) {
        std::cout << "WARNING: Sinograms have different sizes, cannot compare" << std::endl;
        return;
    }

    double maxDiff = 0.0;
    double avgDiff = 0.0;
    int countExceeding = 0;

    for (size_t i = 0; i < sinogramFloat.size(); i++) {
        double diff = std::abs(static_cast<double>(sinogramFloat[i]) -
            static_cast<double>(sinogramDouble[i]));

        avgDiff += diff;
        if (diff > maxDiff) maxDiff = diff;
        if (diff > static_cast<double>(tolerance)) countExceeding++;
    }

    avgDiff /= static_cast<double>(sinogramFloat.size());

    std::cout << "\nPrecision Comparison:" << std::endl;
    std::cout << "  Maximum difference: " << std::scientific << std::setprecision(2) << maxDiff << std::endl;
    std::cout << "  Average difference: " << std::scientific << std::setprecision(2) << avgDiff << std::endl;
    std::cout << "  Pixels exceeding tolerance (" << tolerance << "): "
        << countExceeding << " / " << sinogramFloat.size()
        << " (" << (100.0 * countExceeding / sinogramFloat.size()) << "%)" << std::endl;
}

unsigned char* float_to_byte(int X, int Y, float* in) {

    int p = X * Y; float min = in[p - 1], max = min;
    while (p--) {
        if (min > in[p])min = in[p];
        if (max < in[p]) max = in[p];
    }
    std::cout << " min = " << min << " max = " << max << std::endl;
    unsigned char* out = new unsigned char[X * Y];
    p = X * Y; while (p--)out[p] = (in[p] - min) / (max-min) * 255;
    return out;
}

