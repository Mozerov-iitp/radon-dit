#pragma once
#include <fstream>
#include <string>
#include <vector>
bool load_pfm(const std::string& filename, std::vector<float>& data,
    int& width, int& height) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    // Read header
    std::string magic;
    float scale;
    file >> magic;          // "Pf" for grayscale or "PF" for RGB
    file >> width >> height;
    file >> scale;

    // Skip one whitespace byte (newline or space) after header
    file.ignore(1);

    // Check that it's a single-channel (grayscale) PFM file
    if (magic != "Pf") {
        return false; // Not a single-channel PFM
    }

    // Prepare vector for data
    size_t data_size = width * height;
    data.resize(data_size);

    // Read data (stored bottom-up)
    for (int y = height - 1; y >= 0; --y) {
        file.read(reinterpret_cast<char*>(&data[y * width]), width * sizeof(float));
    }

    // If byte order in file (scale > 0) doesn't match system,
    // additional conversion is needed (scale = -1.0 in our example).
    if (scale > 0.0f) {
        // Swap bytes for big-endian conversion (rarely needed on Windows)
        // Implementation would go here
    }

    return !file.fail();
}

/**
 * Reads a PFM file and returns a float array (caller must delete[])
 * @param filename Path to the PFM file
 * @param out_w Output parameter for width (can be nullptr)
 * @param out_h Output parameter for height (can be nullptr)
 * @return Pointer to float array on success, nullptr on error
 * @warning Caller is responsible for deleting[] the returned pointer
 */
float* read_pfm(const std::string& filename, int* out_w, int* out_h) {
    // Open file directly, avoiding intermediate vector
    std::ifstream file(filename, std::ios::binary);
    if (!file) return nullptr;

    std::string magic;
    int width, height;
    float scale;

    // Read header
    file >> magic >> width >> height >> scale;
    if (magic != "Pf") return nullptr;

    // Skip whitespace
    file.ignore(1);

    // Return dimensions through output parameters
    if (out_w) *out_w = width;
    if (out_h) *out_h = height;

    // Allocate final array immediately
    float* ret = new float[width * height];

    // Read directly into it (bottom-up)
    for (int y = height - 1; y >= 0; --y) {
        file.read(reinterpret_cast<char*>(&ret[y * width]),
            width * sizeof(float));
    }

    // Handle endianness if scale > 0
    if (scale > 0.0f) {
        // Big-endian to little-endian conversion would go here
        // This is rarely needed on x86/x64 systems
    }

    return file.fail() ? (delete[] ret, nullptr) : ret;
}

/**
 * Writes a float array to a PFM (Portable Float Map) file
 * @param filename Output file path
 * @param out Pointer to float data array (size must be w * h)
 * @param w Image width (must be > 0)
 * @param h Image height (must be > 0)
 * @return true on success, false on error
 */
bool write_pfm(const std::string& filename, const float* out, int w, int h) {
    // Validate input parameters
    if (!out || w <= 0 || h <= 0) {
        std::cerr << "Error: Invalid parameters for PFM write" << std::endl;
        return false;
    }

    // Open file in binary mode
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return false;
    }

    try {
        // 1. Write PFM header
        // Magic number "Pf" for single-channel (grayscale) image
        file << "Pf\n";

        // Width and height
        file << w << " " << h << "\n";

        // Scale factor (negative for little-endian)
        // -1.0 indicates little-endian byte order (standard for x86/x64)
        file << "-1.0\n";

        // 2. Write image data
        // PFM stores data BOTTOM-UP
        for (int y = h - 1; y >= 0; --y) {
            const float* row_ptr = out + (y * w);
            file.write(reinterpret_cast<const char*>(row_ptr), w * sizeof(float));
        }

        // Check if any write error occurred
        if (!file.good()) {
            std::cerr << "Error writing data to file" << std::endl;
            return false;
        }
        std::cout << std::endl;
        std::cout << "Radon image successfully written: " << filename
            << " (" << w << "x" << h << ", "
            << w * h * sizeof(float) << " bytes)" << std::endl;
        std::cout << std::endl;
        return true;

    }
    catch (const std::exception& e) {
        std::cerr << "Exception while writing PFM: " << e.what() << std::endl;
        return false;
    }
}