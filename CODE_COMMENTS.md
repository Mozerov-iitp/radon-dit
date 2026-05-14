# Code Comments Summary – Radon-DIT

This document provides a high‑level overview of the C++ implementation of the **Direct Integration Theorem (DIT)** for tomographic image reconstruction.  
For line‑by‑line comments, please refer to the source files (`.h` and `.cpp`).  

> **Corresponding paper:**  
> M. G. Mozerov, *"The Direct Integration Theorem: A Rigorous Framework for Consistent Discrete Solutions of the Inverse Radon Problem"*,  
> arXiv preprint, submitted to IEEE TPAMI (May 2026).  
> Code & data: [https://github.com/Mozerov-iitp/radon-dit/](https://github.com/Mozerov-iitp/radon-dit/)

---

## Project Structure

| File | Description |
|------|-------------|
| `Direct__Integration__Mapping.cpp` | Main program loop, console menu, test drivers (with/without Gaussian noise). |
| `Direct_integration.h` | **Core Direct Integration Method (DIM)** – Fourier domain mapping, angular interpolation, frequency cut‑off. |
| `Direct_RadonT.h` | Radon transform, quality metrics (PSNR, SSIM, MAE), image I/O helpers, noise simulation, error visualisation. |
| `FBP.h` | Filtered Back Projection (ramp / Shepp‑Logan filter), back‑projection loop, cubic interpolation. |
| `inverse_fourier_transform.h` | 1D/2D FFT/IFFT (Cooley–Tukey), bit‑reversal. |
| `Generate_sino.h` | Discrete Radon transform with bicubic interpolation (double precision internally). |
| `pfm_read_write.h`, `stb_image.h`, `stb_image_write.h` | Third‑party image I/O (PFM, PNG). |

---

## Key Data Structures

- `TST_RADON_IMG`  
  Holds sinogram (`float* sngm`), ground truth image (`unsigned char* img`), and dimensions (`int X` = image size, `int A` = number of projection angles). Used as a container for test images.

- `TYPE_ERRORS`  
  Aggregates quality metrics for batch processing:  
  `psnr_img`, `psnr_sng`, `ssim_img`, `var` (ratio of standard deviations).

- Global accumulators: `DIM_ERR`, `FBP_ERR`, `FBP_NORM`, `ROT_ERR`.

---

## Main Functions

### Direct Integration (DIM) – `Direct_integration.h`

- **`dir_intgt_mapping_intrpl(lw, func, in, D, A)`**  
  Core DIM reconstruction. Maps sinogram → Fourier domain.  
  - Iterates over all Fourier coordinates `(u,v)` with zero‑centered indexing.  
  - Computes angle `θ = atan2(v,u)`.  
  - Extracts projection at `θ` using angular interpolation (`func`: 0=nearest, 1=linear, 2=cubic) – functions `GetSngmTet0/1/Cubic`.  
  - Calls `sum_for_point_u_v` to compute the 1D sum over detector index `t`.  
  - Exploits conjugate symmetry for negative frequencies (`out[p] = conj(out[pp])`).  
  - Applies frequency cut‑off and energy rescaling via `kill_frq`.  
  - Returns complex Fourier array (caller must free).

- **`GetSngmTet0/1/Cubic`** – angular interpolation between consecutive projection angles (nearest, linear, cubic Catmull‑Rom).

- **`sum_for_point_u_v(in, u, v, D)`**  
  Computes `Σ_{t} in[t] * exp(-j·ω·t)` with `ω = 2π√(u²+v²)/D`. The input `in` is a 1D projection (complex, but imag part = 0).

- **`kill_frq`** (two overloads)  
  - First overload: circular low‑pass filter (radius `D/2`), then rescales total energy to match original.  
  - Second overload: if `lw > 0`, first applies rectangular cut‑off `|u|,|v| < lw·A/π` before the circular filter, then calls the first overload.

### Filtered Back Projection – `FBP.h`

- **`fbp_reconstruct(sinogram, detector_size, num_angles, filter_type)`**  
  1. Apply 1D FFT to each projection (row).  
  2. Multiply by filter: ramp (`filter_type = 0`) or Shepp‑Logan (`filter_type = 1`).  
  3. Inverse FFT to obtain filtered sinogram.  
  4. For each pixel, back‑project using `reconstruct_pixel` (cubic interpolation in detector coordinate).  
  5. Normalise by number of angles and adjust mean.  
  Returns double image (caller frees).

- **`reconstruct_pixel(filtered_sinogram, pixel_idx, image_size, detector_size, num_angles)`**  
  For a single pixel: loop over all angles, compute detector coordinate `x·cosθ + y·sinθ`, interpolate value in filtered projection, sum.

- **`cubic_interpolate(t, f)`** – Catmull‑Rom spline for 1D interpolation (4 points).

### Radon Transform – `Generate_sino.h`

- **`discreteRadonTransformDoublePtr(image, imageSize, numAngles, numOffsets)`**  
  Forward Radon transform:  
  - For each angle `θ` and offset `s`, integrate along the line `x·cosθ + y·sinθ = s` using **bicubic interpolation** (double precision internally).  
  - Output sinogram as `float*` (caller frees).  
  - Overloaded for input `vector<float>` and `vector<double>`.

- **`bicubicInterpolateDouble(image, width, height, x, y)`**  
  Bicubic interpolation using Catmull‑Rom kernel. Handles boundaries by mirroring.

### Quality Metrics – `Direct_RadonT.h`

- **`get_PSNR`** – multiple overloads for `double* vs unsigned char*`, `unsigned char* vs unsigned char*`, `double* vs float*`, etc. Includes circular mask (only pixels inside inscribed circle).

- **`SSIM`** – full‑reference Structural Similarity Index (constants from original paper).

- **`MAE`** – Mean Absolute Error in percent (relative to 255).

- **`ERRORS`** – writes all metrics to a text file and updates global accumulators (`DIM_ERR`, `FBP_ERR`, `FBP_NORM`).

- **`MEAN_ERRORS`** – prints average metrics after batch processing (e.g., over all 10 test images).

### Noise Simulation & Filtering

- **`add_gaussian_noise(image, size, noise_level_percent)`**  
  Noise standard deviation = `mean_value * noise_level_percent`. Clamped to non‑negative.

- **`gauss_flt(in, sigx, sigy, X, Y, w)`** – 2D Gaussian convolution (separable, mirror boundary). Works for `double*` and `float*`.

- **`add_poisson_noise`, `add_shift_noise`** – additional noise models (available but not used in the main paper).

### Utility & Visualisation

- **`truncate_double_to_byte`** – two overloads:  
  1. `(double* in, unsigned char* result, int T)` – min‑max scaling to [0,255].  
  2. `(double* in, unsigned char* result, unsigned char* gt, unsigned char* err_m, int T)` – also computes absolute error map: `err = |gt - rec| * 10`, clamped to 255.

- **`truncate_double_to_byte_color_error(in, result, gt, err_color, T)`**  
  Produces an **RGB error map**:  
  - Red channel = positive error (rec > gt) amplified ×10.  
  - Blue channel = negative error (gt > rec) amplified ×10.  
  - Green channel = max(red, blue).  
  Used for visualisation of reconstruction errors.

- **`CH3_to_RGB(name, in, X, Y, C)`** – converts planar CH3 buffer to interleaved RGB and saves as PNG.

- **`rotate` / `bi_cube`** – image rotation with bicubic interpolation (used for consistency check of the Radon transform).

- **`mediana_window` / `mediana`** – median filter (not used in main results, but available).

### Console Menu & User Input – `Direct_RadonT.h`

- `get_disk()`, `get_angle_n()`, `get_method()`, `get_noise_n()`, `get_noise_flt()`, `get_interp()`, `get_filter()`, `get_image_n()` – interactive prompts.

- `gen_folder(disk, method)` – builds output directory path based on user choices.

- `waitEnter()` – portable “Press Enter to continue”.

---

## Important Implementation Notes

1. **Image padding** – all input images are zero‑padded to the next power‑of‑two size and centred. Only pixels inside the inscribed circle are used for metrics.

2. **Fourier domain symmetry** – DIM computes only half of the Fourier space (upper half‑plane + positive u on the axis) and uses conjugate symmetry for the rest.

3. **Angular interpolation** – because projections are sampled at discrete angles over 180°, the DIT requires values at arbitrary angles. The code provides nearest, linear, and cubic (Catmull‑Rom) interpolation.

4. **Frequency cut‑off** – `kill_frq` ensures that only frequencies up to the Nyquist limit are kept; energy rescaling prevents loss of contrast.

5. **Noise simulation** – Gaussian noise level is defined as a percentage of the **mean sinogram value**, which is physically meaningful for attenuation coefficients.

6. **Memory management** – all functions that return `new` allocated memory (e.g., `dir_intgt_mapping_intrpl`, `fbp_reconstruct`, `discreteRadonTransformDoublePtr`) require the caller to `delete[]`. The main test functions do this correctly.

---

## Batch Processing & Reporting

When the user selects **“All images”** (`get_image_n() == 10`), the program runs on all 10 test images, accumulates errors in the global structures, and finally prints average metrics via `MEAN_ERRORS()`. Results are also saved per image in `.txt` files.

---

## Citation

If you use this code in your research, please cite the corresponding paper:

```bibtex
@article{Mozerov2026DIT,
  title   = {The Direct Integration Theorem: A Rigorous Framework for Consistent Discrete Solutions of the Inverse Radon Problem},
  author  = {Mikhail G. Mozerov},
  journal = {arXiv preprint},
  year    = {2026},
  note    = {Submitted to IEEE Transactions on Pattern Analysis and Machine Intelligence},
  url     = {https://github.com/Mozerov-iitp/radon-dit/}
}