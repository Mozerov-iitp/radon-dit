
// ============================================================================
// Direct_integration.h
// ============================================================================
// Copyright (c) Mikhail G. Mozerov
// 
// This file is part of the implementation accompanying the paper:
// 
//   Title:  "The Direct Integration Theorem: A Rigorous Framework for Consistent
//           Discrete Solutions of the Inverse Radon Problem"
//   Author: Mikhail G. Mozerov
//   arXiv:  https://arxiv.org/abs/2605.09020 (submitted to IEEE TPAMI)
//   Code:   https://github.com/Mozerov-iitp/radon-dit/
// 
// Permission is hereby granted to use this code for academic and research
// purposes, provided that the original paper is cited.
// ============================================================================
//
// This header implements the core of the Direct Integration Method (DIT)
// for tomographic image reconstruction. The key function
// dir_intgt_mapping_intrpl() directly maps sinogram data into the 2D Fourier
// domain using angular interpolation and a frequency‑domain summation.
// 
// For details, please refer to the above paper, Section III (Direct Integration
// Theorem) and Section IV (Discrete Implementation).
//
// ============================================================================

#pragma once
#include "inverse_fourier_transform.h"
#include "Generate_sino.h"


/**
 * Extracts a projection at angle 'tet' using nearest‑neighbour interpolation.
 * @param in   Sinogram array (angles major, detectors minor).
 * @param D    Number of detectors (image width).
 * @param A    Number of projection angles over 180°.
 * @param tet  Desired angle in radians.
 * @return     Complex array (real part = projection values, imag = 0).
 */

complex_dbl* GetSngmTet0(float* in, int D, int A, double tet) {

	complex_dbl* ret = new complex_dbl[D];
	double da = M_PI / A;
	double a = tet / da;
	int a1 = (int)a;
	int a2 = (a1 + 1) % A;
	double d = a - a1;
	int i_tet = d < 0.5 ? a1 : a2;

	int p = D; while (p--)ret[p] = { (double)in[p + i_tet * D], 0 };

	return ret;
}
/**
 * Linear interpolation between two consecutive projections.
 */
complex_dbl* GetSngmTet1(float* in, int D, int A, double tet) {
	// Allocate output buffer for the interpolated profile
	complex_dbl* ret = new complex_dbl[D];

	// Angular step size for the [0, PI) range
	double da = M_PI / A;

	// 1. Handle negative angles from atan2 [-PI, PI] by shifting them to the positive domain
	bool global_flip = false;
	if (tet < 0.0) {
		tet += M_PI;
		global_flip = true; // Signal that the final spatial profile must be mirrored
	}
	if (tet >= M_PI) {
		tet -= M_PI;
		global_flip = !global_flip;
	}

	// Map the continuous angle to fractional coordinate space
	double a = tet / da;
	int a1 = (int)a;
	double d = a - a1; // Linear interpolation weight

	// 2. Setup the right-hand neighbor angle index and its boundary condition
	int a2 = a1 + 1;
	bool local_flip = false;

	// If the right-hand point hits or exceeds 180 degrees, warp it back to 0 with an axis flip
	if (a2 >= A) {
		a2 -= A;
		local_flip = true;
	}

	// 3. Process the interpolation directly without temporary heap allocations
	for (int p = 0; p < D; p++) {
		// Retrieve values from the left angular neighbor (always inside bounds or normalized)
		double val1 = (double)in[a1 * D + p];

		// Retrieve values from the right angular neighbor with a conditional channel flip
		int target_p = local_flip ? (D - 1 - p) : p;
		double val2 = (double)in[a2 * D + target_p];

		// Perform linear blend
		double vl = val1 * (1.0 - d) + val2 * d;

		// Write to the output array factoring in the global negative angle orientation
		int final_p = global_flip ? (D - 1 - p) : p;
		ret[final_p] = complex_dbl(vl, 0.0);
	}

	return ret;
}
//complex_dbl* GetSngmTet1(float* in, int D, int A, double tet) {
//
//	double* b[2] = { new double [D], new double[D] };
//	complex_dbl* ret = new complex_dbl[D];
//	double da = M_PI / A;
//	double a = tet / da;
//	int a1 = (int)a;
//	int a2 = (a1 + 1) > A - 1 ? A - 1 : a1 + 1;
//	double d = a - a1;
//
//	int p = D; while (p--)
//	{
//		b[0][p] = (double)in[p + a1 * D];
//		b[1][p] = (double)in[p + a2 * D];
//
//	}
//	p = D; while (p--)ret[p] = 
//	{ (double)in[p + a1 * D]*(1.-d) + (double)in[p + a2 * D] * (d), 0};
//	delete[] b[0];
//	delete[] b[1];
//	return ret;
//}
/**
 * Cubic spline interpolation (Catmull‑Rom style) over four adjacent angles.
 */
//complex_dbl* GetSngmTetCubic(float* in, int D, int A, double tet) {
//
//	double* b[4] = { new double[D], new double[D], new double[D], new double[D] };
//	complex_dbl* ret = new complex_dbl[D];
//	double da = M_PI / A;
//	double a = tet / da;
//	int a1 = (int)a;
//	int a_1 = a1 -1 <0 ? 0: a1-1;
//	int a2 = (a1 + 1) > A - 1 ? A - 1 : a1 + 1;
//	int a3 = (a1 + 2) > A - 1 ? A - 1 : a1 + 2;
//	double d = a - a1;
//
//	int p = D; while (p--)
//	{
//		b[0][p] = (double)in[p + a_1 * D];
//		b[1][p] = (double)in[p + a1 * D];
//		b[2][p] = (double)in[p + a2 * D];
//		b[3][p] = (double)in[p + a3 * D];
//
//	}
//	p = D; while (p--) {
//		double bb[4] = { b[0][p],b[1][p], b[2][p], b[2][p] };
//		double vl = cubicInterpolate(bb, d);
//		ret[p] = { vl,0 };
//	}
//	delete[] b[0];
//	delete[] b[1];
//	delete[] b[2];
//	delete[] b[3];
//	return ret;
//}
///**
// * Dispatch function for angular interpolation.
// * @param func  0=nearest, 1=linear, 2=cubic.
// */

complex_dbl* GetSngmTetCubic(float* in, int D, int A, double tet) {
	// Dynamically allocate the resulting complex projection vector
	complex_dbl* ret = new complex_dbl[D];

	// Angular step size of the sparse input sinogram (spanning [0, PI))
	double da = M_PI / A;

	// 1. Normalize negative angles from atan2 [-PI, PI] to the positive [0, 2*PI] domain
	bool global_flip = false;
	if (tet < 0.0) {
		tet += M_PI;
		global_flip = true; // Mark that the detector channel axis must be reversed
	}
	// Handle edge case where tet could slightly exceed or equal PI due to precision
	if (tet >= M_PI) {
		tet -= M_PI;
		global_flip = !global_flip;
	}

	// Convert continuous angle to fractional coordinate index
	double a = tet / da;
	int a1 = (int)a;
	double d = a - a1; // Interpolation weight (fractional part)

	// Temporary arrays to hold data for the 4 neighboring angles
	double* b[4] = { new double[D], new double[D], new double[D], new double[D] };

	// 2. Compute indices for the 4-point cubic stencil with 180-degree Radon boundary conditions
	for (int i = -1; i <= 2; i++) {
		int neighbor_angle = a1 + i;
		bool local_flip = false;

		// Apply anti-periodic wrap-around rules if index steps outside [0, A-1]
		if (neighbor_angle < 0) {
			neighbor_angle = A + neighbor_angle;
			local_flip = true;
		}
		else if (neighbor_angle >= A) {
			neighbor_angle = neighbor_angle - A;
			local_flip = true;
		}

		// Fill the current angle row, account for geometric flip conditions
		for (int p = 0; p < D; p++) {
			// Determine target channel tracking based on the anti-periodic layout
			int target_p = local_flip ? (D - 1 - p) : p;

			// Your sinogram indexing convention layout: [angle * D + offset]
			b[i + 1][p] = (double)in[neighbor_angle * D + target_p];
		}
	}

	// 3. Perform the cubic interpolation and handle global coordinate sign inversion
	for (int p = 0; p < D; p++) {
		// FIXED: Replaced the duplicated b[2][p] with b[3][p]
		double bb[4] = { b[0][p], b[1][p], b[2][p], b[3][p] };
		double vl = cubicInterpolate(bb, d);

		// Map the final value to its appropriate channel based on the initial angle sign
		int final_p = global_flip ? (D - 1 - p) : p;
		ret[final_p] = complex_dbl(vl, 0.0);
	}

	// Clean up temporary memory buffers
	delete[] b[0]; delete[] b[1]; delete[] b[2]; delete[] b[3];

	return ret;
}

complex_dbl* GetSngmTetIntrpl(int func, 
	float* in, int D, int A, double tet) {

	complex_dbl* ret = NULL;
	if (func == 2)return GetSngmTetCubic(in, D, A, tet);
	if (func == 1)return GetSngmTet1(in, D, A, tet);
	return GetSngmTet0(in, D, A, tet);

	return ret;
}
/**
 * Computes the contribution of all frequencies along the line t (detector index)
 * for a given Fourier point (u,v). Used in Direct Integration.
 */
complex_dbl sum_for_point_u_v(complex_dbl* in, int u, int v, int D) {


	double  omg = 2. * M_PI * sqrt(u * u + v * v) / D;
	complex_dbl sum = { 0,0 };
	for (int t = -D / 2 + 1; t < D / 2; t++) {

		complex_dbl omg_t = { 0, -omg * (t) };
		sum += in[t + D / 2] * exp(omg_t);
	}
	return sum;
}
/**
 * Applies a circular low‑pass filter (energy compensation) to the Fourier data.
 * Keeps frequencies inside radius D/2, rescales to preserve total energy.
 */
void kill_frq(complex_dbl* to_kill, int D, int A) {


	double da = M_PI / A; 
	
	double org_frq = 0.0; 
	for (int p = 1; p < D * D; p++)org_frq += std::norm(to_kill[p]);
	for (int i_v = D / 2 - 1; i_v > -D / 2; i_v--)
	{
		for (int i_u = D / 2 - 1; i_u > -D / 2; i_u--) {
			double r = sqrt(i_u * i_u + i_v * i_v);
			int u = i_u < 0 ? i_u + D : i_u;
			int v = i_v < 0 ? i_v + D : i_v;
			int p = u + v * D;
			bool do_it =  r<D/2;
			
			if (!do_it)to_kill[p] = 0;
		
			}

		

	}
	double f_frq = 0.0;
	for (int p = 1; p < D * D; p++)f_frq += std::norm(to_kill[p]);
	double alph = sqrt(org_frq / f_frq);

	for (int p = 1; p < D * D; p++)to_kill[p] *= alph;

}
/**
 * Overloaded version with additional rectangular cut‑off (lw parameter).
 * lw > 0 restricts to |u|<lw*A/π and |v|<lw*A/π before circular cut.
 */
void kill_frq(double lw, complex_dbl* to_kill, int D, int A) {

	if (lw)
	{
		double da = M_PI / A;
		double sq = lw * A / M_PI;
		double R = D / 2 < sq ? D / 2 : sq;
		double org_frq = 0.0;
		for (int p = 1; p < D * D; p++)org_frq += std::norm(to_kill[p]);
		for (int i_v = D / 2 - 1; i_v > -D / 2; i_v--)
		{
			for (int i_u = D / 2 - 1; i_u > -D / 2; i_u--) {
				double ru = abs(i_u);
				double rv = abs(i_v);

				int u = i_u < 0 ? i_u + D : i_u;
				int v = i_v < 0 ? i_v + D : i_v;
				int p = u + v * D;
				bool do_it = ru < R && rv < R;

				if (!do_it)to_kill[p] = 0;

			}



		}
		double f_frq = 0.0;
		for (int p = 1; p < D * D; p++)f_frq += std::norm(to_kill[p]);
		double alph = sqrt(org_frq / f_frq);
		
		for (int p = 1; p < D * D; p++)to_kill[p] *= alph;
		kill_frq(to_kill, D, A);
	}
	 else kill_frq(to_kill, D, A);
}
// ============================================================================
// Function: dir_intgt_mapping_intrpl
// ============================================================================
// Implements the Direct Integration Mapping (DIM) – the core reconstruction
// algorithm from the paper:
//   "The Direct Integration Theorem: A Rigorous Framework for Consistent
//    Discrete Solutions of the Inverse Radon Problem"
//   Mikhail G. Mozerov (arXiv, submitted to IEEE TPAMI)
//
// Brief description:
//   For each point (u,v) in the 2D Fourier domain (zero‑centered), the function
//   determines the corresponding projection angle theta = atan2(v,u). It then
//   extracts the 1D sinogram projection at that angle (using angular
//   interpolation: nearest, linear, or cubic) and computes a weighted sum over
//   detector positions t (via sum_for_point_u_v). The result fills the Fourier
//   cell (u,v). Conjugate symmetry is exploited to reduce computation.
//
// Parameters:
//   lw    – Spectrum correction parameter (frequency cut‑off / energy rescaling)
//   func  – Interpolation type: 0 = nearest neighbour, 1 = linear, 2 = cubic
//   in    – Input sinogram (float array, size A x D)
//   D     – Number of detectors (and output image width)
//   A     – Number of projection angles (over 180°)
//
// Returns:
//   Pointer to a newly allocated array of complex_dbl of size D×D,
//   containing the Fourier domain representation of the reconstructed image.
//   The caller is responsible for deleting it.
// ============================================================================
complex_dbl* dir_intgt_mapping_intrpl(double lw,  int func, float* in, int D, int A) {

	std::cout << "=======  Direct Integration Mapping  ======" << std::endl;
	std::cout << std::endl;
	clock_t start = clock();
	complex_dbl* out = new complex_dbl[D * D];
	memset(out, 0, sizeof(complex_dbl) * D * D);
	complex_dbl* bk = NULL;

	int mod_T = D / 42;
	for (int i_v = D / 2 - 1; i_v > -D / 2; i_v--)
	{
		for (int i_u = D / 2 - 1; i_u > -D / 2; i_u--) {

			int u = i_u < 0 ? i_u + D : i_u;
			int v = i_v < 0 ? i_v + D : i_v;
			int p = u + v * D;

			int uu = -i_u < 0 ? -i_u + D : -i_u;
			int vv = -i_v < 0 ? -i_v + D : -i_v;
			int pp = uu + vv * D;

			bool do_it = (i_v > 0) || ((i_u > 0) && (i_v == 0));
			//--------------------
			double tet = atan2(i_v, i_u);




			if (do_it) {
				bk = GetSngmTetIntrpl(func, in, D, A, tet);
				out[p] = sum_for_point_u_v(bk, i_u, i_v, D);
				delete[] bk;
			}
			else out[p] = { out[pp].real(), -out[pp].imag() };




		}
		if ((i_v + D / 2) % mod_T == 0) { std::cout << char(219); }

	}
	kill_frq(lw, out, D, A);
	//--------------------------
	std::cout << std::endl;
	std::cout << std::endl;
	clock_t finish = clock();
	float time = (float)(((double)(finish - start)) / CLOCKS_PER_SEC);
	std::cout << "================  DIT End  ================" << std::endl;
	std::cout << std::endl;
	std::cout << "        Execution time: " << time << " seconds" << std::endl;

	std::cout << std::endl;
	std::cout << "===========================================" << std::endl;

	return out;
}

