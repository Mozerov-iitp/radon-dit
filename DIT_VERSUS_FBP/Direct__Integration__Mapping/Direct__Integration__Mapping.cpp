//"from sinogram projection images to Cartesian frequency domain"
#define _CRT_SECURE_NO_WARNINGS
#include <string>     
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <conio.h>

#include "pfm_read_write.h"
#include "FBP.h"
#include "inverse_fourier_transform.h"
#include "Direct_integration.h"
#include "Direct_RadonT.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


void  CH3_to_RGB(const char* name, unsigned char* in, int X, int Y, int C) {


    int W = X; int H = Y;
    unsigned char* ret = new unsigned char[W * H * C];
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            int index = (y * W + x) * C;
            int c = 3;
            while (c--)
                ret[index + c] = in[x + y * W + c * W * H];

        }

    stbi_write_png(name,
        X, X, 3, ret, 0);
    delete[] ret;

}


void test_DIM(double lw, int func, const char* name, const char* path, int A) {
    // double lw - — spectrum correction parameter 
    // for incomplete data in the sinogram domain
    std::ofstream file(gen_name(RESULT_TXT, name, path, A));
   
    int X, X_in, Y, Ch;

   unsigned char* imageData = 
       stbi_load(gen_name(IMAGE_IN, name, path, A), &X_in, &Y, &Ch, 1);
   
    //------------------
    TST_RADON_IMG circle = { NULL,NULL};
    test_file_byte_sino(&circle, A, imageData, X_in,  Y, &X);
    write_pfm(gen_name(SINO, name, path, A), circle.sngm, circle.X, circle.A);
    stbi_write_png(gen_name(IMAGE_GT, name, path, A), circle.X, circle.X, 1, circle.img, circle.X);
    //End of generation
    bool Deb = 1; double rot_angle = M_PI / 4;
    if (Deb) {
        double psnr_img, psnr_sino;
        PSNR_rot(&psnr_img, &psnr_sino,
            circle.img, circle.X, A, rot_angle, circle.sngm);
        std::cout <<
            "PSNR of Back and Forth (M_PI / 4) rotation operator = "<< psnr_img << std::endl;

        std::cout << "PSNR of sinogram-to-sinogram of B-F rotation operator = " <<
            psnr_sino <<
            std::endl;
        if (all_images)
        {
            ROT_ERR.psnr_img += psnr_img;
            ROT_ERR.psnr_sng += psnr_sino;
        }
    }
    //____________________
    //Start of reconstruction
    int T, A_180;
    float* sino = read_pfm(gen_name(SINO, name, path, A), &T, &A_180);
    //////// DIM PROCESS ///////////////////////
    complex_dbl* f_uv = 
        dir_intgt_mapping_intrpl(lw, func,circle.sngm, T, A_180);//
    ////////////////////////////////////////////
    //////// INVERSE  FFT //////////////////////
    f_uv[0] = { get_mean_from_sino(circle.sngm,  T, A_180),0 };
    double* rec_img = fft_2D_inv(f_uv, T);
    ////////////////////////////////////////////
    unsigned char* byte_rec_image = new unsigned char[T * T];
    unsigned char* byte_rec_image_err = new unsigned char[T * T];
   truncate_double_to_byte(rec_img,byte_rec_image, circle.img, 
       byte_rec_image_err, T);
    stbi_write_png(gen_name(IMAGE_REC, name, path, A), T, T, 1, byte_rec_image, T);
    stbi_write_png(gen_name(IMAGE_ERR, name, path, A), T, T, 1, byte_rec_image_err, T);
    
    unsigned char* byte_rec_image_err_с = new unsigned char[3 * T * T];
    memset(byte_rec_image_err_с, 0, 3 * T * T);
   
        truncate_double_to_byte_с(rec_img, byte_rec_image,
            circle.img,
            byte_rec_image_err_с, T);
  
    CH3_to_RGB(gen_name(IMAGE_ERR_С, name, path, A), byte_rec_image_err_с, T, T, 3);
    
    ERRORS(file,
        rec_img, circle.img,
        forward_result(X, A, rec_img), circle.sngm,
        T * T, T * A);

    ///----------------------------------------------
    delete[] sino;
    delete[] circle.img;
    delete[] circle.sngm;
    delete[] byte_rec_image;
    delete[] byte_rec_image_err;
    delete[] rec_img;
    stbi_image_free(imageData);
    file.close();
}
void test_add_noise_DIM(int noise_flt, double noise_level_percent, 
    double lw, int func, const char* name, const char* path, int A) {
    // double lw (default = 4) — spectrum correction parameter 
    // for incomplete data in the sinogram domain
    if (noise_level_percent < 0.01)noise_flt = 0;
    std::ofstream file(gen_name(RESULT_TXT_NOISE, name, path, A));

    int X, X_in, Y, Ch;

    unsigned char* imageData =
        stbi_load(gen_name(IMAGE_IN, name, path, A), &X_in, &Y, &Ch, 1);

    //------------------
    TST_RADON_IMG circle = { NULL,NULL };
    test_file_byte_sino(&circle, A, imageData, X_in, Y, &X);
    write_pfm(gen_name(SINO, name, path, A), circle.sngm, circle.X, circle.A);
    stbi_write_png(gen_name(IMAGE_GT, name, path, A), circle.X, circle.X, 1, circle.img, circle.X);
    //End of generation 
    bool Deb = 1; double rot_angle = M_PI / 4;
    if (Deb) {
        double psnr_img, psnr_sino;
        PSNR_rot(&psnr_img, &psnr_sino,
            circle.img, circle.X, A, rot_angle, circle.sngm);
        std::cout <<
            "PSNR of Back and Forth (M_PI / 4) rotation operator = " << psnr_img << std::endl;

        std::cout << "PSNR of sinogram-to-sinogram of B-F rotation operator = " <<
            psnr_sino <<
            std::endl;
        if (all_images)
        {
            ROT_ERR.psnr_img += psnr_img;
            ROT_ERR.psnr_sng += psnr_sino;
        }
    }
    //____________________
    //Start of reconstruction
    int T, A_180;
    float* sino = read_pfm(gen_name(SINO, name, path, A), &T, &A_180);
    add_gaussian_noise(sino, T * A, noise_level_percent);
    // Filtr------------------------
    double sigx = noise_flt == 3 ? noise_level_percent * 25 :
        noise_level_percent * 50;
    double sigy = sigx; int w = 7;
    if (noise_flt == 1 || noise_flt == 3)
        gauss_flt(sino, sigx, sigy, T, A, w);
   
    //////// DIM PROCESS ///////////////////////
    complex_dbl* f_uv =
        dir_intgt_mapping_intrpl(lw, func, sino, T, A_180);//
    ////////////////////////////////////////////
    //////// INVERSE  FFT //////////////////////
    f_uv[0] = { get_mean_from_sino(circle.sngm,  T, A_180),0 };
    double* rec_img = fft_2D_inv(f_uv, T);
    ////////////////////////////////////////////

    if (noise_flt == 2 || noise_flt == 3)
        gauss_flt(rec_img, sigx, sigy, T, T, w);
    unsigned char* byte_rec_image = new unsigned char[T * T];
    unsigned char* byte_rec_image_err = new unsigned char[T * T];
    truncate_double_to_byte(rec_img, byte_rec_image, circle.img,
        byte_rec_image_err, T);
    stbi_write_png(gen_name(IMAGE_REC_NOISE, name, path, A), T, T, 1, byte_rec_image, T);
    stbi_write_png(gen_name(IMAGE_ERR_NOISE, name, path, A), T, T, 1, byte_rec_image_err, T);


    ERRORS(file,
        rec_img, circle.img,
        forward_result(X, A, rec_img), circle.sngm,
        T * T, T * A);
    ///----------------------------------------------
    delete[] sino;
    delete[] circle.img;
    delete[] circle.sngm;
    delete[] byte_rec_image;
    delete[] byte_rec_image_err;
    delete[] rec_img;
    stbi_image_free(imageData);
    file.close();
}
void test_FBP(int norm, int filter_type,  const char* name, const char* path, int A) {

    std::ofstream file(gen_name(RESULT_TXT, name, path, A));

    int X, X_in, Y, Ch;

    unsigned char* imageData =
        stbi_load(gen_name(IMAGE_IN, name, path, A), &X_in, &Y, &Ch, 1);

    //------------------
    TST_RADON_IMG circle = { NULL,NULL };
    test_file_byte_sino(&circle, A, imageData, X_in, Y, &X);
    write_pfm(gen_name(SINO, name, path, A), circle.sngm, circle.X, circle.A);
    stbi_write_png(gen_name(IMAGE_GT, name, path, A), circle.X, circle.X, 1, circle.img, circle.X);
    //End of generation 
    //____________________
    //Start of reconstruction
    int T, A_180;
    float* sino = read_pfm(gen_name(SINO, name, path, A), &T, &A_180);
    
    double* rec_img = fbp_reconstruct(sino, T, A, filter_type);
    if(norm)Norm_with_GT(rec_img, circle.img, T * T);
    ////////////////////////////////////////////
    unsigned char* byte_rec_image = new unsigned char[T * T];
    unsigned char* byte_rec_image_err = new unsigned char[T * T];
    truncate_double_to_byte(rec_img, byte_rec_image, circle.img,
        byte_rec_image_err, T);
    stbi_write_png(gen_name(IMAGE_REC, name, path, A), T, T, 1, byte_rec_image, T);
    stbi_write_png(gen_name(IMAGE_ERR, name, path, A), T, T, 1, byte_rec_image_err, T);
    unsigned char* byte_rec_image_err_с = new unsigned char[3 * T * T];

    truncate_double_to_byte_с(rec_img, byte_rec_image,
        circle.img,
        byte_rec_image_err_с, T);

    CH3_to_RGB(gen_name(IMAGE_ERR_С, name, path, A), byte_rec_image_err_с, T, T, 3);


    ERRORS(file,
        rec_img, circle.img,
        forward_result(X, A, rec_img), circle.sngm,
        T * T, T * A);
    ///----------------------------------------------
    delete[] sino;
    delete[] circle.img;
    delete[] circle.sngm;
    delete[] byte_rec_image;
    delete[] byte_rec_image_err;
    delete[] rec_img;
    stbi_image_free(imageData);
    file.close();
}
void test_add_noise_FBP(int noise_flt, double noise_level_percent, int norm, int filter_type,   const char* name, const char* path, int A) {

    std::ofstream file(gen_name(RESULT_TXT_NOISE, name, path, A));

    int X, X_in, Y, Ch;

    unsigned char* imageData =
        stbi_load(gen_name(IMAGE_IN, name, path, A), &X_in, &Y, &Ch, 1);

    //------------------
    TST_RADON_IMG circle = { NULL,NULL };
    test_file_byte_sino(&circle, A, imageData, X_in, Y, &X);
    write_pfm(gen_name(SINO, name, path, A), circle.sngm, circle.X, circle.A);
    stbi_write_png(gen_name(IMAGE_GT, name, path, A), circle.X, circle.X, 1, circle.img, circle.X);
    //End of generation 
    //____________________
    //Start of reconstruction
    int T, A_180;
    float* sino = read_pfm(gen_name(SINO, name, path, A), &T, &A_180);
    add_gaussian_noise(sino, T * A, noise_level_percent);
    // Filtr------------------------
    double sigx = noise_flt == 3 ? noise_level_percent * 25:
        noise_level_percent * 50;
    double sigy = sigx; int w = 7;
    if (noise_flt == 1 || noise_flt == 3)
        gauss_flt(sino, sigx, sigy, T, A, w);
  
    double* rec_img = fbp_reconstruct(sino, T, A, filter_type);
    if (norm)Norm_with_GT(rec_img, circle.img, T * T);
    ////////////////////////////////////////////
    if (noise_flt == 2 || noise_flt == 3)
        gauss_flt(rec_img, sigx, sigy, T, T, w);
    unsigned char* byte_rec_image = new unsigned char[T * T];
    unsigned char* byte_rec_image_err = new unsigned char[T * T];
    truncate_double_to_byte(rec_img, byte_rec_image, circle.img,
        byte_rec_image_err, T);
    stbi_write_png(gen_name(IMAGE_REC_NOISE, name, path, A), T, T, 1, byte_rec_image, T);
    stbi_write_png(gen_name(IMAGE_ERR_NOISE, name, path, A), T, T, 1, byte_rec_image_err, T);


    ERRORS(file,
        rec_img, circle.img,
        forward_result(X, A, rec_img), circle.sngm,
        T * T, T * A);
    ///----------------------------------------------
    delete[] sino;
    delete[] circle.img;
    delete[] circle.sngm;
    delete[] byte_rec_image;
    delete[] byte_rec_image_err;
    delete[] rec_img;
    stbi_image_free(imageData);
    file.close();
}

int exe_process_noise() {
    DIM_ERR = { 0,0,0,0};
    FBP_ERR = { 0,0,0,0 };
    FBP_NORM = { 0,0,0,0 };
    ROT_ERR = { 0,0,0,0 };
    all_images = 0;
    const char* Images[10] =
    { "Arnolfini", "butterfly", "Hunters",
         "lizard", "seashell", "still-life",
         "disk","Shepp-Logan",
         "Dance", "Time" };
    int disk = get_disk();
    int method = get_method();
    int A = get_angle_n();
    double noise = get_noise_n();
    int noise_flt = get_noise_flt();
    if (method != 3) {
        if (!method) {
            ///Direct integration
            int intr = get_interp();
            int img = get_image_n();
            if (img < 10)
            {
                const char* fld = gen_folder(disk, 0);
                test_add_noise_DIM(noise_flt, noise, 4, intr, Images[img], fld, A);
            }
            else
            {
                for (int i = 0; i < 10; i++)
                {
                    all_images = 1;
                    method_iml = 0;
                    const char* fld = gen_folder(disk, 0);
                    test_add_noise_DIM(noise_flt, noise, 4, intr, Images[i], fld, A);
                }
            }
        }
        else {
            if (method == 1) {
                ///FBP mean
                int flt = 0;
                int img = get_image_n();
                if (img < 10) {
                    const char* fld = gen_folder(disk, 1);
                    test_add_noise_FBP(noise_flt, noise, 0, flt, Images[img], fld, A);
                }
                else
                    for (int i = 0; i < 10; i++)
                    {
                        all_images = 1;
                        method_iml = 1;

                        const char* fld = gen_folder(disk, 1);
                        test_add_noise_FBP(noise_flt, noise, 0, flt, Images[i], fld, A);
                    }

            }
            else {
                ///FBP mean var
                int flt = 0;
                int img = get_image_n();
                if (img < 10) {
                    const char* fld = gen_folder(disk, 1);
                    test_add_noise_FBP(noise_flt, noise, 1, flt, Images[img], fld, A);
                }
                else
                    for (int i = 0; i < 10; i++)
                    {
                        all_images = 1;
                        method_iml = 2;
                        const char* fld = gen_folder(disk, 1);
                        test_add_noise_FBP(noise_flt, noise, 1, flt, Images[i], fld, A);
                    }

            }

        }

    }
    else {
        int intr = 1;
        int img = get_image_n();
        const char* fld[3] = { gen_folder(disk, 0),
                                gen_folder(disk, 1),
                                gen_folder(disk, 2) };
        int flt = 0;
        if (img < 10)
        {
            test_add_noise_DIM(noise_flt, noise, 4, intr, Images[img], fld[0], A);
            test_add_noise_FBP(noise_flt, noise, 0, flt, Images[img], fld[1], A);
            test_add_noise_FBP(noise_flt, noise, 1, flt, Images[img], fld[2], A);
        }

        else for (int i = 0; i < 10; i++)
        {
            all_images = 1;
            method_iml = 0;
            test_add_noise_DIM(noise_flt, noise, 4, intr, Images[i], fld[0], A);
            method_iml = 1;
            test_add_noise_FBP(noise_flt, noise, 0, flt, Images[i], fld[1], A);
            method_iml = 2;
            test_add_noise_FBP(noise_flt, noise, 1, flt, Images[i], fld[2], A);
        }
    }
    if (all_images)MEAN_ERRORS();
    waitEnter();
    return 0;
}

int exe_process() {
    DIM_ERR = { 0,0,0,0 };
    FBP_ERR = { 0,0,0,0 };
    FBP_NORM = { 0,0,0,0 };
    ROT_ERR = { 0,0,0,0 };
    all_images = 0;
    const char* Images[10] =
    { "Arnolfini", "butterfly", "Hunters",
         "lizard", "seashell", "still-life", 
         "disk","Shepp-Logan",
         "Dance", "Time" };
    int disk = get_disk();
    int method = get_method(); 
    int A = get_angle_n();
    if(method!=3){
        if (!method) {
            ///Direct integration
            int intr = get_interp();
            int img = get_image_n();
            if (img < 10)
            {
                const char* fld = gen_folder(disk, 0);
                test_DIM(4, intr, Images[img], fld, A);
            }
            else
            {
                for (int i = 0; i < 10; i++)
                {
                    all_images = 1;
                    method_iml = 0;
                    const char* fld = gen_folder(disk, 0);
                    test_DIM(4, intr, Images[i], fld, A);
                }
            }
        }
        else {
            if (method == 1) {
                ///FBP mean
                int flt = get_filter();
                int img = get_image_n();
                if (img < 10){
                    const char* fld = gen_folder(disk, 1);
                    test_FBP(0, flt, Images[img], fld, A);
                }
                else
                    for (int i = 0; i < 10; i++)
                {
                        all_images = 1;
                        method_iml = 1;
                    const char* fld = gen_folder(disk, 1);
                    test_FBP(0, flt, Images[i], fld, A);
                }

            }
            else {
                ///FBP mean var
                int flt = get_filter();
                int img = get_image_n();
                if (img < 10) {
                    const char* fld = gen_folder(disk, 1);
                    test_FBP(1, flt, Images[img], fld, A);
                }
                else
                    for (int i = 0; i < 10; i++)
                    {
                        all_images = 1;
                        method_iml = 2;
                        const char* fld = gen_folder(disk, 1);
                        test_FBP(1, flt, Images[i], fld, A);
                    }

            }

        }
    
    }
    else{
        int intr = 2;
        int img = get_image_n();
        const char* fld[3] = { gen_folder(disk, 0),
                                gen_folder(disk, 1), 
                                gen_folder(disk, 2)};
        int flt = 1;
        if (img < 10)
        {
            test_DIM(4, intr, Images[img], fld[0], A);
            test_FBP(0, flt, Images[img], fld[1], A);
            test_FBP(1, flt, Images[img], fld[2], A);
        }

        else for (int i = 0; i < 10; i++)
        {
            all_images = 1;
            method_iml = 0;
            test_DIM(4, intr, Images[i], fld[0], A);
            method_iml = 1;
            test_FBP(0, flt, Images[i], fld[1], A);
            method_iml = 2;
            test_FBP(1, flt, Images[i], fld[2], A);
        }
    }
    if(all_images)MEAN_ERRORS();
    waitEnter();
    return 0;
}


int main()
{   
    std::cout << "    This is a demonstration program that models  the Radon  transform" << std::endl;
    std::cout << "and reconstructs  the  original  image   from  a  simulated  sinogram." << std::endl;
    std::cout << "    Reconstruction can be performed using both the proposed method  -" << std::endl;
    std::cout << "Direct Integration and Filtered Back Projection (FBP) for  comparison." << std::endl;

    std::cout << "    The main algorithm does not assume the presence of noise. However,  " << std::endl;
    std::cout << "a version with simulated additive Gaussian noise and its  suppression " << std::endl;
    std::cout << "is  also  implemented." << std::endl;

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "Which version do you choose?" << std::endl;
    std::cout << std::endl;
    std::cout << "Main version (noise-free)              - 0" << std::endl;
    std::cout << "Version with Gaussian noise simulation - 1" << std::endl;
    std::cout << std::endl;
    int noise; std::cin >> noise;
    std::cout << std::endl;
    if (!noise)exe_process();
    else exe_process_noise();
   
}

