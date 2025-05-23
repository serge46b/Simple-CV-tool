#include <stdio.h>
#include <string.h>
#include "funcs.c"

#define MODE_AMOUNT 6
#define FROMAT_AMOUNT 3

enum mode
{
  MODE_MEDIAN,
  MODE_GAUSS,
  MODE_GAUSS_DECONV,
  MODE_EDGES,
  MODE_GRAYSCALE,
  MODE_EQUI
};

enum format
{
  FORMAT_PNG,
  FORMAT_JPG,
  FORMAT_BMP
};
const unsigned char *MODES[MODE_AMOUNT] = {"-median", "-gauss", "-gauss-deconv", "-edges", "-grayscale", "-equalize"};
const unsigned int MODE_ARGC[MODE_AMOUNT] = {1, 2, 2, 2, 0, 0};
const unsigned char *FORMATS[FROMAT_AMOUNT] = {"png", "jpg", "bmp"};

int main(int argc, char *argv[])
{
  if (argc < 1)
  {
    printf("Too few arguments provided. add -h to reveal help");
    return 10;
  }
  if (strcmp(argv[1], "-h") == 0)
  {
    printf("Image Processing Tool\n");
    printf("Usage: %s <input_image> <mode> [mode_args...] <output_image>.<format>\n\n", argv[0]);

    printf("Supported modes:\n");
    printf("  -median <kernel_size>          Apply median filter (noise reduction)\n");
    printf("  -gauss <kernel_size> <sigma>   Apply Gaussian blur\n");
    printf("  -gauss-deconv <k_size> <sigma> Apply Gaussian deconvolution (sharpening)\n");
    printf("  -edges <low_thresh> <high_thresh> Detect edges using Canny algorithm\n");
    printf("  -grayscale                    Convert to grayscale\n");
    printf("  -equalize                     Perform histogram equalization\n\n");

    printf("Supported output formats: png, jpg, bmp\n\n");

    printf("Examples:\n");
    printf("  %s input.jpg -median 5 output.png\n", argv[0]);
    printf("  %s input.png -gauss 7 1.5 output.jpg\n", argv[0]);
    printf("  %s input.bmp -edges 50 150 output.png\n", argv[0]);
    printf("  %s input.jpg -grayscale output.bmp\n", argv[0]);
    printf("  %s input.png -equalize output.png\n\n", argv[0]);

    printf("Note: Kernel sizes should be odd numbers (3, 5, 7, etc.)\n");
    printf("      Threshold values for edge detection typically range 0-255\n");
    return 0;
  }
  if (argc < 3)
  {
    printf("Please, specify mode and output image path. add -h to reveal help");
    return 10;
  }
  int prog_mode = MODE_AMOUNT;
  for (int i = 0; i < MODE_AMOUNT; i++)
    if (strcmp(argv[2], MODES[i]) == 0)
    {
      prog_mode = i;
      break;
    }
  if (prog_mode == MODE_AMOUNT)
  {
    printf("Unknown mode '%s', type -h for help", argv[2]);
    return 12;
  }
  if (argc < 4 + MODE_ARGC[prog_mode])
  {
    printf("Too few arguments for selected mode. Awaited %d arguments. add -h to reveal help", MODE_ARGC[prog_mode]);
    return 10;
  }
  int sv_format = FROMAT_AMOUNT;
  unsigned char *sv_filename = argv[3 + MODE_ARGC[prog_mode]];
  if (strlen(sv_filename) <= 3)
  {
    printf("Save path is incorrect '%s', type -h for help", sv_filename);
    return 13;
  }
  unsigned char *str_format = (unsigned char *)sv_filename + strlen(sv_filename) - 3;
  for (int i = 0; i < FROMAT_AMOUNT; i++)
  {
    if (strcmp(str_format, FORMATS[i]) == 0)
    {
      sv_format = i;
      break;
    }
  }
  if (sv_format == FROMAT_AMOUNT)
  {
    printf("Unknown format '%s', type -h for help", str_format);
    return 13;
  }

  printf("Loading image...\n");
  Image in_img = load_img(argv[1], STBI_default);
  if (in_img.img == NULL)
  {
    printf("Image load failed.\nReason: %s", stbi_failure_reason());
    return 1;
  }
  printf("Image data:\n\twidth: %d\n\theight: %d\n\tchannel_number: %d\n", in_img.x_size, in_img.y_size, in_img.channel_num);

  int ks = 0;
  double sigma = 0;
  Image out_img;
  switch (prog_mode)
  {
  case MODE_MEDIAN:
    ks = atoi(argv[3]);
    printf("Performing median filtering with kernel size = %d\n", ks);
    median_filter(&in_img, ks);
    break;
  case MODE_GAUSS:
    ks = atoi(argv[3]);
    sigma = atof(argv[4]);
    printf("Performing gauss bluring with kernel size = %d and sigma = %lf\n", ks, sigma);
    gauss_filter(&in_img, ks, sigma, 0);
    break;
  case MODE_GAUSS_DECONV:
    ks = atoi(argv[3]);
    sigma = atof(argv[4]);
    printf("Performing gauss deconvolution with kernel size = %d and sigma = %lf\n", ks, sigma);
    gauss_filter(&in_img, ks, sigma, 1);
    break;
  case MODE_EDGES:
    int low = atoi(argv[3]);
    int high = atof(argv[4]);
    printf("Detecting edges\n");
    out_img = edge_detection(in_img, low, high);
    free_img(in_img);
    in_img = out_img;
    break;
  case MODE_GRAYSCALE:
    printf("Converting image to grayscale\n");
    out_img = rgb2gray(in_img);
    free_img(in_img);
    in_img = out_img;
    break;
  case MODE_EQUI:
    printf("equalizing histogram\n");
    equalize_hist(&in_img);
    break;

  default:
    break;
  }
  int e_code = 0;
  switch (sv_format)
  {
  case FORMAT_PNG:
    e_code = save_png(in_img, sv_filename, 0);
    break;
  case FORMAT_JPG:
    e_code = save_jpg(in_img, sv_filename, 10);
    break;

  default:
    break;
  }
  if (e_code == 0)
  {
    printf("Image save failed");
    stbi_image_free(in_img.img);
    return 2;
  }
  printf("Image saved succesfully");
  stbi_image_free(in_img.img);
  return 0;
}