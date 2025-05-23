#include "funcs.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#define STBI_WINDOWS_UTF8
#define STBI_FAILURE_USERMSG
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"
#define STBIW_WINDOWS_UTF8
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"

Image load_img(unsigned char *path, int desired_channel_num)
{
  int x, y, n;
  unsigned char *img = stbi_load(path, &x, &y, &n, desired_channel_num);
  Image in_img;
  in_img.img = img;
  in_img.x_size = x;
  in_img.y_size = y;
  in_img.channel_num = n;
  return in_img;
}

int save_png(Image img, unsigned char *path, int stride_bytes)
{
  return stbi_write_png(path, img.x_size, img.y_size, img.channel_num, img.img, stride_bytes);
}

int save_bmp(Image img, unsigned char *path)
{
  return stbi_write_bmp(path, img.x_size, img.y_size, img.channel_num, img.img);
}

int save_jpg(Image img, unsigned char *path, int quality)
{
  return stbi_write_jpg(path, img.x_size, img.y_size, img.channel_num, img.img, quality);
}

unsigned char get_pixel_grayscale(Image img, int x, int y)
{
  assert(img.channel_num <= 1);
  assert(img.x_size > x);
  assert(img.y_size > y);
  assert(x >= 0);
  assert(y >= 0);
  return img.img[y * img.x_size + x];
}

unsigned char get_pixel(Image img, int x, int y, int c)
{
  assert(img.channel_num > c);
  assert(img.x_size > x);
  assert(img.y_size > y);
  assert(c >= 0);
  assert(x >= 0);
  assert(y >= 0);
  return img.img[y * img.x_size * img.channel_num + x * img.channel_num + c];
}

void set_pixel_grayscale(Image *modifiable_img, int x, int y, unsigned char value)
{
  assert(modifiable_img->channel_num <= 1);
  assert(modifiable_img->x_size > x);
  assert(modifiable_img->y_size > y);
  assert(x >= 0);
  assert(y >= 0);
  modifiable_img->img[y * modifiable_img->x_size + x] = value;
  return;
}

void set_pixel(Image *modifiable_img, int x, int y, int c, unsigned char value)
{
  assert(modifiable_img->channel_num > c);
  assert(modifiable_img->x_size > x);
  assert(modifiable_img->y_size > y);
  assert(c >= 0);
  assert(x >= 0);
  assert(y >= 0);
  modifiable_img->img[y * modifiable_img->x_size * modifiable_img->channel_num + x * modifiable_img->channel_num + c] = value;
  return;
}

Image copy(Image img)
{
  Image copy;
  copy.x_size = img.x_size;
  copy.y_size = img.y_size;
  copy.channel_num = img.channel_num;
  int img_size = sizeof(unsigned char) * img.x_size * img.y_size * img.channel_num;
  copy.img = malloc(img_size);
  memcpy(copy.img, img.img, img_size);
  return copy;
}

void fill_zeros(Image *modifiable_img)
{
  assert(modifiable_img->channel_num >= 3);
  for (int x = 0; x < modifiable_img->x_size; x++)
    for (int y = 0; y < modifiable_img->y_size; y++)
      for (int c = 0; c < 3; c++)
        set_pixel(modifiable_img, x, y, c, 0);
}

void free_img(Image img)
{
  free(img.img);
}

void median_filter(Image *modifiable_img, int kernel_size)
{
  assert(modifiable_img->channel_num >= 3);
  assert(kernel_size > 0);
  Image unmodified_image = copy(*modifiable_img);
  int kernal_size_t = kernel_size * kernel_size;
  unsigned char kernel[kernal_size_t];
  for (int x = 0; x < modifiable_img->x_size; x++)
    for (int y = 0; y < modifiable_img->y_size; y++)
      for (int c = 0; c < 3; c++)
      {
        memset(kernel, 0, kernal_size_t);
        for (int k_x = 0; k_x < kernel_size; k_x++)
        {
          int img_x = x - kernel_size / 2 + k_x;
          if (img_x < 0)
            img_x = -img_x;
          if (img_x >= modifiable_img->x_size)
            img_x -= img_x - modifiable_img->x_size + 1;
          for (int k_y = 0; k_y < kernel_size; k_y++)
          {
            int img_y = y - kernel_size / 2 + k_y;
            if (img_y < 0)
              img_y = -img_y;
            if (img_y >= modifiable_img->y_size)
              img_y -= img_y - modifiable_img->y_size + 1;
            kernel[k_x * kernel_size + k_y] = get_pixel(unmodified_image, img_x, img_y, c);
          }
        }
        unsigned char new_val = median_uchar(kernel, kernal_size_t);
        set_pixel(modifiable_img, x, y, c, new_val);
      }
  free_img(unmodified_image);
  return;
}

int compare_uchar(const void *a, const void *b)
{
  unsigned char arg1 = *(const unsigned char *)a;
  unsigned char arg2 = *(const unsigned char *)b;
  if (arg1 < arg2)
    return -1;
  if (arg1 > arg2)
    return 1;
  return 0;
}

unsigned char median_uchar(unsigned char *array, int size)
{
  if (size == 0)
    return 0;
  qsort(array, size, sizeof(unsigned char), compare_uchar);
  unsigned char result;
  if (size % 2 == 1)
    result = array[size / 2];
  else
    result = array[size / 2 - 1] / 2 + array[size / 2] / 2;
  return result;
}

void convolve(double modifiable_img[], int x_size, int y_size, double *kernel, int kernel_size, char edge_processing_mode)
{
  double *unmodifiet_mat = malloc(x_size * y_size * sizeof(double));
  memcpy(unmodifiet_mat, modifiable_img, x_size * y_size * sizeof(double));
  assert(kernel_size > 0);
  int kernal_size_t = kernel_size * kernel_size;
  for (int x = 0; x < x_size; x++)
    for (int y = 0; y < y_size; y++)
    {
      double result = 0;
      for (int k_x = 0; k_x < kernel_size; k_x++)
      {
        int img_x = x - kernel_size / 2 + k_x;
        if (img_x < 0)
        {
          if (edge_processing_mode == 0)
            continue;
          img_x = -img_x;
        }
        if (img_x >= x_size)
        {
          if (edge_processing_mode == 0)
            continue;
          img_x -= img_x - x_size + 1;
        }
        for (int k_y = 0; k_y < kernel_size; k_y++)
        {
          int img_y = y - kernel_size / 2 + k_y;
          if (img_y < 0)
          {
            if (edge_processing_mode == 0)
              continue;
            img_y = -img_y;
          }
          if (img_y >= y_size)
          {
            if (edge_processing_mode == 0)
              continue;
            img_y -= img_y - y_size + 1;
          }
          result += kernel[(kernel_size - k_x - 1) * kernel_size + (kernel_size - k_y - 1)] * unmodifiet_mat[img_y * x_size + img_x];
        }
      }
      modifiable_img[y * x_size + x] = result;
    }
  free(unmodifiet_mat);
}

void convolve_img(Image *modifiable_img, double *kernel, int kernel_size, char edge_processing_mode)
{
  // edge_processing modes:
  // 0 - zero everything beyond the edge
  // 1 - mirror content near edges
  Image unmodified_image = copy(*modifiable_img);
  double *result = malloc(modifiable_img->y_size * modifiable_img->x_size * sizeof(double));
  for (int c = 0; c < (modifiable_img->channel_num >= 3 ? 3 : 1); c++)
  {
    extract_channel(*modifiable_img, c, result);
    convolve(result, modifiable_img->x_size, modifiable_img->y_size, kernel, kernel_size, edge_processing_mode);
    for (int x = 0; x < modifiable_img->x_size; x++)
      for (int y = 0; y < modifiable_img->y_size; y++)
      {
        double res = result[y * modifiable_img->x_size + x];
        unsigned char new_val = res > 255 ? 255 : (res < 0 ? 0 : (unsigned char)res);
        if (modifiable_img->channel_num >= 3)
          set_pixel(modifiable_img, x, y, c, new_val);
        else
          set_pixel_grayscale(modifiable_img, x, y, new_val);
      }
  }
  free(result);
  return;
}

void gauss_filter(Image *modifiable_img, int kernel_size, double sigma, char is_deconvolve)
{
  double *kernel = (double *)malloc(sizeof(double) * kernel_size * kernel_size);
  double sum = 0;
  for (int k_x = 0; k_x < kernel_size; k_x++)
    for (int k_y = 0; k_y < kernel_size; k_y++)
    {
      int dx = k_x - kernel_size / 2;
      int dy = k_y - kernel_size / 2;
      double value = exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
      value /= (2 * 3.14 * sigma * sigma);
      kernel[k_y * kernel_size + k_x] = value;
      sum += value;
    }
  for (int i = 0; i < kernel_size * kernel_size; i++)
    kernel[i] /= sum;
  if (is_deconvolve == 0)
    convolve_img(modifiable_img, kernel, kernel_size, 1);
  else
  {
    Image unmodified_img = copy(*modifiable_img);
    deconvolve(unmodified_img, kernel, kernel_size, 10, modifiable_img);
    free_img(unmodified_img);
  }
  free(kernel);
  return;
}

void deconvolve(Image in_img, double *kernel, int kernel_size, int iterations, Image *out_img)
{

  assert(kernel_size > 0);
  assert(in_img.x_size == out_img->x_size);
  assert(in_img.y_size == out_img->y_size);
  assert(in_img.channel_num >= 3);
  assert(in_img.channel_num == out_img->channel_num);
  fill_zeros(out_img);
  float *estimate = malloc(in_img.x_size * in_img.y_size * sizeof(float));
  float *reblurred = malloc(in_img.x_size * in_img.y_size * sizeof(float));
  float *ratio = malloc(in_img.x_size * in_img.y_size * sizeof(float));
  for (int c = 0; c < 3; c++)
  {
    for (int x = 0; x < in_img.x_size; x++)
      for (int y = 0; y < in_img.y_size; y++)
        estimate[y * in_img.x_size + x] = (float)get_pixel(in_img, x, y, c);
    int pad = kernel_size / 2;
    for (int iter = 0; iter < iterations; iter++)
    {
      for (int y = 0; y < in_img.y_size; y++)
        for (int x = 0; x < in_img.x_size; x++)
        {
          reblurred[y * in_img.x_size + x] = 0;
          for (int ky = 0; ky < kernel_size; ky++)
            for (int kx = 0; kx < kernel_size; kx++)
            {
              int img_y = y + ky - pad;
              int img_x = x + kx - pad;

              if (img_y >= 0 && img_y < in_img.y_size && img_x >= 0 && img_x < in_img.x_size)

                reblurred[y * in_img.x_size + x] += estimate[img_y * in_img.x_size + img_x] *
                                                    kernel[ky * kernel_size + kx];
            }
        }
      for (int x = 0; x < in_img.x_size; x++)
        for (int y = 0; y < in_img.y_size; y++)
          ratio[y * in_img.x_size + x] = (float)get_pixel(in_img, x, y, c) / (reblurred[y * in_img.x_size + x] + 1e-6f);
      for (int y = 0; y < in_img.y_size; y++)
        for (int x = 0; x < in_img.x_size; x++)
        {
          float correction = 0;
          for (int ky = 0; ky < kernel_size; ky++)
            for (int kx = 0; kx < kernel_size; kx++)
            {
              int img_y = y + ky - pad;
              int img_x = x + kx - pad;

              if (img_y >= 0 && img_y < in_img.y_size && img_x >= 0 && img_x < in_img.x_size)

                correction += ratio[img_y * in_img.x_size + img_x] *
                              kernel[ky * kernel_size + kx];
            }

          estimate[y * in_img.x_size + x] *= correction;
        }
    }
    for (int x = 0; x < in_img.x_size; x++)
      for (int y = 0; y < in_img.y_size; y++)
      {
        float est = estimate[y * in_img.x_size + x];
        set_pixel(out_img, x, y, c, est < 0 ? 0 : (est > 255 ? 255 : (unsigned char)est));
      }
  }
  free(estimate);
  free(reblurred);
  free(ratio);
  return;
}

Image rgb2gray(Image in_img)
{
  assert(in_img.channel_num >= 3);
  Image out_img;
  out_img.x_size = in_img.x_size;
  out_img.y_size = in_img.y_size;
  out_img.channel_num = 1;
  out_img.img = malloc(in_img.x_size * in_img.y_size * sizeof(unsigned char));
  for (int x = 0; x < in_img.x_size; x++)
    for (int y = 0; y < in_img.y_size; y++)
    {
      float res = 0;
      res += (float)get_pixel(in_img, x, y, 0) * 0.299;
      res += (float)get_pixel(in_img, x, y, 1) * 0.587;
      res += (float)get_pixel(in_img, x, y, 2) * 0.114;
      res = res < 0 ? 0 : (res > 255 ? 255 : res);
      set_pixel_grayscale(&out_img, x, y, (unsigned char)res);
    }
  return out_img;
}

void extract_channel(Image img, int c_num, double *mat)
{
  assert(img.channel_num >= 3);
  for (int x = 0; x < img.x_size; x++)
    for (int y = 0; y < img.y_size; y++)
      mat[y * img.x_size + x] = (double)get_pixel(img, x, y, c_num);
}

Image extract_channel_img(Image img, int c_num)
{
  assert(img.channel_num >= 3);
  Image out_img;
  out_img.img = malloc(img.x_size * img.y_size * sizeof(unsigned char));
  out_img.x_size = img.x_size;
  out_img.y_size = img.y_size;
  out_img.channel_num = 1;
  for (int x = 0; x < img.x_size; x++)
    for (int y = 0; y < img.y_size; y++)
      set_pixel_grayscale(&out_img, x, y, get_pixel(img, x, y, c_num));
  return out_img;
}

void compose_img(Image r, Image g, Image b, Image *out_buffer)
{
  assert(out_buffer->channel_num >= 3);
  assert(r.channel_num == 1);
  assert(g.channel_num == 1);
  assert(b.channel_num == 1);
  assert(r.x_size == g.x_size);
  assert(g.x_size == b.x_size);
  assert(r.y_size == g.y_size);
  assert(g.y_size == b.y_size);
  out_buffer->x_size = r.x_size;
  out_buffer->y_size = r.y_size;
  for (int x = 0; x < out_buffer->x_size; x++)
    for (int y = 0; y < out_buffer->y_size; y++)
    {
      set_pixel(out_buffer, x, y, 0, get_pixel_grayscale(r, x, y));
      set_pixel(out_buffer, x, y, 1, get_pixel_grayscale(g, x, y));
      set_pixel(out_buffer, x, y, 2, get_pixel_grayscale(b, x, y));
    }
  return;
}

void img2mat(Image img, double *mat)
{
  assert(img.channel_num == 1);
  for (int x = 0; x < img.x_size; x++)
    for (int y = 0; y < img.y_size; y++)
    {
      double val = (double)get_pixel_grayscale(img, x, y);
      mat[y * img.x_size + x] = val;
    }
}

void mat2img(double *mat, Image *img)
{
  for (int x = 0; x < img->x_size; x++)
    for (int y = 0; y < img->y_size; y++)
    {
      double px_val = mat[y * img->x_size + x];
      set_pixel_grayscale(img, x, y, (px_val < 0 ? 0 : (px_val > 255 ? 255 : (unsigned char)px_val)));
    }
}

Image add(Image img1, Image img2)
{
  assert(img1.x_size == img2.x_size);
  assert(img1.y_size == img2.y_size);
  assert(img1.channel_num == img2.channel_num);
  Image out;
  out.x_size = img1.x_size;
  out.y_size = img1.y_size;
  out.channel_num = img1.channel_num;
  out.img = malloc(out.x_size * out.y_size * sizeof(unsigned char));
  for (int x = 0; x < img1.x_size; x++)
    for (int y = 0; y < img1.y_size; y++)
      for (int c = 0; c < (img1.channel_num >= 3 ? 3 : 1); c++)
      {
        int sm_res = 0;
        if (img1.channel_num >= 3)
          sm_res = (int)get_pixel(img1, x, y, c) + (int)get_pixel(img2, x, y, c);
        else
          sm_res = (int)get_pixel_grayscale(img1, x, y) + (int)get_pixel_grayscale(img2, x, y);
        if (img1.channel_num >= 3)
          set_pixel(&out, x, y, c, sm_res < 0 ? 0 : (sm_res > 255 ? 255 : (unsigned char)sm_res));
        else
          set_pixel_grayscale(&out, x, y, sm_res < 0 ? 0 : (sm_res > 255 ? 255 : (unsigned char)sm_res));
      }
  return out;
}

Image edge_detection(Image in_img, unsigned int low, unsigned int high)
{
  double kernel_x[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  double kernel_y[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
  const long long IMG_SIZE = in_img.y_size * in_img.x_size;
  Image gray_img = rgb2gray(in_img);
  double *result_x = malloc(IMG_SIZE * sizeof(double));
  img2mat(gray_img, result_x);
  convolve(result_x, in_img.x_size, in_img.y_size, kernel_x, 3, 1);
  double *result_y = malloc(IMG_SIZE * sizeof(double));
  img2mat(gray_img, result_y);
  convolve(result_y, in_img.x_size, in_img.y_size, kernel_y, 3, 1);
  double *angles = malloc(IMG_SIZE * sizeof(double));
  double *grad_res = malloc(IMG_SIZE * sizeof(double));
  for (int i = 0; i < IMG_SIZE; i++)
  {
    angles[i] = atan(result_y[i] / (result_x[i] == 0 ? 1e-10 : result_x[i]));
    double res = sqrt(result_x[i] * result_x[i] + result_y[i] * result_y[i]);
    grad_res[i] = res > 255 ? 255 : (res < 0 ? 0 : res);
  }

  for (int x = 1; x < in_img.x_size - 1; x++)
    for (int y = 1; y < in_img.y_size - 1; y++)
    {
      double act_angle = angles[y * in_img.x_size + x];
      double act_value = grad_res[y * in_img.x_size + x];
      double prev, next;
      if (act_angle > 0 && act_angle <= 3.14 / 4)
      {
        prev = grad_res[(y - 1) * in_img.x_size + (x + 1)];
        next = grad_res[(y + 1) * in_img.x_size + (x - 1)];
      }
      else if (act_angle > 3.14 / 4 && act_angle <= 3.14 / 2)
      {
        prev = grad_res[(y - 1) * in_img.x_size + x];
        next = grad_res[(y + 1) * in_img.x_size + x];
      }
      else if (act_angle > -3.14 / 4 && act_angle <= 0)
      {
        prev = grad_res[(y - 1) * in_img.x_size + (x - 1)];
        next = grad_res[(y + 1) * in_img.x_size + (x + 1)];
      }
      else if (act_angle == 0 || (act_angle > -3.14 / 2 && act_angle <= -3.14 / 4))
      {
        prev = grad_res[y * in_img.x_size + (x - 1)];
        next = grad_res[y * in_img.x_size + (x + 1)];
      }
      if (act_value < prev || act_value < next)
        grad_res[y * in_img.x_size + x] = 0;
    }

  for (int x = 1; x < in_img.x_size - 1; x++)
    for (int y = 1; y < in_img.y_size - 1; y++)
    {
      double act_value = grad_res[y * in_img.x_size + x];
      char changed = 0;
      if (act_value < low)
        grad_res[y * in_img.x_size + x] = 0;
      else if (act_value > high)
        grad_res[y * in_img.x_size + x] = 255;
      else
      {
        for (int m = -1; m <= 1; ++m)
        {
          for (int n = -1; n <= 1; ++n)
          {
            if (m == 0 && n == 0)
              continue;
            if (grad_res[(y + n) * in_img.x_size + x + m] > high)
            {
              grad_res[y * in_img.x_size + x] = 255;
              changed = 1;
              break;
            }
          }
          if (changed == 1)
            break;
        }
        if (changed == 0)
          grad_res[y * in_img.x_size + x] = 0;
      }
    }

  Image grad = {malloc(IMG_SIZE * sizeof(unsigned char)),
                in_img.x_size,
                grad.y_size = in_img.y_size,
                grad.channel_num = 1};
  mat2img(grad_res, &grad);
  free_img(gray_img);
  return grad;
}

void calc_hist(Image grayscale_img, double hist[256], double cdf[256])
{
  assert(hist != NULL);
  memset(hist, 0, 256 * sizeof(double));
  if (cdf != NULL)
    memset(cdf, 0, 256 * sizeof(double));
  for (int x = 0; x < grayscale_img.x_size; x++)
    for (int y = 0; y < grayscale_img.y_size; y++)
      hist[get_pixel_grayscale(grayscale_img, x, y)]++;
  for (int i = 0; i < 256; i++)
  {
    hist[i] /= grayscale_img.x_size * grayscale_img.y_size;
    if (cdf == NULL)
      continue;
    if (hist[i] != 0 && cdf[0] == 0)
      cdf[0] = hist[i];
    if (i == 0)
      cdf[i] = hist[i];
    else
      cdf[i] = cdf[i - 1] + hist[i];
  }
  return;
}

void equalize_hist(Image *modifiable_img)
{
  if (modifiable_img->channel_num >= 3)
  {
    equalize_hist_RGB(modifiable_img);
    return;
  }
  assert(modifiable_img->channel_num == 1);
  double hist[256];
  double cdf[256];
  calc_hist(*modifiable_img, hist, cdf);
  for (int x = 0; x < modifiable_img->x_size; x++)
    for (int y = 0; y < modifiable_img->y_size; y++)
    {
      unsigned char v = get_pixel_grayscale(*modifiable_img, x, y);
      double tmp = (cdf[v] - cdf[0]) / (1. - cdf[0]);
      unsigned char new_v = (unsigned char)(tmp * 255);
      set_pixel_grayscale(modifiable_img, x, y, new_v);
    }
}

void equalize_hist_RGB(Image *modifiable_img)
{
  Image r = extract_channel_img(*modifiable_img, 0);
  equalize_hist(&r);
  Image g = extract_channel_img(*modifiable_img, 1);
  equalize_hist(&g);
  Image b = extract_channel_img(*modifiable_img, 2);
  equalize_hist(&b);
  compose_img(r, g, b, modifiable_img);
  free_img(r);
  free_img(g);
  free_img(b);
  return;
}