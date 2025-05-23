#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
typedef struct
{
  unsigned char *img;
  int x_size;
  int y_size;
  int channel_num;
} Image;
/**
 * @brief Loads an image from a file
 * @param path Path to the image file
 * @param desired_channel_num Desired number of color channels (0 for auto)
 * @return Image structure containing the loaded image data
 */
Image load_img(unsigned char *path, int desired_channel_num);

/**
 * @brief Saves an image as a PNG file
 * @param img Image to save
 * @param path Output file path
 * @param stride_bytes Number of bytes between rows in the image
 * @return Non-zero on success, zero on failure
 */
int save_png(Image img, unsigned char *path, int stride_bytes);

/**
 * @brief Saves an image as a BMP file
 * @param img Image to save
 * @param path Output file path
 * @return Non-zero on success, zero on failure
 */
int save_bmp(Image img, unsigned char *path);

/**
 * @brief Saves an image as a JPG file
 * @param img Image to save
 * @param path Output file path
 * @param quality JPEG quality (1-100)
 * @return Non-zero on success, zero on failure
 */
int save_jpg(Image img, unsigned char *path, int quality);

/**
 * @brief Gets a grayscale pixel value from an image
 * @param img Image to read from
 * @param x X coordinate
 * @param y Y coordinate
 * @return Pixel value (0-255)
 */
unsigned char get_pixel_grayscale(Image img, int x, int y);

/**
 * @brief Gets a pixel value from a specific channel of an image
 * @param img Image to read from
 * @param x X coordinate
 * @param y Y coordinate
 * @param c Channel index
 * @return Pixel value (0-255)
 */
unsigned char get_pixel(Image img, int x, int y, int c);

/**
 * @brief Sets a grayscale pixel value in an image
 * @param modifiable_img Pointer to image to modify
 * @param x X coordinate
 * @param y Y coordinate
 * @param value New pixel value (0-255)
 */
void set_pixel_grayscale(Image *modifiable_img, int x, int y, unsigned char value);

/**
 * @brief Sets a pixel value in a specific channel of an image
 * @param modifiable_img Pointer to image to modify
 * @param x X coordinate
 * @param y Y coordinate
 * @param c Channel index
 * @param value New pixel value (0-255)
 */
void set_pixel(Image *modifiable_img, int x, int y, int c, unsigned char value);

/**
 * @brief Creates a deep copy of an image
 * @param img Image to copy
 * @return New image with copied data
 */
Image copy(Image img);

/**
 * @brief Fills the first 3 channels of an image with zeros
 * @param modifiable_img Pointer to image to modify
 */
void fill_zeros(Image *modifiable_img);

/**
 * @brief Frees the memory allocated for an image
 * @param img Image to free
 */
void free_img(Image img);

/**
 * @brief Applies median filter to an image
 * @param modifiable_img Pointer to image to modify
 * @param kernel_size Size of the kernel (must be odd)
 */
void median_filter(Image *modifiable_img, int kernel_size);

/**
 * @brief Comparison function for unsigned chars (used in qsort)
 * @param a First value to compare
 * @param b Second value to compare
 * @return -1 if a < b, 1 if a > b, 0 if equal
 */
int compare_uchar(const void *a, const void *b);

/**
 * @brief Calculates median of an array of unsigned chars
 * @param array Array of values
 * @param size Size of the array
 * @return Median value
 */
unsigned char median_uchar(unsigned char *array, int size);

/**
 * @brief Performs convolution on a 2D matrix
 * @param modifiable_img Matrix to convolve (modified in-place)
 * @param x_size Width of the matrix
 * @param y_size Height of the matrix
 * @param kernel Convolution kernel
 * @param kernel_size Size of the kernel
 * @param edge_processing_mode 0 for zero padding, 1 for mirroring
 */
void convolve(double modifiable_img[], int x_size, int y_size, double *kernel, int kernel_size, char edge_processing_mode);

/**
 * @brief Performs convolution on an image
 * @param modifiable_img Pointer to image to modify
 * @param kernel Convolution kernel
 * @param kernel_size Size of the kernel
 * @param edge_processing_mode 0 for zero padding, 1 for mirroring
 */
void convolve_img(Image *modifiable_img, double *kernel, int kernel_size, char edge_processing_mode);

/**
 * @brief Applies Gaussian filter to an image
 * @param modifiable_img Pointer to image to modify
 * @param kernel_size Size of the kernel
 * @param sigma Standard deviation for Gaussian
 * @param is_deconvolve 0 for normal convolution, 1 for deconvolution
 */
void gauss_filter(Image *modifiable_img, int kernel_size, double sigma, char is_deconvolve);

/**
 * @brief Performs deconvolution on an image
 * @param in_img Input image
 * @param kernel Blur kernel
 * @param kernel_size Size of the kernel
 * @param iterations Number of iterations to perform
 * @param out_img Output image (modified in-place)
 */
void deconvolve(Image in_img, double *kernel, int kernel_size, int iterations, Image *out_img);

/**
 * @brief Converts an RGB image to grayscale
 * @param in_img Input RGB image
 * @return Grayscale image
 */
Image rgb2gray(Image in_img);

/**
 * @brief Extracts a single channel from an image into a matrix
 * @param img Input image
 * @param c_num Channel index to extract
 * @param mat Output matrix (must be pre-allocated)
 */
void extract_channel(Image img, int c_num, double *mat);
/**
 * Extracts a single channel from a multi-channel image and returns it as a grayscale image.
 *
 * @param img The input image with at least 3 channels (asserts img.channel_num >= 3)
 * @param c_num The channel number to extract (0 for Red, 1 for Green, 2 for Blue in RGB image)
 * @return Image A new single-channel grayscale image containing only the specified channel
 *
 * @note The function allocates memory for the output image which must be freed by the caller
 * @warning The input image must have at least 3 channels or the function will assert
 */
Image extract_channel_img(Image img, int c_num);
/**
 * Composes a multi-channel image from three single-channel grayscale images (R, G, B).
 *
 * @param r Red channel as a single-channel grayscale image
 * @param g Green channel as a single-channel grayscale image
 * @param b Blue channel as a single-channel grayscale image
 * @param out_buffer Pointer to the output image buffer (must have >= 3 channels)
 *
 * @note All input images must be single-channel and have identical dimensions
 * @warning The function asserts that:
 *          - out_buffer has >= 3 channels
 *          - All input images are single-channel
 *          - All input images have matching dimensions
 *          - out_buffer is pre-allocated with sufficient capacity
 */
void compose_img(Image r, Image g, Image b, Image *out_buffer);
/**
 * @brief Converts a grayscale image to a matrix
 * @param img Input grayscale image
 * @param mat Output matrix (must be pre-allocated)
 */
void img2mat(Image img, double *mat);

/**
 * @brief Converts a matrix to a grayscale image
 * @param mat Input matrix
 * @param img Output image (modified in-place)
 */
void mat2img(double *mat, Image *img);

/**
 * @brief Adds two images together
 * @param img1 First image
 * @param img2 Second image
 * @return New image containing the sum
 */
Image add(Image img1, Image img2);

/**
 * @brief Performs edge detection using Canny algorithm
 * @param in_img Input image
 * @param low Low threshold for hysteresis
 * @param high High threshold for hysteresis
 * @return Image with detected edges
 */
Image edge_detection(Image in_img, unsigned int low, unsigned int high);
/**
 * Calculates the normalized histogram and (optionally) the cumulative distribution function (CDF)
 * of a grayscale image.
 *
 * @param grayscale_img Input single-channel grayscale image
 * @param hist Pre-allocated array of 256 doubles to store the normalized histogram
 *             (counts normalized by total pixels, so values range [0,1])
 * @param cdf Pre-allocated array of 256 doubles to store the CDF, or NULL if not needed
 *            (values range [0,1] if provided)
 *
 * @note The function will:
 *       - Normalize the histogram by dividing by total pixels
 *       - Compute CDF if cdf array is provided (NULL skips CDF calculation)
 *       - Initialize both arrays to zero before computation
 * @warning Asserts that hist is not NULL
 */
void calc_hist(Image grayscale_img, double hist[256], double cdf[256]);
/**
 * Performs histogram equalization on a grayscale or RGB image in-place.
 * For RGB images, delegates to equalize_hist_RGB().
 *
 * @param modifiable_img Pointer to the image to be modified (either grayscale or RGB)
 *                       For RGB: processes each channel independently in HSV color space
 *                       For grayscale: performs standard histogram equalization
 *
 * @note The function:
 *       - Modifies the input image directly (in-place operation)
 *       - Handles both single-channel and multi-channel images
 *       - For RGB, converts to HSV space, equalizes V channel, then converts back
 * @warning Asserts the image is single-channel when processing grayscale
 */
void equalize_hist(Image *modifiable_img);
/**
 * Helper function that performs histogram equalization on an RGB image in-place
 * by converting to HSV color space and equalizing only the Value (V) channel.
 *
 * @param modifiable_img Pointer to an RGB image (asserts channel_num >= 3)
 *                       The image will be modified directly
 *
 * @note The function:
 *       - Preserves hue and saturation while equalizing value
 *       - Converts RGB→HSV and back HSV→RGB
 *       - Uses equalize_hist() for the V channel processing
 * @warning Expects a 3+ channel image (typically RGB)
 * @warning Modifies the input image directly
 */
void equalize_hist_RGB(Image *modifiable_img);
#endif