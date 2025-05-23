# Image Processing Console Tool

A lightweight command-line utility for applying various image processing operations.

## Features

- **Filters**:
  - Median filter (noise reduction)
  - Gaussian blur
  - Gaussian deconvolution (sharpening)
  - Edge detection (Canny-like)
- **Color Operations**:
  - Grayscale conversion
  - Histogram equalization
- **Format Support**:
  - Input: Any format supported by stb_image
  - Output: PNG, JPG, BMP

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/image-processor.git
   cd image-processor
   ```

2. Compile (requires C compiler):

   ```bash
   gcc main.c -o image_processor -lm -std=c99
   ```

## Usage

```bash
./image_processor <input_image> <mode> [mode_args...] <output_image>.<format>
```

### Modes

| Mode            | Arguments                    | Description                         |
| --------------- | ---------------------------- | ----------------------------------- |
| `-median`       | `<kernel_size>`              | Median filter for noise reduction   |
| `-gauss`        | `<k_size> <sigma>`           | Gaussian blur                       |
| `-gauss-deconv` | `<k_size> <sigma>`           | Gaussian deconvolution (sharpening) |
| `-edges`        | `<low_thresh> <high_thresh>` | Edge detection                      |
| `-grayscale`    | (none)                       | Converts to grayscale               |
| `-equilize`     | (none)                       | Histogram equalization              |

### Examples

```bash
# Apply median filter with 5x5 kernel
./image_processor input.jpg -median 5 output.png

# Gaussian blur with 7x7 kernel (σ=1.5)
./image_processor input.png -gauss 7 1.5 output.jpg

# Edge detection with thresholds 50/150
./image_processor input.bmp -edges 50 150 output.png

# Convert to grayscale
./image_processor input.jpg -grayscale output.bmp

# Histogram equalization
./image_processor input.png -equilize output.png
```

more examples could be found in **examples** folder

## Technical Notes

- Kernel sizes must be **odd integers** (3, 5, 7...)
- Edge thresholds range: **0-255**
- For RGB images, histogram equalization works in HSV color space (V channel only)

## Dependencies

- C standard library
- [stb_image](https://github.com/nothings/stb) (included)
- Math library (`-lm` flag when compiling)

## License

[MIT License](LICENSE)
