#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Accelerate/Accelerate.h>


#define SIGMA 0.6f     // Standard deviation for Gaussian blur
#define TLOW  0.3f     // (Not used in this snippet)
#define THIGH 0.8f     // (Not used in this snippet)

// Function prototypes for PGM I/O
int readPGM(const char *filename, uint8_t **data, int *width, int *height);
int writePGM(const char *filename, uint8_t *data, int width, int height);

// Create a one-dimensional Gaussian kernel
void createGaussianKernel(float sigma, float **kernel, int *kernelLength) {
    int length = 1 + 2 * (int)ceil(2.5 * sigma);
    *kernelLength = length;
    *kernel = (float *)malloc(length * sizeof(float));
    int center = length / 2;
    float sum = 0.0f;
    for (int i = 0; i < length; i++) {
        float x = (float)(i - center);
        float value = expf(-0.5f * (x * x) / (sigma * sigma)) / (sigma * sqrtf(2 * M_PI));
        (*kernel)[i] = value;
        sum += value;
    }
    // Normalize kernel so that the sum equals 1.0
    for (int i = 0; i < length; i++) {
        (*kernel)[i] /= sum;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s input.pgm output.pgm\n", argv[0]);
        return -1;
    }

    uint8_t *inputData = NULL;
    int width, height;
    if (readPGM(argv[1], &inputData, &width, &height) != 0) {
        fprintf(stderr, "Error reading image %s\n", argv[1]);
        return -1;
    }
    int pixelCount = width * height;

    // Create vImage_Buffers for processing.
    // We'll convert the 8-bit image to float for processing.
    vImage_Buffer srcBuffer, floatBuffer, blurBuffer, dxBuffer, dyBuffer, magBuffer, outputBuffer;
    srcBuffer.data     = inputData;
    srcBuffer.width    = width;
    srcBuffer.height   = height;
    srcBuffer.rowBytes = width * sizeof(uint8_t);

    float *floatData = (float *)malloc(pixelCount * sizeof(float));
    float *blurData  = (float *)malloc(pixelCount * sizeof(float));
    float *dxData    = (float *)malloc(pixelCount * sizeof(float));
    float *dyData    = (float *)malloc(pixelCount * sizeof(float));
    float *magData   = (float *)malloc(pixelCount * sizeof(float));

    floatBuffer.data     = floatData;
    floatBuffer.width    = width;
    floatBuffer.height   = height;
    floatBuffer.rowBytes = width * sizeof(float);

    blurBuffer.data     = blurData;
    blurBuffer.width    = width;
    blurBuffer.height   = height;
    blurBuffer.rowBytes = width * sizeof(float);

    dxBuffer.data     = dxData;
    dxBuffer.width    = width;
    dxBuffer.height   = height;
    dxBuffer.rowBytes = width * sizeof(float);

    dyBuffer.data     = dyData;
    dyBuffer.width    = width;
    dyBuffer.height   = height;
    dyBuffer.rowBytes = width * sizeof(float);

    magBuffer.data     = magData;
    magBuffer.width    = width;
    magBuffer.height   = height;
    magBuffer.rowBytes = width * sizeof(float);

    // Convert input (Planar8) to PlanarF
    vImageConvert_Planar8toPlanarF(&srcBuffer, &floatBuffer, 0.0f, 255.0f, kvImageNoFlags);

    // --- Gaussian Blur ---
    // We perform separable convolution: first horizontal, then vertical.
    float *gaussKernel = NULL;
    int kernelLength = 0;
    createGaussianKernel(SIGMA, &gaussKernel, &kernelLength);

    // Temporary buffer for horizontal pass
    float *tempData = (float *)malloc(pixelCount * sizeof(float));
    vImage_Buffer tempBuffer;
    tempBuffer.data     = tempData;
    tempBuffer.width    = width;
    tempBuffer.height   = height;
    tempBuffer.rowBytes = width * sizeof(float);

    // Horizontal convolution: kernel of size [kernelLength x 1]
    vImageConvolve_PlanarF(&floatBuffer, &tempBuffer, NULL, 0, 0,
                           gaussKernel, kernelLength, 1, 0.0f, kvImageEdgeExtend);
    // Vertical convolution: kernel of size [1 x kernelLength]
    vImageConvolve_PlanarF(&tempBuffer, &blurBuffer, NULL, 0, 0,
                           gaussKernel, 1, kernelLength, 0.0f, kvImageEdgeExtend);

    free(gaussKernel);
    free(tempData);

    // --- Compute Image Derivatives ---
    // We use a simple [-1, 0, 1] kernel for both x and y derivatives.
    float kernel_dx[3] = { -1.0f, 0.0f, 1.0f };
    float kernel_dy[3] = { -1.0f, 0.0f, 1.0f };

    // Horizontal derivative (dx)
    vImageConvolve_PlanarF(&blurBuffer, &dxBuffer, NULL, 0, 0,
                           kernel_dx, 3, 1, 0.0f, kvImageEdgeExtend);
    // Vertical derivative (dy)
    vImageConvolve_PlanarF(&blurBuffer, &dyBuffer, NULL, 0, 0,
                           kernel_dy, 1, 3, 0.0f, kvImageEdgeExtend);

    // --- Compute Gradient Magnitude ---
    // Use vDSP to compute: magnitude = sqrt(dx^2 + dy^2)
    float *dxSq = (float *)malloc(pixelCount * sizeof(float));
    float *dySq = (float *)malloc(pixelCount * sizeof(float));
    vDSP_vsq(dxData, 1, dxSq, 1, pixelCount);  // dxSq = dx^2
    vDSP_vsq(dyData, 1, dySq, 1, pixelCount);  // dySq = dy^2
    vDSP_vadd(dxSq, 1, dySq, 1, magData, 1, pixelCount); // magData = dxSq + dySq
    int n = pixelCount;
        vvsqrtf(magData, magData, &n);
       // magData = sqrt(dxSq+dySq)
    free(dxSq);
    free(dySq);

    // (In a full Canny edge detector you would now perform non-maximal suppression
    // and hysteresis thresholding. These steps are not as easily vectorized and
    // are left as further work.)

    // --- Convert Result to 8-bit for Output ---
    // Normalize the gradient magnitude to [0, 255]
    uint8_t *outputData = (uint8_t *)malloc(pixelCount * sizeof(uint8_t));
    outputBuffer.data     = outputData;
    outputBuffer.width    = width;
    outputBuffer.height   = height;
    outputBuffer.rowBytes = width * sizeof(uint8_t);

    float maxVal;
    vDSP_maxv(magData, 1, &maxVal, pixelCount);
    if (maxVal < 1e-6f) { maxVal = 1.0f; }
    vDSP_vsdiv(magData, 1, &maxVal, magData, 1, pixelCount);
    float scale = 255.0f;
    vDSP_vsmul(magData, 1, &scale, magData, 1, pixelCount);
    vImageConvert_PlanarFtoPlanar8(&magBuffer, &outputBuffer, 0.0f, 255.0f, kvImageNoFlags);

    if (writePGM(argv[2], outputData, width, height) != 0) {
        fprintf(stderr, "Error writing output image %s\n", argv[2]);
    }

    // Clean up
    free(inputData);
    free(floatData);
    free(blurData);
    free(dxData);
    free(dyData);
    free(magData);
    free(outputData);

    return 0;
}

// --- PGM I/O Functions ---

// Read a binary PGM (P5) file.
int readPGM(const char *filename, uint8_t **data, int *width, int *height) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("fopen");
        return -1;
    }
    char header[3];
    if (fscanf(fp, "%2s", header) != 1) {
        fclose(fp);
        return -1;
    }
    if (strcmp(header, "P5") != 0) {
        fprintf(stderr, "Not a P5 PGM file.\n");
        fclose(fp);
        return -1;
    }
    // Skip comments
    int c = fgetc(fp);
    while (c == '#') {
        while (fgetc(fp) != '\n');
        c = fgetc(fp);
    }
    ungetc(c, fp);
    int maxval;
    if (fscanf(fp, "%d %d %d", width, height, &maxval) != 3) {
        fclose(fp);
        return -1;
    }
    fgetc(fp); // consume newline
    int size = (*width) * (*height);
    *data = (uint8_t *)malloc(size * sizeof(uint8_t));
    if (fread(*data, sizeof(uint8_t), size, fp) != (size_t)size) {
        fclose(fp);
        free(*data);
        return -1;
    }
    fclose(fp);
    return 0;
}

// Write a binary PGM (P5) file.
int writePGM(const char *filename, uint8_t *data, int width, int height) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("fopen");
        return -1;
    }
    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    int size = width * height;
    if (fwrite(data, sizeof(uint8_t), size, fp) != (size_t)size) {
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}
