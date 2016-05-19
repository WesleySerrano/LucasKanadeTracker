/**
@file utils.h
Application auxiliary functions

@brief Feature tracking using Lucas-Kanade algorithm

@date 18 May 2016

@author Wesley Serrano
*/

#include "CImg.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cassert>

using namespace cimg_library;
using namespace std;

typedef struct _point
{
  double x;
  double y;
} Point;

typedef struct _vector2
{
  double x;
  double y;
} Vector2;

/**
  Returns the maximum of three values of the colors channels of a pixel

  @param r The value of the red channel
  @param g The value of the green channel
  @param b The value of the blue channel

  @return The maximum of the color channels
*/
double MaxTone(double r,double g,double b)
{
  double aux = max(r,g);
  return max(aux,b);
}

/**
  Returns the minimum of three values of the colors channels of a pixel

  @param r The value of the red channel
  @param g The value of the green channel
  @param b The value of the blue channel

  @return The minimum of the color channels
*/
float MinTone(float r,float g,float b)
{
  float aux = min(r,g);
  return min(aux,b);
}

/**
  Returns the average of three values of the colors channels of a pixel

  @param r The value of the red channel
  @param g The value of the green channel
  @param b The value of the blue channel

  @return The average of the color channels
*/
float AvgTone(float r,float g,float b)
{
  return (r + g + b)/3.0;
}

/**
  Creates a gray scale version of an image

  @param image Image to be transformed
  @param scale Indicates whether or not the RGB channels values should be in the range [0, 1] or [0, 255] (as integer values)


  @return Image in gray scale
*/
CImg<double> makeImageGray(CImg<double> image, bool scale = true)
{
  const bool scaleFactor = scale ? 255 : 1;
  CImg<double> result = image;

  cimg_forXY(result, x, y)
  {
    const double r = result(x, y, 0, 0)/scaleFactor;
    const double g = result(x, y, 0, 1)/scaleFactor;
    const double b = result(x, y, 0, 2)/scaleFactor;

    result(x, y, 0, 0) = 0.299 * r + 0.587 * g + 0.114 * b;
    result(x, y, 0, 1) = 0.299 * r + 0.587 * g + 0.114 * b;
    result(x, y, 0, 2) = 0.299 * r + 0.587 * g + 0.114 * b;
  }

  return result;
}


/**
  Creates a window to display a image given as input parameter
  with green dots on the points given as the other parameter

  @param image Image to be marked and shown on a window
  @param pointsToMark The coordinates of the points to be marked on the image
*/
void markPointsOnImage(CImg<double> image, vector<Point> pointsToMark)
{
  CImg<double> points = image;

  for(vector<Point>::iterator it = pointsToMark.begin(); it != pointsToMark.end(); ++it)
  {
    const double x = it->x;
    const double y = it->y;

    points(x, y, 0, 0) = 0;
    points(x, y, 0, 1) = 255;
    points(x, y, 0, 2) = 0;
  }

  displayImage(points);
}

/**
  Creates a window to display a image given as input parameter

  @param image Image to be shown on a window
*/

void displayImage(CImg<double> image)
{
  CImgDisplay display(image, "Image");

  while (!display.is_closed()) 
  {
    display.wait();
  }
}

/**
  Creates a window to display a vector of images given as input parameter.
  Displays one image at a time

  @param images The vector containing the images to be shown on a window
*/
void displayImages(vector< CImg<double> > images)
{
  for(int i = 0; i < images.size(); i++) displayImage(images[i]);
}

/**
  Sums the rgb channels of two images and save the result on another image.
  If one of the values is greater than 255, it is changed to 255.

  @param image1 First image to be used
  @param image2 Second image to be used 

  @return Image containing the result of the sum
*/
CImg<double> add(CImg<double> image1, CImg<double> image2)
{
   const int WIDTH1 = image1.width(), HEIGHT1 = image1.height();
   const int WIDTH2 = image2.width(), HEIGHT2 = image2.height();

   //assert(WIDTH1 == WIDTH2 && HEIGHT1 == HEIGHT2);

   CImg<double> newImage(WIDTH1,HEIGHT1,1,3,0);
   
   cimg_forXY(image1,x,y)
   {
     const double r1 = image1(x,y,0,0), g1 = image1(x,y,0,1), b1 = image1(x,y,0,2);
     const double r2 = image2(x,y,0,0), g2 = image2(x,y,0,1), b2 = image2(x,y,0,2);

     newImage(x,y,0,0) = r1 + r2;
     newImage(x,y,0,1) = g1 + g2;
     newImage(x,y,0,2) = b1 + b2;
   }

   return newImage;
}

/**
  Subtracts the rgb channels of two images and save the result on another image.
  For every color channel of every pixel, the absolute value of the difference is saved. 

  @param image1 First image to be used
  @param image2 Second image to be used 

  @return Image containing the result of the subtraction
*/
CImg<double> subtract(CImg<double> image1, CImg<double> image2)
{
   const int WIDTH1 = image1.width(), HEIGHT1 = image1.height();
   const int WIDTH2 = image2.width(), HEIGHT2 = image2.height();

   //assert(WIDTH1 == WIDTH2 && HEIGHT1 == HEIGHT2);


   CImg<double> newImage(WIDTH1,HEIGHT1,1,3,0);
   
   cimg_forXY(image1,x,y)
   {
     const double r1 = image1(x,y,0,0), g1 = image1(x,y,0,1), b1 = image1(x,y,0,2);
     const double r2 = image2(x,y,0,0), g2 = image2(x,y,0,1), b2 = image2(x,y,0,2);

     newImage(x,y,0,0) = abs(r1 - r2);
     newImage(x,y,0,1) = abs(g1 - g2);
     newImage(x,y,0,2) = abs(b1 - b2);
   }

   return newImage;
}

CImg<double> unidimensionalConvolution(CImg<double> image, double* convolutionMatrix, const int matrixSize, const bool vertical, const bool grayScale = false)
{
   const int offset = (int) matrixSize/2;
   const int WIDTH = image.width(), HEIGHT = image.height();
   CImg<double> resultImage(WIDTH,HEIGHT,1, grayScale ? 1 : 3,0);

   cimg_forXY(image,x,y)
   {
     double imageR, imageG, imageB;
     double imageC;

     double r = 0.0, g = 0.0, b = 0.0;
     double c = 0.0;

     for(int i = -offset; i <= offset; i++)
     {
      const int pixelY = vertical ? y + i : y;
      const int matrixCell = i + offset;
      bool pixelOutOfImage;
      
      const int pixelX = vertical? x : x + i;

      if((pixelY < 0 || pixelY >= HEIGHT) || (pixelX < 0 || pixelX >= WIDTH)) pixelOutOfImage = true;
      else pixelOutOfImage = false;

      if(grayScale)
      {
        imageC = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
        c += imageC * convolutionMatrix[matrixCell];
      }
      else
      {
       imageR = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
       imageG = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,1); 
       imageB = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,2); 

       r += imageR * convolutionMatrix[matrixCell];
       g += imageG * convolutionMatrix[matrixCell];
       b += imageB * convolutionMatrix[matrixCell];
      }
      
     }

     if(grayScale) resultImage(x,y,0,0) = c;
     else
     { 
      resultImage(x,y,0,0) = r;
      resultImage(x,y,0,1) = g;
      resultImage(x,y,0,2) = b;
     }
   }   

   return resultImage;
}

/**
  Convolves an image with a matrix
  @param image Image to apply to convolution
  @param convolutionMatrix The convolution matrix to be applied. Must be square.
  @param matrixSize The square matrix dimension
  @param grayScale Indicates Whether or not the image to be convolved is a gray scale image

  @return The convolved image
*/
CImg<double> convolve(CImg<double> image, double** convolutionMatrix, const int matrixSize, const bool grayScale = false)
{
   const int offset = (int) matrixSize/2;
   const int WIDTH = image.width(), HEIGHT = image.height();
   CImg<double> resultImage(WIDTH,HEIGHT,1, grayScale ? 1 : 3,0);

   cimg_forXY(image,x,y)
   {
     double imageR, imageG, imageB;
     double imageC;

     double r = 0.0, g = 0.0, b = 0.0;
     double c = 0.0;

     for(int i = -offset; i <= offset; i++)
     {
      const int pixelY = y + i;
      const int row = i + offset;
      bool pixelOutOfImage;

      for(int j = -offset; j <= offset; j++)
      {
        const int pixelX = x + j;
        const int column = j + offset;

        if((pixelY < 0 || pixelY >= HEIGHT) || (pixelX < 0 || pixelX >= WIDTH)) pixelOutOfImage = true;
        else pixelOutOfImage = false;

        if(grayScale)
        {
          imageC = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
          c += imageC * convolutionMatrix[row][column];
        }
        else
        {
         imageR = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
         imageG = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,1); 
         imageB = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,2); 

         r += imageR * convolutionMatrix[row][column];
         g += imageG * convolutionMatrix[row][column];
         b += imageB * convolutionMatrix[row][column];
        }
      }
     }

     if(grayScale) resultImage(x,y,0,0) = c;
     else
     { 
      resultImage(x,y,0,0) = r;
      resultImage(x,y,0,1) = g;
      resultImage(x,y,0,2) = b;
     }
   }   

   return resultImage;
}

/**
  Creates a Gaussian filter matrix
  @param filterSize The matrix dimension. Must be odd
  @param delta The delta parameter of the gaussian function
  @param correction Indicates whether or not the matrix values should be adjusted

  @return The filter matrix
*/
double** makeGaussianFilter(const int filterSize, const double delta, const bool correction = false)
{
   assert(filterSize%2 == 1);

   double** filterMatrix;
   filterMatrix = new double*[filterSize];
   for (int i = 0; i < filterSize; i++) filterMatrix[i] = new double[filterSize];

   double sum = 0.0;

  for (int i = 0; i < filterSize; i++)
  {
    for (int j = 0; j < filterSize; j++)
    {
      int i2 = i-(filterSize-2);
      int j2 = j-(filterSize-2);
      double value = (exp((-(i2*i2+j2*j2))/(2*delta*delta)))/(2*M_PI*delta*delta);
      filterMatrix[j][i] = value;
      sum += value;
    }
  }

  if(correction) sum = sum/2.0;

  for (int i = 0; i < filterSize; i++)
  {
    for (int j = 0; j < filterSize; j++)
    {
      double value2 = filterMatrix[j][i]/sum;
      filterMatrix[j][i] = value2;
    }
  }

   return filterMatrix;
}

double** makeConstantGaussianFilter(const double correction)
{
   double** filterMatrix;
   filterMatrix = new double*[5];
   for (int i = 0; i < 5; i++) filterMatrix[i] = new double[5];
   double filterWeight = correction ? 8.0 : 16.0;
    filterWeight *= filterWeight;

   filterMatrix[0][0] = 1.0/filterWeight; filterMatrix[0][1] = 4.0/filterWeight; filterMatrix[0][2] = 6.0/filterWeight; filterMatrix[0][3] = 4.0/filterWeight; filterMatrix[0][4] = 1.0/filterWeight;
   filterMatrix[1][0] = 4.0/filterWeight; filterMatrix[1][1] = 16.0/filterWeight; filterMatrix[1][2] = 24.0/filterWeight; filterMatrix[1][3] = 16.0/filterWeight; filterMatrix[1][4] = 4.0/filterWeight;
   filterMatrix[2][0] = 6.0/filterWeight; filterMatrix[2][1] = 24.0/filterWeight; filterMatrix[2][2] = 36.0/filterWeight; filterMatrix[2][3] = 24.0/filterWeight; filterMatrix[2][4] = 6.0/filterWeight;
   filterMatrix[3][0] = 4.0/filterWeight; filterMatrix[3][1] = 16.0/filterWeight; filterMatrix[3][2] = 24.0/filterWeight; filterMatrix[3][3] = 16.0/filterWeight; filterMatrix[3][4] = 4.0/filterWeight;
   filterMatrix[2][0] = 1.0/filterWeight; filterMatrix[4][1] = 4.0/filterWeight; filterMatrix[4][2] = 6.0/filterWeight; filterMatrix[4][3] = 4.0/filterWeight; filterMatrix[4][4] = 1.0/filterWeight;

   return filterMatrix;
}

double* makeConstantUnidimensionalGaussianFilter(const double correction)
{
   double* filterMatrix;
   filterMatrix = new double[5];

   double filterWeight = correction ? 8.0 : 16.0;

   filterMatrix[0] = 1.0/filterWeight; filterMatrix[1] = 4.0/filterWeight; filterMatrix[2] = 6.0/filterWeight; filterMatrix[3] = 4.0/filterWeight; filterMatrix[4] = 1.0/filterWeight;

   return filterMatrix;
}

/**
  Applies a 3x3 Gaussian filter on a image
  @param image image to apply the filter
  @param correction Wether a correction factor should or not be applied to the result
  @param delta The delta parameter in the gaussian filter generation

  @return The original image blurred
*/

CImg<double> blur3x3(CImg<double> image, const bool correction = false, const double delta = 1.5, const bool grayScale = false)
{
   const int filterSize = 3;

   double** convolutionMatrix;
   convolutionMatrix = makeGaussianFilter(filterSize, delta);
   
   return convolve(image, convolutionMatrix, filterSize, grayScale);
}

/**
  Applies a 5x5 Gaussian filter on a image
  @param image image to apply the filter
  @param correction Wether a correction factor should or not be applied to the result
  @param delta The delta parameter in the gaussian filter generation

  @return The original image blurred
*/

CImg<double> blur5x5(CImg<double> image, const bool correction = false, const double delta = 1.5, const bool grayScale = false)
{   
   const int filterSize = 5;

   double** convolutionMatrix;
   convolutionMatrix = makeGaussianFilter(filterSize, delta, correction);

   //convolutionMatrix = makeConstantGaussianFilter(correction);

   return convolve(image, convolutionMatrix, filterSize, grayScale);
   
   /*double* convolutionMatrix;
   convolutionMatrix = makeConstantUnidimensionalGaussianFilter(correction);

   CImg<double> verticallyFiltered = unidimensionalConvolution(image, convolutionMatrix, filterSize, true, grayScale);

   CImg<double> horizontallyFiltered = unidimensionalConvolution(verticallyFiltered, convolutionMatrix, filterSize, false, grayScale);
   

   return horizontallyFiltered;*/
}

/**
  Applies a Gaussian filter on a image. The filter matrix size depends on a parameter. If the filter size is not 3 or 5, it does nothing.
  @param image The image to apply the filter
  @param filterSize The filter matrix size.
  @param correction Wether a correction factor should or not be applied to the result
  @param delta The delta parameter in the gaussian filter generation

  @return The original image blurred
*/
CImg<double> blur(CImg<double> image, const int filterSize, const bool correction = false, const double delta = 2.0, const bool grayScale = false)
{
  if(filterSize == 3) return blur3x3(image, correction, delta, grayScale);
  else if(filterSize == 5) return blur5x5(image, correction, delta, grayScale);
}

/**
  Double
  @param image First image to be expanded

  @return The expanded image
*/
CImg<double> expandImage(CImg<double> image)
{
   const int WIDTH = image.width(), HEIGHT = image.height();
   
   CImg<double> expandedImage(2*WIDTH, 2*HEIGHT,1,3,0);

   cimg_forXY(expandedImage,x,y)
   {
     if(x%2 == 0 && y%2 == 0)
     {
       const int px = x/2, py = y/2;
       double r = image(px,py,0,0), g = image(px,py,0,1), b = image(px,py,0,2);   
       
       expandedImage(x,y,0,0) = r;
       expandedImage(x,y,0,1) = g;
       expandedImage(x,y,0,2) = b;
     }
   }

   return expandedImage;
}

/**
  Reduces a image by half
  @param image First image to be reduced

  @return The reduced image
*/
CImg<double> reduceImage(CImg<double> image, const bool grayScale = false)
{
   const int WIDTH = image.width(), HEIGHT = image.height();
   
   CImg<double> reducedImage(WIDTH%2 == 0 ? WIDTH/2 : (WIDTH + 1)/2,HEIGHT%2 == 0 ? HEIGHT/2 : (HEIGHT + 1)/2,1, grayScale ? 1 : 3,0);

   cimg_forXY(image,x,y)
   {
     if(x%2 == 0 && y%2 == 0)
     {
       if(grayScale)
       {
        double c = image(x,y,0,0);
        reducedImage(x/2,y/2,0,0) = c;        
       }
       else
       {
        double r = image(x,y,0,0), g = image(x,y,0,1), b = image(x,y,0,2);   
       
        reducedImage(x/2,y/2,0,0) = r;
        reducedImage(x/2,y/2,0,1) = g;
        reducedImage(x/2,y/2,0,2) = b;
       }
     }
   }

   return reducedImage;
}

/**
  Makes a gaussian pyramid of images from an original image. 
  
  @param imageFilePath The path to the image
  @param numberOfImages Output parameter to save the number of generated images from the original
  @param gaussianFilterSize The size of the gaussian filter to be used

  @return Image containing all the pyramid images 
*/
vector< CImg<double> > makeGaussianPyramid(CImg<double> image, const int gaussianFilterSize, const bool grayScale = false)
{
  vector< CImg<double> > pyramid;

  const int MINIMUM_WIDTH = 32;
  const double DELTA = 2.0;
  const int width = image.width(), height = image.height();
    
  const int numberOfImages = (int)(log2(width) - log2(MINIMUM_WIDTH));
  
  CImg<double> firstBlur = blur(image, gaussianFilterSize, false, DELTA, grayScale);
  pyramid.push_back(grayScale ? firstBlur : image);

  //CImg<double> blurredImage = blur(pyramid[0], gaussianFilterSize, false, DELTA, grayScale);

  //CImg<double> silhouetteImage = subtract(pyramid[0], blurredImage);

  for(int i = 1; i < numberOfImages; i++)
  { 
    CImg<double>  blurredImage = blur(pyramid[i - 1], gaussianFilterSize, false, DELTA, grayScale);
    
    pyramid.push_back(reduceImage(blurredImage, grayScale));
  }

  return pyramid;
}
