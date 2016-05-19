/**
@file main.cpp
Application main file. Receives user input and runs the tracker algorithm

@brief Feature tracking using Lucas-Kanade algorithm

@date 18 May 2016

@author Wesley Serrano
*/

#include "utils.h"

/**
  Calculates the inverse of a 2x2 matrix

  @param matrix The matrix to be inverted

  @return The inverse matrix
*/
double** invertMatrix(double** matrix)
{
  const int matrixSize = 2;
  double** inverse;
  inverse = new double*[matrixSize];
  for (int i = 0; i < matrixSize; i++)
  {
    inverse[i] = new double[matrixSize];
  }

  const double determinant = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);

  inverse[0][0] = matrix[1][1]/determinant;
  inverse[0][1] = -matrix[0][1]/determinant;
  inverse[1][0] = -matrix[1][0]/determinant;
  inverse[1][1] = matrix[0][0]/determinant;

  return inverse;
}

/**
  Calculates the matrix (...)

  @param gradientX The horizontal gradient for each pixel of an image
  @param gradientY The vertical gradient for each pixel of an image
  @param x The x coordinate of the pixel
  @param y The y coordinate of the pixel
  @param offset The distance from the pixel to its furthest horizontal (and vertical) neighbours

  @return The matrix (...)
*/
double** calculateMatrix(CImg<double> gradientX, CImg<double> gradientY, int x, int y, int offset)
{
  const int matrixSize = 2;
  double** matrix;
  matrix = new double*[matrixSize];
  for (int i = 0; i < matrixSize; i++)
  {
    matrix[i] = new double[matrixSize];
    for (int j = 0; j < matrixSize; j++)
    {
      matrix[i][j] = 0;
    }
  }

  for(int i = -offset; i <= offset; i++)
  {
    const int py = y + i;
    for(int j = -offset; j <= offset; j++)
    {
      const int px = x + j;
      //const double delta = sqrt(2.0);
      //const double squareDelta = delta * delta;
      const double weight = 1.0;//(exp((-(i*i+j*j))/(2*squareDelta)))/(2*squareDelta*M_PI);
      const double cX = gradientX(px, py, 0, 0);
      const double cY = gradientY(px, py, 0, 0);

      matrix[0][0] += weight * (cX * cX);
      matrix[1][0] += weight * (cX * cY);
      matrix[0][1] += weight * (cX * cY);
      matrix[1][1] += weight * (cY * cY);
    }
  }

  return matrix;
}

/**
  Calculates the vector (...)

  @param gradientX The horizontal gradient for each pixel of an image
  @param gradientY The vertical gradient for each pixel of an image
  @param timeGradient The time gradient for each pixel of an image
  @param x The x coordinate of the pixel
  @param y The y coordinate of the pixel
  @param offset The distance from the pixel to its furthest horizontal (and vertical) neighbours

  @return The vector (...)
*/
double* calculateVector(CImg<double> gradientX, CImg<double> gradientY, CImg<double> timeGradient, int x, int y, int offset)
{
  const int vectorSize = 2;
  double* vector;
  vector = new double[vectorSize];
  for (int i = 0; i < vectorSize; i++)
  {
    vector[i] = 0;    
  }

  for(int i = -offset; i <= offset; i++)
  {
    const int py = y + i;
    for(int j = -offset; j <= offset; j++)
    {
      const int px = x + j;
      const double cX = gradientX(px, py, 0, 0);
      const double cY = gradientY(px, py, 0, 0);
      const double cT = timeGradient(px, py, 0, 0);

      vector[0] += (cX * cT);
      vector[1] += (cY * cT);
    }
  }

  vector[0] = -vector[0];
  vector[1] = -vector[1];

  return vector;
}

/**
  Multiplies a 2x2 matrix by a 2-dimensional vector

  @param matrix The matrix
  @param vector The bidimensional vector

  @return The result of the multiplication
*/
Vector2 multiplyVectorByMatrix(double** matrix, double* vector)
{
  Vector2 result;

  result.x = matrix[0][0]*vector[0] + matrix[0][1]*vector[1];
  result.y = matrix[1][0]*vector[0] + matrix[1][1]*vector[1];

  return result;
}

/**
  Calculates smallest eigenvalue of a 2x2 matrix

  @param matrix The matrix which eigenvalues will be calculated

  @return The smallest eigenvalue of the input matrix
*/
double minimunEigenValue(double** matrix)
{
  const double trace = matrix[0][0] + matrix[1][1];
  const double determinant = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);

  const double delta = (trace * trace) - (4.0 * determinant);
  const double lambda1 = (trace + sqrt(delta))/2.0;
  const double lambda2 = (trace - sqrt(delta))/2.0;

  if(lambda1 > 0 && lambda2 > 0) return min(lambda1, lambda2);
  else return 0;
}

/**
  Calculates the horizontal and vertical gradientes of an image

  @param image The image which the gradients will be calculated

  @return The imagens containing the gradients of the image
*/
vector< CImg<double> > makeGradients(CImg<double> image)
{
  const int WIDTH = image.width(), HEIGHT = image.height(); 
  CImg<double> resultX(WIDTH, HEIGHT, 1, 1, 0), resultY(WIDTH, HEIGHT, 1, 1, 0);
  vector< CImg<double> > gradients;

   cimg_forXY(image,x,y)
   {
     double c1, c2;
     
     if(x > 0 && x < WIDTH - 1)
     {
       c1 = image(x + 1, y, 0, 0), c2 = image(x - 1, y, 0, 0);

       resultX(x, y, 0, 0) = (c1 - c2)/2.0;
     }
     else
     {
       resultX(x,y,0,0) = 0;
     }

     if(y > 0 && y < HEIGHT - 1)
     {
       c1 = image(x, y + 1, 0, 0), c2 = image(x, y - 1, 0, 0);

       resultY(x, y, 0, 0) = (c1 - c2)/2.0;
     }
     else
     {
       resultY(x,y,0,0) = 0;
     }  
   }

  gradients.push_back(resultX);
  gradients.push_back(resultY);

  return gradients;
}

/**
  Calculates the time gradient of an image

  @param image1 The image which gradient will be calculated
  @param image2 The next image in the frames sequence

  @return The image containing the time gradient of each pixel
*/
CImg<double> makeTimeGradient(CImg<double> image1, CImg<double> image2)
{
  CImg<double> result = subtract(image2, image1);

  return makeImageGray(result);
}

/**
  Find the points to be tracked by the Lucas Kanade algorithm

  @param imageGradients The horizontal and vertical gradients of an image
  @param trackerFilterSize The size of the windows to filter the points found

  @return The points to be tracked during the algorithm
*/
vector<Point> findTrackPoints(vector< CImg<double> > imageGradients, const int trackerFilterSize)
{
  CImg<double> gradientX = imageGradients[0], gradientY = imageGradients[1];
  const int WIDTH = gradientX.width(), HEIGHT = gradientX.height();
  const int offset = (int) trackerFilterSize/2;
  double maxLambda;
  double **lambdaMatrix;
  vector<Point> pointsToTrack;

  lambdaMatrix = new double*[HEIGHT];
  for (int i = 0; i < WIDTH; i++)
  {
    lambdaMatrix[i] = new double[WIDTH];
  }

  cimg_forXY(gradientX, x, y)
  {
    if((x > (offset - 1)) && (x < (WIDTH - offset)) && (y > (offset - 1)) && (y < (HEIGHT - offset)))
    {
      double** matrix;
      matrix = calculateMatrix(gradientX, gradientY, x, y, offset);

      const double lambda = minimunEigenValue(matrix, x, y);

      lambdaMatrix[y][x] = lambda;

      if(x == offset && y == offset)
      {
        maxLambda = lambda;
      }
      else
      {
        if(lambda > maxLambda) maxLambda = lambda;
      }
    }
  }

  for (int i = offset; i < HEIGHT - offset; i++)
  {
    for (int j = offset; j < WIDTH - offset; j++)
    {
      bool maxInNeighbourhood = true;

      for(int k = -offset; k <= offset; k++)
      {
       int py = i + k;
       for(int l = -offset; l <= offset; l++)
       {
         int px = j + l;

         if( lambdaMatrix[py][px] > lambdaMatrix[i][j])
         {
          maxInNeighbourhood = false;
          break;
         }
       }
       if(!maxInNeighbourhood) break;
      }

      if(lambdaMatrix[i][j] > 0.1*maxLambda && maxInNeighbourhood)
      {
        Point newPointToTrack;

        newPointToTrack.y = i;
        newPointToTrack.x = j;

        pointsToTrack.push_back(newPointToTrack);
      }
    }
  }

  return pointsToTrack;
}

/**
  Find the flow vectors of the input points given

  @param pointsToTrack The points to be tracked during the algorithm
  @param imageGradients The horizontal and vertical gradients of an image
  @param trackerFilterSize The size of the windows to filter the points found

  @return The flow vectors for each point of the input
*/
vector<Vector2> calculateOpticalFlow(vector<Point> pointsToTrack, vector< CImg<double> > imageGradients, CImg<double> timeGradient, const int trackerFilterSize)
{
  const int offset = (int) trackerFilterSize/2;
  
  CImg<double> gradientX = imageGradients[0], gradientY = imageGradients[1];
  vector<Vector2> opticalFlows;
  double averageVx = 0, averageVy = 0;
  
  for(vector<Point>::iterator it = pointsToTrack.begin(); it != pointsToTrack.end(); ++it)
  {
    const double x = it->x;
    const double y = it->y;

    double **matrix, **inverse, *vector;
    matrix = calculateMatrix(gradientX, gradientY, x, y, offset);

    inverse = invertMatrix(matrix);
    vector = calculateVector(gradientX, gradientY, timeGradient, x, y, offset);

    Vector2 opticalFlow = multiplyVectorByMatrix(inverse, vector);

    averageVx += opticalFlow.x;
    averageVy += opticalFlow.y;

    opticalFlows.push_back(opticalFlow);
  }

  cout << averageVx/pointsToTrack.size() << " " << averageVy/pointsToTrack.size() << endl;

  return opticalFlows;
}

void pyramidalOpticalFlow(vector<Point> pointsToTrack, vector< CImg<double> > imageOneGradients, vector< CImg<double> > imageTwoGradients, const int trackerFilterSize)
{
  const int NUMBER_OF_IMAGES = imageOneGradients.size();

  for(int i = NUMBER_OF_IMAGES - 1; i >= 0; i++)
  {

  }
}

/**
  Find the points corresponding to the last level of th pyramid

  @param pointsToTrack The points to be tracked during the algorithm
  @param numberOfLevels The number of levels of the pyramid

  @return The corresponding track points in the pyramid last level
*/
vector<Point> findLastLevelPoints(vector<Point> pointsToTrack, int numberOfLevels)
{
  const int p = numberOfLevels - 1;

  vector<Point> newPoints;

  for(vector<Point>::iterator it = pointsToTrack.begin(); it != pointsToTrack.end(); ++it)
  {
    const double x = it->x;
    const double y = it->y;

    Point newPoint;

    newPoint.x = x/pow(2.0,p);
    newPoint.y = y/pow(2.0,p);

    newPoints.push_back(newPoint);
  }

  return newPoints;
}

int main (int argc, char **argv)
{
   if(argc < 3)
   {
   	 cout << "Insufficient arguments\n";
     cout << "You must type: " << argv[0] << " <first image file path> <second image file path>";
   	 exit(-1);
   }

  const int gaussianFilterSize = 5, trackerFilterSize = 3;
  CImg<double> frame1(argv[1]), frame2(argv[2]);

  frame1 = makeImageGray(frame1);
  frame2 = makeImageGray(frame2);
  const int WIDTH = frame1.width(), HEIGHT = frame2.height();

  vector< CImg<double> > gradientsFrame1 = makeGradients(frame1);
  vector< CImg<double> > gradientsFrame2 = makeGradients(frame2);

  CImg<double> timeGradient = makeTimeGradient(frame1, frame2);

  vector< CImg<double> > gaussianPyramidFrame1 = makeGaussianPyramid(frame1, gaussianFilterSize);
  vector< CImg<double> > gaussianPyramidFrame2 = makeGaussianPyramid(frame2, gaussianFilterSize);

  vector<Point> pointsToTrack = findTrackPoints(gradientsFrame1, trackerFilterSize);

  vector<Point> newPointsToTrack = findLastLevelPoints(pointsToTrack, gaussianPyramidFrame1.size());
  
  pyramidalOpticalFlow(newPointsToTrack, gaussianPyramidFrame1, gaussianPyramidFrame2, trackerFilterSize);

  return 0;
}