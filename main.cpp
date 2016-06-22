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
  Calculates the spatial gradient matrix

  @param gradientX The horizontal gradient for each pixel of an image
  @param gradientY The vertical gradient for each pixel of an image
  @param x The x coordinate of the pixel
  @param y The y coordinate of the pixel
  @param offset The distance from the pixel to its furthest horizontal (and vertical) neighbours

  @return The spatial gradient matrix
*/
double** calculateMatrix(CImg<double> gradientX, CImg<double> gradientY, int x, int y, int offset)
{
  
  const int matrixSize = 2;
  const int HEIGHT = gradientX.height(), WIDTH = gradientX.width();
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
//assert(px >= 0); assert(px <= WIDTH - 1);
//assert(py >= 0); assert(py <= HEIGHT - 1);
      //const double delta = sqrt(2.0);
      //const double squareDelta = delta * delta;
      const double weight = 1.0;//(exp((-(i*i+j*j))/(2*squareDelta)))/(2*squareDelta*M_PI);
      double cX, cY;
      
      if(px >= 0 && px <= (WIDTH - 1) && py >= 0 && py <= (HEIGHT - 1))
      {
        cX = gradientX(px, py, 0, 0);
        cY = gradientY(px, py, 0, 0);
      }
      else 
      {
        cX = 0.0;
        cY = 0.0;
      }

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
  Find the points corresponding to the given level of the pyramid
  @param pointsToTrack The points to be tracked during the algorithm
  @param levels The pyramid level index (lowest is 0)
  @return The corresponding track points in the pyramid last level
*/
vector<Point> findCorrespondingLevelPoints(vector<Point> pointsToTrack, int level)
{
  vector<Point> newPoints;

  for(vector<Point>::iterator it = pointsToTrack.begin(); it != pointsToTrack.end(); ++it)
  {
    const double x = it->x;
    const double y = it->y;

    Point newPoint;

    newPoint.x = x/pow(2.0,level);
    newPoint.y = y/pow(2.0,level);

    newPoints.push_back(newPoint);
  }

  return newPoints;
}

/**
  Find the flows corresponding to the given level of the pyramid
  @param flows The flows found during the algorithm
  @param levels The pyramid level index (lowest is 0)
  @return The corresponding track points in the pyramid last level
*/
vector<Vector2> findCorrespondingLevelFlows(vector<Vector2> flows, int level)
{
  vector<Vector2> newFlows;

  for(vector<Vector2>::iterator it = flows.begin(); it != flows.end(); ++it)
  {
    const double x = it->x;
    const double y = it->y;

    Vector2 newFlow;

    newFlow.x = x/pow(2.0,level);
    newFlow.y = y/pow(2.0,level);

    newFlows.push_back(newFlow);
  }

  return newFlows;
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
CImg<double> makeTimeGradient(CImg<double> image1, CImg<double> image2, Point pointToTrack, Vector2 lastLevelFlow, Vector2 iterationGuess, const int offset)
{
  const int WIDTH = image1.width(), HEIGHT = image1.height();
  const int WIDTH2 = image2.width(), HEIGHT2 = image2.height(); 
  CImg<double> result(WIDTH, HEIGHT, 1, 1, 0);

  cimg_forXY(result, x, y)
  {
    result(x,y,0,0) = 0;
  }
  
  int x = pointToTrack.x, y = pointToTrack.y;

  for(int i = -offset; i <= offset; i++)
  {
    const int py = y + i;
    for(int j = -offset; j <= offset; j++)
    {
      const int px = x + j;
      if(px >= 0 && px <= (WIDTH - 1) && py >= 0 && py <= (HEIGHT - 1))
      {
        int px2 = px + lastLevelFlow.x + iterationGuess.x;
        int py2 = py + lastLevelFlow.y + iterationGuess.y;

        if(px2 < 0) px2 = 0;
        else if(px2 >= WIDTH) px2 = WIDTH - 1;
        
        if(py2 < 0) py2 = 0;
        else if(py2 >= HEIGHT) py2 = HEIGHT - 1;

        result(px,py,0,0) = abs(image1(px, py, 0, 0) - image2(px2, py2, 0, 0));
      }
    }
  }

   return result;  
}

/**
  Find the points to be tracked by the Lucas Kanade algorithm

  @param imageGradients The horizontal and vertical gradients of an image
  @param trackerFilterSize The size of the windows to filter the points found

  @return The points to be tracked during the algorithm
*/
vector<Point> findTrackPoints(vector< CImg<double> > imageGradients, const int trackerFilterSize, Point firstBoundaryPoint, Point secondBoundaryPoint)
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
    for (int j = 0; j < WIDTH; j++) lambdaMatrix[i][j] = 0;
  }
  
  const int minX = min(firstBoundaryPoint.x, secondBoundaryPoint.x);
  const int minY = min(firstBoundaryPoint.y, secondBoundaryPoint.y);
  const int maxX = max(firstBoundaryPoint.x, secondBoundaryPoint.x);
  const int maxY = max(firstBoundaryPoint.y, secondBoundaryPoint.y);

  for(int y = minY + 1; y < maxY; y++)
  {
    for(int x = minX + 1; x < maxX; x++)
    {
      if((x > (offset - 1)) && (x < (WIDTH - offset)) && (y > (offset - 1)) && (y < (HEIGHT - offset)))
      {
        double** matrix;
        matrix = calculateMatrix(gradientX, gradientY, x, y, offset);

        const double lambda = minimunEigenValue(matrix);

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
  }

  for (int i = minY + offset + 1; i < maxY - offset - 1; i++)
  {
    for (int j = minX + offset + 1; j < maxX - offset - 1; j++)
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
  Find the flow vector of the input point given

  @param pointToTrack The point to be tracked during the algorithm
  @param imageGradients The horizontal and vertical gradients of an image
  @param trackerFilterSize The size of the windows to filter the points found

  @return The flow vectors for each point of the input
*/
Vector2 calculateOpticalFlow(Point pointToTrack, Vector2 lastLevelFlow, vector< CImg<double> > imageGradients, CImg<double> frame1, CImg<double> frame2, const int offset)
{
  const int NUMBER_OF_ITERATIONS = 1;

  Vector2 opticalFlow;
  Vector2 iterationGuess;
  iterationGuess.x = 0.0;
  iterationGuess.y = 0.0;
  
  CImg<double> gradientX = imageGradients[0], gradientY = imageGradients[1];
  
  const double x = pointToTrack.x;
  const double y = pointToTrack.y;
  double **matrix, **inverse, *vector;

  matrix = calculateMatrix(gradientX, gradientY, ceil(x), ceil(y), offset);

  inverse = invertMatrix(matrix);

  for(int i = 1; i <= NUMBER_OF_ITERATIONS; i++)
  {
    CImg<double> timeGradient = makeTimeGradient(frame1, frame2, pointToTrack, lastLevelFlow, iterationGuess, offset);
    vector = calculateVector(gradientX, gradientY, timeGradient, ceil(x), ceil(y), offset);

    opticalFlow = multiplyVectorByMatrix(inverse, vector);

    iterationGuess.x = iterationGuess.x + opticalFlow.x;
    iterationGuess.y = iterationGuess.y + opticalFlow.y;

    if(isnan(opticalFlow.x)) opticalFlow.x = 0.0;
    if(isnan(opticalFlow.y)) opticalFlow.y = 0.0;
  }

  opticalFlow = iterationGuess;

  return opticalFlow;
}

/**
  Find the flow vectors of the input points given

  @param pointsToTrack The points to be tracked during the algorithm
  @param imageGradients The horizontal and vertical gradients of an image
  @param trackerFilterSize The size of the windows to filter the points found

  @return The flow vectors for each point of the input
*/
vector<Vector2> calculateOpticalFlow(vector<Point> pointsToTrack, vector< CImg<double> > imageGradients, vector<Vector2> lastLevelFlows, CImg<double> frame1, CImg<double> frame2, const int trackerFilterSize)
{
  const int offset = (int) trackerFilterSize/2;

  vector<Vector2> opticalFlows;
  double averageVx = 0, averageVy = 0;

  for(int j = 0; j < pointsToTrack.size(); j++)
  {
    Point point = pointsToTrack[j];
    Vector2 lastLevelFlow = lastLevelFlows[j];
    Vector2 opticalFlow = calculateOpticalFlow(point, lastLevelFlow, imageGradients, frame1, frame2, offset);

    averageVx += opticalFlow.x;
    averageVy += opticalFlow.y;

    opticalFlows.push_back(opticalFlow);
  }

  //cout << averageVx/pointsToTrack.size() << " " << averageVy/pointsToTrack.size() << endl;

  return opticalFlows;
}

vector<Vector2> pyramidalOpticalFlow(vector<Point> points, vector< CImg<double> > imageOnePyramid, vector< CImg<double> > imageTwoPyramid, const int gaussianFilterSize, const int trackerFilterSize)
{
  const int NUMBER_OF_IMAGES = imageOnePyramid.size();
  
  vector<Vector2> opticalFlows, aproximations;
  vector<Point> pointsToTrack;

  for(int i = 0; i < points.size(); i++)
  {
    Vector2 v;
    v.x = 0;
    v.y = 0;
    aproximations.push_back(v);
  }

  pointsToTrack = findCorrespondingLevelPoints(points, NUMBER_OF_IMAGES - 1);
  for(int i = NUMBER_OF_IMAGES - 1; i >= 0; i--)
  {
    cout << "Calculating flow for level " << i << endl;
    CImg<double> frame1 = imageOnePyramid[i], frame2 = imageTwoPyramid[i];

    cout << "Making gradients ";
    vector< CImg<double> > gradientsFrame1 = makeGradients(frame1);
    cout << "Finished\n";

    cout << "Finding flows ";
    opticalFlows =  calculateOpticalFlow(pointsToTrack, gradientsFrame1, aproximations, frame1, frame2, trackerFilterSize);

    for(int j = 0; j < pointsToTrack.size(); j++)
    {
      pointsToTrack[j].x = 2.0*(pointsToTrack[j].x + opticalFlows[j].x);
      pointsToTrack[j].y = 2.0*(pointsToTrack[j].y + opticalFlows[j].y);

      aproximations[j].x = 2.0*(aproximations[j].x + opticalFlows[j].x);
      aproximations[j].y = 2.0*(aproximations[j].y + opticalFlows[j].y);
    }
    cout << "Finished\n";
  }

  opticalFlows = aproximations;
  double averageVx = 0.0, averageVy = 0.0;

  for(int i = 0; i < opticalFlows.size(); i++)
  {
    averageVx += opticalFlows[i].x;
    averageVy += opticalFlows[i].y;
  }
  cout << averageVx/opticalFlows.size() << " " << averageVy/opticalFlows.size() << endl;

  return opticalFlows;
}

int main (int argc, char **argv)
{
   if(argc < 3)
   {
   	 cout << "Insufficient arguments\n";
     cout << "You must type: " << argv[0] << " <number of images> <images format> <files names prefix> (Optional)";
   	 exit(-1);
   }

  const int gaussianFilterSize = 5, trackerFilterSize = 3;
  const int NUMBER_OF_IMAGES = atoi(argv[1]);
  vector< CImg<double> > frames;

  for(int i = 1; i <= NUMBER_OF_IMAGES; i++)
  { 
    std::stringstream str;
    str << i;
    char fileName[100];
    if(argc >= 4)
    {
      strcpy(fileName, argv[3]);
      strcat(fileName, str.str().c_str());
    }
    else strcpy(fileName, str.str().c_str());
    strcat(fileName, ".");
    strcat(fileName, argv[2]);

    CImg<double> frame(fileName);
    frames.push_back(frame);
  }

  CImg<double> frame1 = makeImageGray(frames[0]);
  CImg<double> frame2 = makeImageGray(frames[1]);
  const int WIDTH = frame1.width(), HEIGHT = frame1.height();
  cout << "Making gradients\n";
  clock_t begin_time = clock();
  vector< CImg<double> > gradientsFrame1 = makeGradients(frame1);
  vector< CImg<double> > gradientsFrame2 = makeGradients(frame2);

  cout << "Finished in " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds\n";

  cout << "Making gaussian pyramids\n";
  begin_time = clock();
  vector< CImg<double> > gaussianPyramidFrame1 = makeGaussianPyramid(frame1, gaussianFilterSize);
  vector< CImg<double> > gaussianPyramidFrame2 = makeGaussianPyramid(frame2, gaussianFilterSize);
  cout << "Finished in " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds\n";

  vector< Point > boundary = getBoundary(frame1);

  Point firstBoundaryPoint = boundary[0], secondBoundaryPoint = boundary[1];

  cout << "Finding first points\n";
  begin_time = clock();
  vector<Point> pointsToTrack = findTrackPoints(gradientsFrame1, trackerFilterSize, firstBoundaryPoint, secondBoundaryPoint);
  cout << "Finished in " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds\n";

  vector<Vector2> opticalFlows = pyramidalOpticalFlow(pointsToTrack,  gaussianPyramidFrame1, gaussianPyramidFrame2, gaussianFilterSize, trackerFilterSize);
  
  markPointsOnImage(frame1, pointsToTrack);
  int x1 = -1, y1 = -1, x0 = WIDTH, y0 = HEIGHT;

  for(int j = 0; j < pointsToTrack.size(); j++)
  {      
    if(pointsToTrack[j].x >= 0 && pointsToTrack[j].x < WIDTH)
    {
      if(pointsToTrack[j].x < x0) x0 = pointsToTrack[j].x;
      if(pointsToTrack[j].x > x1) x1 = pointsToTrack[j].x;
    }

    if(pointsToTrack[j].y >= 0 && pointsToTrack[j].y < HEIGHT)
    {
      if(pointsToTrack[j].y < y0) y0 = pointsToTrack[j].y;
      if(pointsToTrack[j].y > y1) y1 = pointsToTrack[j].y;
    }
  }
  displayImageWithRectangle(frame1, x0, y0, x1, y1);
  
  for(int i = 1; i < NUMBER_OF_IMAGES - 1; i++)
  {
    cout << "Frame " << i - 1 << " to " << i << endl;
    int x1 = -1, y1 = -1, x0 = WIDTH, y0 = HEIGHT;
    for(int j = 0; j < pointsToTrack.size(); j++)
    {      
      pointsToTrack[j].x += opticalFlows[j].x;
      pointsToTrack[j].y += opticalFlows[j].y;

      if(pointsToTrack[j].x >= 0 && pointsToTrack[j].x < WIDTH)
      {
        if(pointsToTrack[j].x < x0) x0 = pointsToTrack[j].x;
        if(pointsToTrack[j].x > x1) x1 = pointsToTrack[j].x;
      }

      if(pointsToTrack[j].y >= 0 && pointsToTrack[j].y < HEIGHT)
      {
        if(pointsToTrack[j].y < y0) y0 = pointsToTrack[j].y;
        if(pointsToTrack[j].y > y1) y1 = pointsToTrack[j].y;
      }
    }

    frame1 = frame2;
    frame2 = makeImageGray(frames[i + 1]);
    
    displayImageWithRectangle(frame1, x0, y0, x1, y1);
    gradientsFrame1 = gradientsFrame2;
    gradientsFrame2 = makeGradients(frame2);

    gaussianPyramidFrame1 = gaussianPyramidFrame2;
    gaussianPyramidFrame2 = makeGaussianPyramid(frame2, gaussianFilterSize);

    opticalFlows = pyramidalOpticalFlow(pointsToTrack,  gaussianPyramidFrame1, gaussianPyramidFrame2, gaussianFilterSize, trackerFilterSize);
  }
  
  x1 = -1; y1 = -1; x0 = WIDTH; y0 = HEIGHT;

  for(int j = 0; j < pointsToTrack.size(); j++)
  {      
    pointsToTrack[j].x += opticalFlows[j].x;
    pointsToTrack[j].y += opticalFlows[j].y;

    if(pointsToTrack[j].x >= 0 && pointsToTrack[j].x < WIDTH)
    {
      if(pointsToTrack[j].x < x0) x0 = pointsToTrack[j].x;
      if(pointsToTrack[j].x > x1) x1 = pointsToTrack[j].x;
    }

    if(pointsToTrack[j].y >= 0 && pointsToTrack[j].y < HEIGHT)
    {
      if(pointsToTrack[j].y < y0) y0 = pointsToTrack[j].y;
      if(pointsToTrack[j].y > y1) y1 = pointsToTrack[j].y;
    }
  }

  //markPointsOnImage(frame2, pointsToTrack);
   displayImageWithRectangle(frame1, x0, y0, x1, y1);

  return 0;
}