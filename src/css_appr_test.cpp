#include <iostream>
#include <math.h>
using namespace std;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/opencv.hpp"

#include "./engine/EdgeProcessor.h"
#include "./engine/CurvatureScaleSpace.h"
#include "./util/Visualization.h"

// Approximation
#include "./engine/ApprCSS.h"

int main(int argc, const char* argv[])
{
	// Open image
	cv::Mat image;
	if (argc > 1)
		image = cv::imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE); 

	// Check for valid input
	if (!image.data)
	{
		std::cout << "Could not find or open image." << endl;
		return -1;
	}
	
	// Read edges
	EdgeProcessor edgeProcessor;
	std::vector<Edge*> edgeVector;
	std::vector<InterestPoint*> pointVector;
	edgeProcessor.getEdgesAndInterestPoints(image, &edgeVector, &pointVector);

	// print number of edges and ambiguous points
	std::cout << "Found " << edgeVector.size() << " edges and " << pointVector.size() << " ambiguous points" << std::endl;
	
	// Curvature Scale Space Calculation and Tests (Max)
	std::cout << "************************************************************************" << std::endl;
	clock_t t3 = clock();
	CSSApproximation cssappr;
	cssappr.calculateCSSAppr(&edgeVector);
	clock_t t4 = clock();
	std::cout << "calculateScaleSpaces (Approximation) took: " << double(t4 - t3) / CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << "************************************************************************" << std::endl;
	
	cout << "Goodbye." << endl;
	return 0;
}
