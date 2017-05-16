#include <iostream>
#include <math.h>
using namespace std;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

// Image processing and visualization
#include "./engine/EdgeProcessor.h"
#include "./util/Visualization.h"

int main(int argc, const char* argv[])
{
	// Open image
	cv::Mat image;
	if (argc > 1)
	/*
	 * open image at given location
	 */
	{
		image = cv::imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	}
	else
	/*
	 * open defaut image if no location was given
	 */
	{
		image = cv::imread(
				"../testimages/ethz_binary/ETHZ_APPLELOGO_BINARY/crystal.tif",
				CV_LOAD_IMAGE_GRAYSCALE);
	}

	// Check for valid input
	if (!image.data)
	{
		cout << "Could not find or open image." << endl;
		return -1;
	}
	
	// Display original image
	/*
	cv::namedWindow("Original Image", CV_WINDOW_AUTOSIZE);
	cv::imshow("Original Image", image);
	*/

	// Read edges
	EdgeProcessor edgeProcessor;

	//TODO: maybe use unordered set of C++11
	std::vector<Edge*> edgeVector;
	std::vector<InterestPoint*> pointVector;

	// Process image and check required time
	clock_t t1 = clock();
	edgeProcessor.getEdgesAndInterestPoints(image, &edgeVector, &pointVector);
	//edgeProcessor.closeLoops(&edgeVector);
	clock_t t2 = clock();
	std::cout << "getEdgesAndInterestPoints took: " << double(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;

	// Print number of edges and ambiguous points
	std::cout << "Found " << edgeVector.size() << " edges and "
		<< pointVector.size() << " ambiguous points" << std::endl;
		
	std::cout << "Length of first edge is " << edgeVector[0]->length() << "px" << std::endl;
	
    // Visualization (write SVG)
    Visualization visualization;
    //visualization.validateEdges(image, &edgeVector, &pointVector);
    visualization.saveSVG(image, &edgeVector, 1);

	// End of program
	cout << "Goodbye." << endl;
	return 0;
}
