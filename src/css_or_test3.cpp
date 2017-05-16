#include <iostream>
#include <math.h>
using namespace std;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "./engine/EdgeProcessor.h"
#include "./engine/CurvatureScaleSpace.h"
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

	// Check time for processing
	clock_t t1 = clock();
	edgeProcessor.getEdgesAndInterestPoints(image, &edgeVector, &pointVector);
	clock_t t2 = clock();
	std::cout << "getEdgesAndInterestPoints took: " << double(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;

	// print number of edges and ambiguous points
	std::cout << "Found " << edgeVector.size() << " edges and "
		<< pointVector.size() << " ambiguous points" << std::endl;
	
    // Visualization
    Visualization visualization;
    visualization.validateEdges(image, &edgeVector, &pointVector);
    visualization.saveSVG(image, &edgeVector, 1);
    
    // Curvature Scale Space calculation
	CurvatureScaleSpace css;
	clock_t t3 = clock();
	css.calculateScaleSpaces(&edgeVector);
	clock_t t4 = clock();
	std::cout << "calculateScaleSpaces took: " << double(t4 - t3) / CLOCKS_PER_SEC << "s" << std::endl;
	
	// draw and display resulting image
	cv::Mat imageResult(image.rows, image.cols, CV_8UC3, 0.0);
	for (int i = 0; i < edgeVector.size(); i++)
	{
		visualization.generateColorWheel(edgeVector.at(i)->getEdgePoints().size());
		for (int j = 0; j < edgeVector.at(i)->getEdgePoints().size(); j++)
		{
			//cv::circle(imageResult, edgeVector.at(i)->getEdgePoints()[j], 0, color, 1, 8, 0);
			cv::circle(imageResult, edgeVector.at(i)->getEdgePoints()[j], 0, visualization.getColorWheelValueBGR(j), 1, 8, 0);
			
			// mark first element
			if (j == 0)
			{
				cv::circle(imageResult, edgeVector.at(i)->getEdgePoints()[j], 0, cv::Scalar(255, 255, 255), 1, 8, 0);
			}
		}
	}

	cv::namedWindow("Output Image", CV_WINDOW_AUTOSIZE);
	cv::imshow("Output Image", imageResult);
	cv::imwrite("./postProcessedEdgeImage3.png", imageResult);
	
	cv::waitKey(0);
	cout << "Goodbye." << endl;
	return 0;
}
