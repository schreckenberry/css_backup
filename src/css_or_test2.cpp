#include <iostream>
#include <math.h>
using namespace std;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "./engine/EdgeProcessor.h"

cv::RNG rng(12345);

int main(int argc, const char* argv[])
{
	// 1. Open image
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
	cv::namedWindow("Original_image", CV_WINDOW_AUTOSIZE);
	cv::imshow("Original_image", image);

	// Program structure
	// 2. Read edges
	EdgeProcessor edgeProcessor;

	//TODO maybe use unordered set of C++11
	std::vector<Edge*> edgeVector;
	std::vector<InterestPoint*> pointVector;

	// Check time for processing
	clock_t t1 = clock();
	edgeProcessor.getEdgesAndInterestPoints(image, &edgeVector, &pointVector);
	clock_t t2 = clock();
	double elapsed_secs = double(t2 - t1) / CLOCKS_PER_SEC;
	std::cout << "getEdgesAndInterestPoints took: " << elapsed_secs << "s" << std::endl;
	
	// ================================================================================================
	// Draw edges (edgeVector) without postprocessing and count shared pixels
	cv::Mat unPostProcessedEdgeImage = cv::Mat::zeros(image.size(), CV_8UC3);
	cv::Mat sharedPixels(image.size(), CV_8UC1, 0.0);

	for (int i = 0; i < edgeVector.size(); i++)
	{
		//cv::Scalar color = cv::Scalar(0, 0, 255); //BGR
		cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));

		for (int j = 0; j < edgeVector.at(i)->getEdgePoints().size(); j++)
		{
			// draw single pixels of edge i consisting of j points
			cv::circle(unPostProcessedEdgeImage, edgeVector.at(i)->getEdgePoints()[j], 0, color, 1, 8, 0);
			// count shared pixels
	        sharedPixels.at<uchar>(edgeVector.at(i)->getEdgePoints()[j])++;
		}
	}
	
	// find shared pixels and mark in unPostProcessedEdgeImage
	int numberSharedPixels = 0;
	for (int y = 0; y < sharedPixels.rows; y++)
	{
		for (int x = 0; x < sharedPixels.cols; x++)
		{
			// pixel is shared
			if (sharedPixels.at<uchar>(y, x) > 1.0)
			{
				numberSharedPixels++;
				cout << "#Edges at (" << y << ", " << x << ") = " << (int) sharedPixels.at<uchar>(y, x) << endl;
				cv::circle(unPostProcessedEdgeImage, cv::Point(x, y), 0, cv::Scalar(0, 0, 255), 2, 8, 0);
			}
		}
	}
	cout << "numberSharedPixels = " << numberSharedPixels << endl;
	//cv::imwrite("/home/markus01/Desktop/output.png", unPostProcessedEdgeImage);
	
	// Display drawing
	cv::namedWindow("unPostProcessedEdgeImage", CV_WINDOW_AUTOSIZE);
	imshow("unPostProcessedEdgeImage", unPostProcessedEdgeImage);
	// ================================================================================================

	cv::waitKey(0);
	cout << "Goodbye." << endl;
	return 0;
}
