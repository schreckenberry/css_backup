#include <iostream>
#include <math.h>
using namespace std;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "./engine/EdgeProcessor.h"
#include "./engine/EdgeAnalyzer.h"

cv::RNG rng(12345);
cv::Mat edgeImage;
cv::Mat postProcessedEdgeImage;
bool verbose = false;

void callbackFunction(int event, int x, int y, int flags, void* userdata)
{
    if  ( event == cv::EVENT_LBUTTONDOWN )
    {
        //cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if  ( event == cv::EVENT_RBUTTONDOWN )
    {
        cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
        cv::imwrite( "./edgeImage.png", edgeImage );
        cv::imwrite( "./postProcessedEdgeImage.png", postProcessedEdgeImage );
        cout << "Images saved" << endl;
    }
    else if  ( event == cv::EVENT_MBUTTONDOWN )
    {
        //cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if ( event == cv::EVENT_MOUSEMOVE )
    {
        //cout << "Mouse move over the window - position (" << x << ", " << y << ")" << endl;
    }
}

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
	std::cout << "getEdgesAndInterestPoints took: " << elapsed_secs << std::endl;
/*
	//Hu Moments
	double* hu;
	hu = edgeVector[0]->getHuMoments(5, 10.0);

	for(int i = 0; i < 7; i++)
	{
		std::cout << hu[i] << std::endl;
	}

*/
	if(verbose){
		// print edges
		std::cout << "---------RESULT---------" << std::endl;
		for (int i = 0; i < edgeVector.size(); i++)
		{
			std::cout << "Edge index: " << i << std::endl;
			std::cout << "Edge start direction: " << edgeVector[i]->getStartDirection(5) << " rad" << std::endl;
			std::cout << "Edge end direction: " << edgeVector[i]->getEndDirection(5) << " rad" << std::endl;

			for (int j = 0; j < edgeVector[i]->getEdgePoints().size(); j++)
			{
				std::cout << "[" << edgeVector[i]->getEdgePoints()[j].x << ", " << edgeVector[i]->getEdgePoints()[j].y << "]";
			}
			std::cout << std::endl;
			std::cout << "Number of following Edges: " << edgeVector[i]->getNumberFollowingEdges() << std::endl;
		}
	}

	std::cout << "Total number of Edges found: " << edgeVector.size() << std::endl;
	std::cout << "Total number of InterestPoints found: " << pointVector.size() << std::endl;


	// draw edges
	edgeImage = cv::Mat::zeros(image.size(), CV_8UC3);
	for (int i = 0; i < edgeVector.size(); i++)
	{
		//cv::Scalar color = cv::Scalar(0, 0, 255); //BGR
		cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));

		for (int j = 0; j < edgeVector.at(i)->getEdgePoints().size(); j++)
		{
			cv::circle(edgeImage, edgeVector.at(i)->getEdgePoints()[j], 0, color, 1, 8, 0);
		}
	}


	// draw points
	for (int k = 0; k < pointVector.size(); k++)
	{
		cv::Scalar color = cv::Scalar(0, 0, 255); //BGR
		cv::line(edgeImage,
		         cv::Point(pointVector.at(k)->getPoint().x - 2,
		         pointVector.at(k)->getPoint().y),
		         cv::Point(pointVector.at(k)->getPoint().x + 2,
			     pointVector.at(k)->getPoint().y), color, 1, 8, 0);
			     cv::line(edgeImage,
			     cv::Point(pointVector.at(k)->getPoint().x,
			     pointVector.at(k)->getPoint().y - 2),
			     cv::Point(pointVector.at(k)->getPoint().x,
			     pointVector.at(k)->getPoint().y + 2), color, 1, 8, 0);

		//cv::circle(postProcessedEdge, pointVector.at(k)->getPoint(), 0, color, 2, 8, 0);
	}


	// postprocess edges
	int T = 20;
	int lengthForDirection = 5;
	std::cout << "Postprocessing started " << std::endl;
	edgeProcessor.postProcessEdges(&edgeVector, T, lengthForDirection, M_PI);
	std::cout << "Postprocessing finished " << std::endl;
	if(verbose){
	std::cout << "---------RESULT---------" << std::endl;
		for (int i = 0; i < edgeVector.size(); i++)
		{
			std::cout << "Edge index: " << i << std::endl;
			std::cout << "Edge start direction: " << edgeVector[i]->getStartDirection(5) << " rad" << std::endl;
			std::cout << "Edge end direction: " << edgeVector[i]->getEndDirection(5) << " rad" << std::endl;

			for (int j = 0; j < edgeVector[i]->getEdgePoints().size(); j++)
			{
				std::cout << "[" << edgeVector[i]->getEdgePoints()[j].x << ", " << edgeVector[i]->getEdgePoints()[j].y << "]";
			}
			std::cout << std::endl;
			std::cout << "Number of Edge Followers: " << edgeVector[i]->getNumberFollowingEdges() << std::endl;
		}
	}

	std::cout << "Total number of Edges found: " << edgeVector.size() << std::endl;
	std::cout << "Total number of InterestPoints found: " << pointVector.size() << std::endl;

	// 3&4. analyze edges
	// EdgeAnalyzer analyzeEdges(vector<Edge>) returns vector<InterestPoint>
	EdgeAnalyzer edgeAnalyzer;
	std::vector<InterestPoint*> maximaVector;
	std::vector<InterestPoint*> minimaVector;
	std::vector<InterestPoint*> zeroVector;
	clock_t t3 = clock();
	edgeAnalyzer.calculateScaleSpaces(&edgeVector, &maximaVector, &minimaVector, &zeroVector);
	clock_t t4 = clock();
	std::cout << "calculateScaleSpaces took: " << double(t4 - t3) / CLOCKS_PER_SEC << std::endl;


	// draw edges
	postProcessedEdgeImage = cv::Mat::zeros(image.size(), CV_8UC3);
	for (int i = 0; i < edgeVector.size(); i++)
	{
		//cv::Scalar color = cv::Scalar(0, 0, 255); //BGR
		cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));

		for (int j = 0; j < edgeVector.at(i)->getEdgePoints().size(); j++)
		{
			cv::circle(postProcessedEdgeImage, edgeVector.at(i)->getEdgePoints()[j], 0, color, 1, 8, 0);
		}
	}

	std::cout << "edgeVector.at(0).size() = " << edgeVector.at(0)->length() << std::endl;

	// draw points
	for (int k = 0; k < pointVector.size(); k++)
	{
		cv::Scalar color = cv::Scalar(0, 0, 255); //BGR

		cv::circle(postProcessedEdgeImage, pointVector.at(k)->getPoint(), 0, color, 0, 8, 0);
	}
	for (int k = 0; k < maximaVector.size(); k++)
	{
		cv::Scalar color = cv::Scalar(0, 255, 0); //BGR

		cv::circle(postProcessedEdgeImage, maximaVector.at(k)->getPoint(), 2, color, 0, 8, 0);
	}
	for (int k = 0; k < minimaVector.size(); k++)
	{
		cv::Scalar color = cv::Scalar(0, 0, 255); //BGR

		cv::circle(postProcessedEdgeImage, minimaVector.at(k)->getPoint(), 2, color, 0, 8, 0);
	}
	for (int k = 0; k < zeroVector.size(); k++)
	{
		cv::Scalar color = cv::Scalar(255, 255, 255); //BGR

		cv::circle(postProcessedEdgeImage, zeroVector.at(k)->getPoint(), 2, color, 0, 8, 0);
	}

	// 5. match edges
	// EdgeMatcher matchEdges(vector<InterestPoint> obj, vector<InterestPoint> image) returns vector<position>

	// Display
	cv::namedWindow("Image display", CV_WINDOW_AUTOSIZE);
	cv::imshow("Image display", image);

	// Display
	cv::namedWindow("Edges and InterestPoints", CV_WINDOW_AUTOSIZE);
	//set the callback function for any mouse event
	cv::setMouseCallback("Edges and InterestPoints", callbackFunction, NULL);
	imshow("Edges and InterestPoints", edgeImage);

	// Display
	cv::namedWindow("Post Processed Edges and InterestPoints", CV_WINDOW_AUTOSIZE);
	//set the callback function for any mouse event
	cv::setMouseCallback("Post Processed Edges and InterestPoints", callbackFunction, NULL);
	imshow("Post Processed Edges and InterestPoints", postProcessedEdgeImage);

	cv::waitKey(0);

	cout << "Goodbye." << endl;
	return 0;
}
