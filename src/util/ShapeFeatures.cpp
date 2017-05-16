#include "ShapeFeatures.h"

#include "Visualization.h"
Visualization vis;

int testFrag = -1; // set -1 to avoid test

ShapeFeatures::ShapeFeatures()
{
	//std::cout << "ShapeFeatures object created" << std::endl;
}

void ShapeFeatures::fillFragmentArea(cv::Mat &fragArea, std::vector<cv::Point> &fragPts)
{
	// determine top left corner
	int xMin = 100000; int xMax = -100000;
	int yMin = 100000; int yMax = -100000;
	
	for (size_t i = 0; i < fragPts.size(); i++)
	{
		// xMin
		if (fragPts[i].x < xMin)
		{
			xMin = fragPts[i].x;
		}
		// yMin
		if (fragPts[i].y < yMin)
		{
			yMin = fragPts[i].y;
		}
		// xMax
		if (fragPts[i].x > xMax)
		{
			xMax = fragPts[i].x;
		}		
		// yMax
		if (fragPts[i].y > yMax)
		{
			yMax = fragPts[i].y;
		}
	}
	
	// align fragment to origin on top left
	for (size_t i = 0; i < fragPts.size(); i++)
	{
		fragPts[i].x = fragPts[i].x - xMin;
		fragPts[i].y = fragPts[i].y - yMin;
	}
	
	// create new image where filled area can be drawn
	int cols = xMax - xMin + 1;
	int rows = yMax - yMin + 1;
	
	fragArea = cv::Mat::zeros(rows, cols, CV_8UC1);

	const int npts = static_cast<int> (fragPts.size());
	const cv::Point* pts = &fragPts[0];
	
	// fill area where fragPts (pts) define vertices of polygon
	cv::fillPoly(fragArea, &pts, &npts, 1, cv::Scalar(255), 8);
	
	// visualize result (nuFragment to choose fragment to be visualized)
	static int nuFragment = -1; nuFragment++;
	if (nuFragment == testFrag)
	{
		vis.saveCvDataSVG(fragArea, fragPts);
	}
}

void ShapeFeatures::fragAreaToVector(cv::Mat &fragArea, std::vector<cv::Point> &fragAreaPts)
{
	fragAreaPts.clear();
	
	// check every pixel and store area in vector
	for (size_t x = 0; x < fragArea.cols; x++)
	{
		for (size_t y = 0; y < fragArea.rows; y++)
		{
			if (fragArea.at<uchar>(y, x) > 0)
			{
				fragAreaPts.push_back(cv::Point(x, y));
			}
		}
	}
}

float ShapeFeatures::correctRotation(float rot, std::vector<cv::Point> &fragPts, int uCoordFrag, float radius, cv::Point2f mean, bool flat)
{
	// determine position of local extremum
	cv::Point extPnt(fragPts[uCoordFrag]);
	
	float angleBase;
	std::vector<float> angles;
	if (flat)
	{
		// determine 2 possible angles
		angles.push_back(rot);
		angles.push_back(rot + CV_PI);

		angleBase = CV_PI;
	}
	else
	{
		// determine 4 possible angles
		angles.push_back(rot);
		angles.push_back(rot + CV_PI/2.0);
		angles.push_back(rot + CV_PI);
		angles.push_back(rot + 3.0*CV_PI/2.0);

		angleBase = CV_PI/2.0;
	}
	
	// determine angle of vector from extPnt to pcaMean
	float refAngle = atan2((mean.y-extPnt.y), (mean.x-extPnt.x));
	
	// determine which angles have the smallest difference
	float minDiff = 100000; int iMin = 0;
	for (size_t i = 0; i < angles.size(); i++)
	{
		float diff = fabs(refAngle-angles[i]);
		
		// take circular distance
		if (diff > CV_PI)
		{
			diff = fabs(2.0*CV_PI - diff);
		}
		
		if (diff <= minDiff)
		{
			minDiff = diff;
			iMin = i;
		}
	}

	static int nuFragment = -1; nuFragment++;
	if (nuFragment == testFrag)
	{
		std::cout << "(extPnt.x, extPnt.y) = (" << extPnt.x << ", " << extPnt.y << ")" << std::endl;
		std::cout << "(mean.x, mean.y) = (" << mean.x << ", " << mean.y << ")" << std::endl;
		std::cout << "radius = " << radius << std::endl;
	}
	
	return (rot + iMin*angleBase);
}

cv::Point2f ShapeFeatures::calcAreaMean(std::vector<cv::Point> &fragPts, int nuAreaPts, cv::Point2f pcaMean)
{
	// calculate mean of the fragment points
	int nuFragPts = (int) fragPts.size();
	float sumX = 0, sumY = 0;
	float meanX, meanY;
	
	for (int i = 0; i < nuFragPts; i++)
	{
		sumX += fragPts[i].x;
		sumY += fragPts[i].y;
	}
	
	// substract mean of pcaMean which is calculated for all points (fragment + area)
	meanX = (pcaMean.x * nuAreaPts - sumX) / (1.0 * (nuAreaPts - nuFragPts));
	meanY = (pcaMean.y * nuAreaPts - sumY) / (1.0 * (nuAreaPts - nuFragPts));
	
	// return result
	return cv::Point2f(meanX, meanY);
}

float ShapeFeatures::calcOrientation(std::vector<cv::Point> &fragPts, int uCoordFrag, float radius, float kappa)
{
	// crop and align fragment and fill fragment area
	cv::Mat fragArea;
	fillFragmentArea(fragArea, fragPts);
	
	// convert area to vector structure
	std::vector<cv::Point> fragAreaPts;
	fragAreaToVector(fragArea, fragAreaPts);

	// save the data in cv::Mat data structure
	int sz = fragAreaPts.size();
	cv::Mat data(sz, 2, CV_64FC1);
	
	for (int i = 0; i < sz; i++)
	{
		data.at<double>(i, 0) = fragAreaPts[i].x;
		data.at<double>(i, 1) = fragAreaPts[i].y;
	}
	
	// Perform PCA Analysis
	cv::PCA pca(data, cv::Mat(), CV_PCA_DATA_AS_ROW);
	
	// store mean value
	float meanX = pca.mean.at<double>(0, 0);
	float meanY = pca.mean.at<double>(0, 1);
	cv::Point2f pcaMean(meanX, meanY);

	// store eigenvalues
	float lambda1 = pca.eigenvalues.at<double>(0, 0);
	float lambda2 = pca.eigenvalues.at<double>(0, 1);
	
	// store eigenvectors
	float v1x = pca.eigenvectors.at<double>(0, 0);
	float v1y = pca.eigenvectors.at<double>(0, 1);
	
	float v2x = pca.eigenvectors.at<double>(1, 0);
	float v2y = pca.eigenvectors.at<double>(1, 1);

	//std::cout << "[v1x, v1y] = " << "[" << v1x << ", " << v1y << "]" << " - lambda1 = " << lambda1 << std::endl;
	//std::cout << "[v2x, v2y] = " << "[" << v2x << ", " << v2y << "]" << " - lambda2 = " << lambda2 << std::endl;
	
	// choose correct axis
	if ((lambda1 > 50 * lambda2) && (fabs(kappa) < 0.02))
	{
		// analyse the two shorter axis
		return correctRotation(atan2(v2y, v2x), fragPts, uCoordFrag, radius, pcaMean, true);
	}
	else if ((lambda2 > 50 * lambda1) && (fabs(kappa) < 0.02))
	{
		// analyse the two shorter axis
		return correctRotation(atan2(v1y, v1x), fragPts, uCoordFrag, radius, pcaMean, true);
	}
	else
	{
		// analyse all four axis
		return correctRotation(atan2(v1y, v1x), fragPts, uCoordFrag, radius, pcaMean, false);
	}
}

std::vector<cv::Point> ShapeFeatures::thinEdgeFragment(std::vector<cv::Point> &fragPts)
{
	int noPts = fragPts.size();
	
	std::vector<int> mappingVector;
	mappingVector.resize(noPts, 1); // fill with ones

	for (int i = 0; i < noPts; i++)
	{
		// only if more than 2 following pixels
		if ((noPts-(i+1)) >= 2)
		{
			// check maximal absolute distances to next two points
			int dx1 = abs(fragPts[i].x-fragPts[i+1].x);
			int dy1 = abs(fragPts[i].y-fragPts[i+1].y);
			
			int dx2 = abs(fragPts[i].x-fragPts[i+2].x);
			int dy2 = abs(fragPts[i].y-fragPts[i+2].y);
			
			if (dx1 < 2 && dy1 < 2 && dx2 < 2 && dy2 < 2)
			{
				// one of the points can be deleted if not already done
				if ((mappingVector[i] == 1) && (mappingVector[i+1] == 1) && (mappingVector[i+2] == 1))
				{
					mappingVector[i+1] = 0;
				}
			}
		}
	}
	
	// now create new vector with thinned edge (bis this deleting single points is avoided)
	std::vector<cv::Point> thinnedEdgeFragment;
	for (int i = 0; i < noPts; i++)
	{
		if (mappingVector[i] == 1)
		{
			thinnedEdgeFragment.push_back(fragPts[i]);
		}
	}
	
	// plot values before and after edge thinning
	
	std::cout << "x_frag = [";
	for (int i = 0; i < thinnedEdgeFragment.size(); i++)
	{
		std::cout << thinnedEdgeFragment[i].x << "; ";
	}
	std::cout << "]" << std::endl;

	std::cout << "y_frag = [";
	for (int i = 0; i < thinnedEdgeFragment.size(); i++)
	{
		std::cout << thinnedEdgeFragment[i].y << "; ";
	}
	std::cout << "]" << std::endl;
	
	return thinnedEdgeFragment;
}

ShapeFeatures::~ShapeFeatures()
{
	// destructor
}




