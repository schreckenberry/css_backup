#ifndef SHAPE_FEATURES_H
#define SHAPE_FEATURES_H

#include <iostream>

#include "../util/Edge.h"

class ShapeFeatures
{
private:
	
public:

	ShapeFeatures();		//!< constructor
		
	/*
	 * Calculate characteristic orientation for an edge fragment
	 * @param fragPts: points of the fragment
	 * @param uCoordFrag: Center point of curvature extremum
	 * @param radius: radius of characteristic scale circle
	 * @param kappa: curvature at curvature extremum on lowest considered scale
	 */
	float calcOrientation(std::vector<cv::Point> &fragPts, int uCoordFrag, float radius, float kappa);
	
	/*
	 * Calculate corrected rotation so that it is inside fragment
	 * @param rot: inout orientation
	 * @param fragPts: points of the fragment
	 * @param uCoordFrag: Center point of curvature extremum
	 * @param radius: radius of characteristic scale circle
	 * @param areaMean: mean position of area without fragment points
	 * @param flat: flag signalizes very flat segment 
	 */
	float correctRotation(float rot, std::vector<cv::Point> &fragPts, int uCoordFrag, float radius, cv::Point2f areaMean, bool flat = false);
	
	/*
	 * Calculate mean of the filled area without fragment points
	 * @param fragPts: points of the fragment
	 * @param nuAreaPts: number of area without fragment points
	 * @param mean: mean position calculated by PCA
	 */
	cv::Point2f calcAreaMean(std::vector<cv::Point> &fragPts, int nuAreaPts, cv::Point2f mean);
	
	/*
	 * Calculate thinned edge fragment for given fragment
	 * @param fragPoints: points of the fragment
	 * @return: thinned edge fragment
	 */
	std::vector<cv::Point> thinEdgeFragment(std::vector<cv::Point> &fragPts);
	
	/*
	 * Crop and fill the area of a fragment for feature calculation
	 * @param fragArea: output image storing the filled fragment
	 * @param fragPoints: points of the fragment
	 */
	void fillFragmentArea(cv::Mat &fragArea, std::vector<cv::Point> &fragPts);
	
	/*
	 * Write area from image to vector structure
	 * @param fragArea: area image of a fragment
	 * @param fragAreaPts: area points of the fragment in vector structure
	 */
	void fragAreaToVector(cv::Mat &fragArea, std::vector<cv::Point> &fragAreaPts);

	~ShapeFeatures();		//!< destructor
};

#endif // SHAPE_FEATURES_H
