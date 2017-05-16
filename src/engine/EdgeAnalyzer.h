#ifndef EDGE_ANALYZER_H
#define EDGE_ANALYZER_H

#include <vector>
// Util
#include "../util/Edge.h"
#include "../util/InterestPoint.h"
// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#define SQRT_2_PI 2.50662827463100050241



class EdgeAnalyzer
{
public:
	void calculateScaleSpaces(std::vector<Edge*>* v,
			std::vector<InterestPoint*>* maximaVector,
			std::vector<InterestPoint*>* minimaVector,
			std::vector<InterestPoint*>* zeroVector);

private:
	void gauss(int x, float sigma, float * output, int normalized);
	void gauss1st(int x, float sigma, float * output, int normalized);
	void gauss2nd(int x, float sigma, float * output, int normalized);
	void calculateInterestPoints(Edge* e,
			std::vector<InterestPoint*>* maximaVector,
			std::vector<InterestPoint*>* minimaVector,
			std::vector<InterestPoint*>* zeroVector);
	Edge* resizeEdge(Edge* e, int maxEdgePoints);
	void convolution(int * input, float * output, int iLength, float * kernel, int kLength);
	inline void calcKappa(float * xCoorU, float * yCoorU, float * xCoorUU, float * yCoorUU, int length, float * output);
	// Maxima
	std::vector<int> findLocalMaxima(float * input, int length, float threshold);
	void getMaximaMax(std::vector<int>* maximaPts,
			std::vector<float>* maximaSigma,
			std::vector< std::vector<int> >* maxima,
			float threshold);
	void findOriginMax(std::vector<int>* maximaPts,
			std::vector<float>* maximaSigma,
			std::vector<int> low,
			std::vector<int> start, int index);
	// Minima
	std::vector<int> findLocalMinima(float * input, int length, float threshold);
	void getMinimaMax(std::vector<int>* minimaPts,
			std::vector<float>* minimaSigma,
			std::vector< std::vector<int> >* minima,
			float threshold);
	void findOriginMin(std::vector<int>* minimaPts,
			std::vector<float>* minimaSigma,
			std::vector<int> low,
			std::vector<int> start, int index);
	// Zero Crossings
	std::vector<int> findZeroCrossings(float * input, int length);
	void getZeroMax(std::vector<int>* zeroPts,
			std::vector<float>* zeroSigma,
			std::vector< std::vector<int> >* zero,
			float threshold);
	void findSimilar(std::vector<int>* zeroPts,
			std::vector<float>* zeroSigma,
			std::vector<int> indices, int diff, int index);


	float max(float * input, int length);
	float min(float * input, int length);
	int checkIfExists(std::vector<int> arr, int x);
	int findClosest(std::vector<int> low, int x);
	int indexOfSmallestElement(std::vector<int> array);

};

#endif // EDGE_ANALYZER_H
