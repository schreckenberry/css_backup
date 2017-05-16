#ifndef APPRCSS_H
#define APPRCSS_H

#include <vector>
// Util
#include "../util/Edge.h"
#include "../util/InterestPoint.h"
#include "../util/CssPoint.h"
// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

class CSSApproximation
{
private:
	std::vector<int> uMinima;		//!< vector to store u-coordinates of curvature minima for current sigma
	std::vector<int> uMaxima;		//!< vector to store u-coordinates of curvature maxima for current sigma
	float thresLocalMinMax;			//!< threshold which was calculated for selection of local minima and maxima
	
	/*
	 *  vectors with signatures of local minima and maxima in the CSS
	 */
	std::vector< std::vector<CssPoint> > tracedMaxSignatures;
	std::vector< std::vector<CssPoint> > tracedMinSignatures;
	
	/*
	 *  corresponding characteristic rotations (one for every signature in corresponding order)
	 *  in radians
	 */
	std::vector<float> rotMinSignatures;
	std::vector<float> rotMaxSignatures;
	
	std::vector<float> paddedIntConX;	//Save the padded integralcontour of X for re-use
	std::vector<float> paddedIntConY;	//Save the padded integralcontour of X for re-use

public:
	 /*
	  * Constructor
	  */
	CSSApproximation();
	
	void calculateCSSAppr(std::vector<Edge*>* v);

	void calculateCSSAInterestPoints(Edge* e);

	inline void calcKappa(std::vector<float>& xCoorU, std::vector<float>& yCoorU, std::vector<float>& xCoorUU,
					std::vector<float>& yCoorUU, int length, std::vector<float>& output);
	
	/**
	 * Calculate integralcontour
	 * @param sig
	 * @param res
	 */
	void calculateIntegralContour(std::vector<float>& sig, std::vector<float>& res);
	
	/**
	 * Calculate positive_modulo
	 * @param i 
	 * @param n
	 */
	inline int positive_modulo(int i, int n);
	
	/**
	 * Approximation of convolution of a signal with a gaussian function (FAGS Algorithm)
	 * @param sig input signal for convolution
	 * @param sigma standard deviation of the gaussian
	 * @param nFlt 
	 * @param norm normalization true/false
	 * @param res result of the approximation
	 */
	void approximateConvGauss(std::vector<int>& sig, float sigma, int nFlt, bool norm, std::vector<float>& res);
	
	/**
	 * Approximation of convolution of a signal with the first derivative of a gaussian function
	 * Loosly based on SURF -> Approximation with two rects
	 * @param sig input signal for convolution
	 * @param sigma standard deviation of the gaussian
	 * @param res result of the approximation
	 */
	void approximateConvGauss1st(std::vector<int>& sig, float sigma, std::vector<float>& res);
	
	/**
	 * Approximation of convolution of a signal with the second derivative of a gaussian function
	 * Loosly based on SURF -> Approximation with three rects
	 * @param sig: input signal for convolution
	 * @param sigma: standard deviation of the gaussian
	 * @param res: result of the approximation
	 */
	void approximateConvGauss2nd(std::vector<int>& sig, float sigma, std::vector<float>& res);
	
	/**
	 * Approximation of convolution of a signal with a gaussian function (FAGS Algorithm) and its derivatives
	 * Idea: (f*g)'=f'*g=f*g', where ' denotes the derivative and * the convolution
	 * @param sig input signal for convolution
	 * @param sigma standard deviation of the gaussian
	 * @param nFlt 
	 * @param norm normalization true/false
	 * @param res result of the approximation of the convolution f*g
	 * @param res1st result of the approximation of the convolution f'*g
	 * @param res2nd result of the approximation of the convolution f''*g
	 */
	
	void approximateConvGaussAll(std::vector<int>& sig, float sigma, int nFlt, bool norm,
					       std::vector<float>& res, std::vector<float>& res1st, std::vector<float>& res2nd);
	
	
	/**
	 * Subfunction for the approximations
	 * Filters a signal with a rectangle/box nFlt-times -> Approximation of the gaussian filtering
	 * Makes use of integral contours
	 * @param sig input signal
	 * @param nFlt
	 * @param sigma
	 * @param res
	 */
	void rectFiltering(std::vector<int>& sig, float sigma, int nFlt, std::vector<float>& res);
	
	/**
	 * Subfunction for the approximation (2nd derivative)
	 * @param sig input signal
	 * @param winSize size of window for filtering
	 * @param res result
	 */
	void movingAveraging(std::vector<float>& sig, int winSize, std::vector<float>& res);
	
	/**
	 * Subfunction for padding the signal
	 * @param sig input signal
	 * @param paddedSig padded signal
	 * @param padSize number of bins which are added before and behind the signal
	 */
	void padSignal(std::vector<int>& sig, std::vector<float>& paddedSig, int padSize);
	bool intConXCalc;
	bool intConYCalc;
	
	/**
	 * And the show must go on...
	 */
	void findLocalMinMax(std::vector<float>& kappa, int length, float sigma, bool isClosed);
	void traceMinMaxSignatures(std::vector<float>& kappa, int length, float sigma, bool isClosed);
	void traceMinMaxSignaturesV2(std::vector<float>& kappa, int length, float sigma, bool isClosed);
	void calculateShapeFeatures(Edge* e);
	void deleteShortSignatures(int minLength);
	
	inline float calcThreshold(float sigma);
	
	/**
	 * method destructor
	 */
	~CSSApproximation();
};

#endif // APPRCSS_H
