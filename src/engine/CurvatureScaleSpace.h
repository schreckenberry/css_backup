#ifndef CURVATURE_SCALE_SPACE_H
#define CURVATURE_SCALE_SPACE_H

#include <vector>
// Util
#include "../util/Edge.h"
#include "../util/InterestPoint.h"
#include "../util/CssPoint.h"
// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#define SQRT_2_PI 2.50662827463100050241

class CurvatureScaleSpace
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

public:
	/*
	 * constructor
	 */
	CurvatureScaleSpace();
	
	void calculateScaleSpaces(std::vector<Edge*>* v);

	void calculateCSSInterestPoints(Edge* e);

	//TODO: test and make output to vector
    /** Calculates values of the gaussian first derivatve for a given standard deviation
     *  @param xCoorU       first derivative of x coordinates
     *  @param yCoorU       first derivative of y coordinates
     *  @param xCoorUU      second derivative of x coordinates
     *  @param yCoorUU      second derivative of y coordinates
     *  @param length       length of input arrays //TODO move into function
     *  @param output       array which stores the result
     */
     // inline, d. h. compiler kopiert funktion direkt an jeweilige stelle (statt sprung)
	inline void calcKappa(std::vector<float>* xCoorU, std::vector<float>* yCoorU,
			std::vector<float>* xCoorUU, std::vector<float>* yCoorUU, int length, float* output);

	/*
	 * Calculates values of the gaussian for a given standard deviation
	 * with appropriate length (e. g. 4.5 * sigma)
	 * @param sigma: value for the standard deviation
	 * @param output: vector with result
	 * @param normalized: 0 or 1, if result should be normalized
	 */
	void gauss(float sigma, std::vector<float>* output, bool normalized);

	/*
	 * Calculates values of the MIRRORED 1st derivate of gaussian for a given
	 * standard deviation with appropriate length (e. g. 4.5 * sigma)
	 * @param sigma: value for the standard deviation
	 * @param output: vector with result
	 * @param normalized: 0 or 1, if result should be normalized
	 */
	void gauss1st(float sigma, std::vector<float>* output, bool normalized);

	/*
	 * Calculates values of the MIRRORED 2nd derivate of gaussian for a given
	 * standard deviation with appropriate length (e.g. 4.5 * sigma)
	 * @param sigma: value for the standard deviation
	 * @param output: vector with result
	 * @param normalized: 0 or 1, if result should be normalized
	 */
	void gauss2nd(float sigma, std::vector<float>* output, bool normalized);
	
	/* 
	 * Convolution with constant padding.
	 * NOTE: actually it is a correlation so feed rotated signals (gauss1st and gauss2nd are mirrord)
	 * @param sig: input signal for convolution
	 * @param kern: kernel for convolution (must be odd)
	 * @param res: result of convolution with size of sig
	 */
	void convolve(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res);

	/*
	 * Convolution without padding.
	 * NOTE: actually it is a correlation so feed 180 deg rotated signals
	 * (derivated gaussians are rotated here)
	 * @param sig: input signal for convolution
	 * @param kern: kernel for convolution (must be odd)
	 * @param res: result of convolution with size of sig
	 */
	void convolveClosedContour1(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res);
	void convolveClosedContour2(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res);
	void convolveClosedContour3(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res);
	void convolveClosedContour4(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res);
	
	/* Print estimated status of convolution.
	 */
	void printConvolutionStatus(float currentSigma);

	/* Determines local minima and maxima for the given input array.
	 * Coordinates are stored in member variables uMinima and uMaxima.
	 * @param kappa: pointer to array with current values of kappa
	 * @param length: length of array kappa
	 * @param sigma: current sigma for which kappa was calculated
	 * @param isClosed: defines if contour is closed
	 */
	void findLocalMinMax(float* kappa, int length, float sigma, bool isClosed);
	
	/* Calculate threshold based on current sigma.
	 * @param sigma: current sigma
	 */
	float calcThreshold(float sigma);
	
	/*
	 * Trace the signatures of local minima and maxima in the CSS
	 * @param kappa: pointer to array with current values of kappa
	 * @param length: length of array kappa
	 * @param sigma: current sigma for which kappa was calculated
	 * @param isClosed: defines if contour is closed
	 */
	void traceMinMaxSignatures(float* kappa, int length, float sigma, bool isClosed);
	
	/*
	 * Delete all signatures which store less than minLength datapoints
	 * @param minLength: Minimal length
	 */
	void deleteShortSignatures(int minLength = 2);
	
	/*
	 * Calculate Features for the Contour Fragments
	 * @param e: current edge
	 */
	void calculateShapeFeatures(Edge* e);
	
	/*
	 * method destructor
	 */
	~CurvatureScaleSpace();
};

#endif // CURVATURE_SCALE_SPACE_H
