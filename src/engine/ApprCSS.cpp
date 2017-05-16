#include "ApprCSS.h"
#include <math.h>
#include <iostream>
//for GUI stuff can be deleted if no imshow() is used
//#include <opencv2/highgui/highgui.hpp>
//using namespace cv;

#include "./../util/Visualization.h"
Visualization visual_t;

#include "./../util/ShapeFeatures.h"
ShapeFeatures shapeFeat_t;

#include "./../engine/CurvatureScaleSpace.h"
CurvatureScaleSpace css_t;

float lF = 4.5;					// lengthFactor * sigma samples at each side of max of gaussian
float sigmaM = 200.0; 				// max. sigma
float sigmaStepSize = 2.0; 			// sigma step size, default 2.0
float sigmaSt = 8.0; 				// min. sigma, default 8.0
float maxSigmaDis = 3.0 * sigmaStepSize; 	// tracing of extrema signatures (see traceMinMaxSignatures)
int minimumSigmaSteps = 8; 			// delete signatures shorter than N sigma steps, default 8
const float sF = 0.3; 				// scale factor for radius


float maxPadSize = std::ceil(5.5*sigmaM);
static int *adressX = NULL;
static int *adressY = NULL;


CSSApproximation::CSSApproximation()
{
	//std::cout << "CSSApproximation object created" << std::endl;
}

void CSSApproximation::calculateCSSAppr(std::vector< Edge* >* v)
{
	// delete old tracing results
	tracedMinSignatures.clear();
	tracedMaxSignatures.clear();
	
	rotMinSignatures.clear();
	rotMaxSignatures.clear();
	
	for (std::vector<Edge*>::iterator it = v->begin(); it != v->end(); ++it)
	{
		calculateCSSAInterestPoints((*it));
	}
}

void CSSApproximation::calculateCSSAInterestPoints(Edge* e)
{
	int eLength = int(e->length()); // length of edge

	// Get X and Y values seperately
	std::vector<int> X;
	std::vector<int> Y;
	X = e->getEdgeXCoordinatesVector();
	Y = e->getEdgeYCoordinatesVector();
	
	std::vector<float> _xAppr;
	std::vector<float> _xApprU;
	std::vector<float> _xApprUU;
	std::vector<float> _yAppr;
	std::vector<float> _yApprU;
	std::vector<float> _yApprUU;
	std::vector<float> _kappaA;
	_kappaA.resize(eLength);

	bool normState = true;
	/************************************************/
// 	std::vector<float> gauss;
// 	float sigma2=sigmaSt;
// 	while(sigma2<=sigmaM)
// 	{
// 		css_t.gauss1st(sigma2, &gauss, normState);
// 		css_t.convolveClosedContour4(&X, &gauss, &_xApprU);
// 		css_t.convolveClosedContour4(&Y, &gauss, &_yApprU);
// 		css_t.gauss2nd(sigma2, &gauss, normState);
// 		css_t.convolveClosedContour4(&X, &gauss, &_xApprUU);
// 		css_t.convolveClosedContour4(&Y, &gauss, &_yApprUU);
// 		calcKappa(_xApprU, _yApprU, _xApprUU, _yApprUU, eLength, _kappaA);
// 		
// 		visual_t.saveCSVConvolution(_yApprU, 1);
// 		visual_t.saveCSVConvolution(_yApprUU, 2);
// 		visual_t.saveCSVKappa(_kappaA, eLength, sigma2);
// 		sigma2+=sigmaStepSize;
// 	}
	/************************************************/
	
	intConXCalc = false;
	intConYCalc = false;
	
	float testScale = 8;
	
	float sigma=sigmaSt;
	while(sigma<=sigmaM)
	{
// 		std::cout << "sigma=" << sigma << std::endl;
		approximateConvGauss1st(X, sigma, _xApprU);
		approximateConvGauss2nd(X, sigma, _xApprUU);
		approximateConvGauss1st(Y, sigma, _yApprU);
		approximateConvGauss2nd(Y, sigma, _yApprUU);
		
		calcKappa(_xApprU, _yApprU, _xApprUU, _yApprUU, eLength, _kappaA);
		
// 		visual_t.saveCSVApproximation(_yApprU, 1);
// 		visual_t.saveCSVApproximation(_yApprUU, 2);
// 		visual_t.saveCSVKappaAppr(_kappaA, eLength, sigma);

		findLocalMinMax(_kappaA, eLength, sigma,e->getIsLoop());
		traceMinMaxSignaturesV2(_kappaA, eLength, sigma,e->getIsLoop());
		
		visual_t.saveCSVMinMax(&uMinima, &uMaxima, sigma);
		
		/***********************************************************/
		//NOTE: Development purposes
// 		if (sigma == testScale)
// 		{
// 			std::vector<float> _xSmoothed;
// 			std::vector<float> _ySmoothed;
// 			
// 			approximateConvGauss(X, sigma, 5, true, _xSmoothed);
// 			approximateConvGauss(Y, sigma, 5, true, _ySmoothed);
// 			
// 			Convert to Edge object for drawing gauss filtered contours
// 			std::vector<cv::Point> edgePoints;
// 			cv::Point p;
// 			for (int i = 0; i < eLength; i++)
// 			{
// 				p.x = std::round(_xSmoothed[i]);
// 				p.y = std::round(_ySmoothed[i]);
// 				std::cout << "p = (" << _xSmoothed[i] << ", " << _ySmoothed[i] << ")" << std::endl;
// 				
// 				edgePoints.push_back(p);
// 			}
// 			Edge edgeGauss(edgePoints);
// 			
// 			write SVG file of colored smoothed curve
// 			visual_t.saveSVGcolorWheel(&edgeGauss, sigma);
// 
// 			write Matlab file		
// 			visual_t.saveMatlabKappaVec(_kappaA, eLength, sigma, uMinima, uMaxima, thresLocalMinMax);
// 			std::cout << "sigmaPlotMatlab = " << sigma << std::endl;
// 			std::cout << "thresLocalMinMax = " << thresLocalMinMax << std::endl;
// 		}
		/***********************************************************/
		
		sigma+=sigmaStepSize;
	}
	
	std::cout << "tracedMinSignatures.size() = " << tracedMinSignatures.size() << std::endl;
	std::cout << "tracedMaxSignatures.size() = " << tracedMaxSignatures.size() << std::endl;
// 	
	// delete signatures shorter than N sigma sigmaStepSize
	deleteShortSignatures(minimumSigmaSteps);
	
	// Tracing finished here - tracedMin/MaxSignatures are cleared for every edge
	calculateShapeFeatures(e);
	
	visual_t.generateColorWheel(eLength);
	visual_t.saveCSVTracedMinMax(tracedMinSignatures, tracedMaxSignatures);
	visual_t.saveSvgCssResult(e, tracedMinSignatures, tracedMaxSignatures, rotMinSignatures, rotMaxSignatures, sF);
}

void CSSApproximation::calcKappa(std::vector<float>& xCoorU, std::vector<float>& yCoorU, std::vector<float>& xCoorUU,
					std::vector<float>& yCoorUU, int length, std::vector<float>& output)
{
	output.resize(length);
	for (int i=0; i<length; i++)
		output[i] = (xCoorU[i]*yCoorUU[i] - xCoorUU[i]*yCoorU[i]) / 
					std::pow( std::pow(xCoorU[i], 2.0) + std::pow(yCoorU[i], 2.0), 1.5);
}

void CSSApproximation::rectFiltering(std::vector<int>& sig, float sigma, int nFlt, std::vector<float>& res)
{
	/*
	 * Calculating the paramters
	 * wIdeal: ideal width of the rect
	 * wl: first odd valued integer less than wIdeal
	 * wu: next odd value > wl
	 * For derivation of m refer to paper "Fast Almost Gaussian Smoothing"
	 * Calculation of actual sigma optional
	 */
	float wIdeal = std::sqrt(12*sigma*sigma/nFlt+1);
	int wl = std::floor(wIdeal);
	if(!(wl%2))
		wl=wl-1;
	int wu = wl+2;
	int mFlt=round((12*sigma*sigma - nFlt*wl*wl - 4*nFlt*wl - 3*nFlt)/(-4*wl - 4));
	//float sigmAct = std::sqrt((mFlt*wl*wl + (nFlt-mFlt)*wu*wu - nFlt)/12);
	//std::cout << "sigma=" << sigma << "; wl=" << wl << "; m=" << mFlt << "; sigmaAct=" << sigmAct << std::endl;
	
	/*
	 * Some variables and vectors
	 */
	int sigSize = sig.size();
	std::vector<float> intConTmp;
	/*
	 * Resizing of the result vector res
	 */
	res.resize(sigSize);	
	
	//TODO: Better way for init?		
	for(int i=0; i<sigSize; i++)
		res[i]=sig[i];
	
	/*
	 * Filtering with width=wl mFlt-times
	 * up/down are the lower/upper borders of the rect
	 */
	int up = (wl-1)/2;
	int down = up+1;
	if(up!=0)
	{
		for(int j=0; j<mFlt; j++)
		{
			calculateIntegralContour(res, intConTmp);
			for(int i=down; i<(sigSize-up); i++)
				res[i]=(intConTmp[i+up]-intConTmp[i-down])/wl;
		}
	}
	/*
	 * Filtering with width=wu (nFlt-mFlt)-times
	 * up/down are the lower/upper borders of the rect
	 */
	up=(wu-1)/2;
	down=up+1;
	if(up!=0)
	{
		for(int j=0; j<(nFlt-mFlt); j++)
		{
			calculateIntegralContour(res, intConTmp);
			for(int i=down; i<(sigSize-up); i++)
				res[i]=(intConTmp[i+up]-intConTmp[i-down])/wu;
		}
	}
}

void CSSApproximation::approximateConvGauss(std::vector<int>& sig, float sigma, int nFlt, bool norm, std::vector<float>& res)
{
	// TODO: Smaller Padding Size?! See Paper Fast Almost Gaussian Smoothing(?)
	int sigSize=sig.size();
	int maxPadSize = 1.5*sigmaM;
	//float wIdeal = std::sqrt(12*sigma*sigma/nFlt+1);
	//int maxPadSize = std::round(nFlt*(wIdeal-1)/2);
	
	std::vector<int> paddedSig;
	std::vector<float> paddedRes;
	res.resize(sigSize);
	
	paddedSig.clear();
	int index = 0;
	
	// pad beginning of signal
	for (int i = 0; i < maxPadSize; i++)
	{
		index = positive_modulo(sigSize-maxPadSize+i, sigSize);
		paddedSig.push_back(sig[index]);
		//std::cout << "index beginning = " << index << std::endl;
	}
	// write signal
	for (int i = 0; i < sigSize; i++)
	{
		paddedSig.push_back(sig[i]);
	}
	// pad end of signal
	for (int i = 0; i < maxPadSize; i++)
	{
		index = positive_modulo(sigSize+i, sigSize);
		paddedSig.push_back(sig[index]);
		//std::cout << "index end = " << index << std::endl;
	}
	
	rectFiltering(paddedSig, sigma, nFlt, paddedRes);
	
	for(int i=0; i<sigSize; i++)
		res[i] = paddedRes[i+maxPadSize];
}

void CSSApproximation::approximateConvGauss1st(std::vector<int>& sig, float sigma, std::vector<float>& res)
{
	//NOTE: Implementation like SURF
	int w = std::round(5*sigma);
	if (w%2 == 0)
		w--;
	int whalf=(w-1)/2;
	
	// Padding (executed only once...)
	int sigSize = sig.size();
	int offset = maxPadSize-w;
	
	static std::vector<float> paddedSig;
	std::vector<float> paddedRes;
	static std::vector<float> intCon;
	
	const int longSigSize = sigSize+2*maxPadSize;
	paddedRes.resize(longSigSize);
	res.resize(sigSize);
	
	if(!intConXCalc)
	{
		adressX = &sig[0];
		paddedSig.clear();
		padSignal(sig, paddedSig, maxPadSize);
		calculateIntegralContour(paddedSig, paddedIntConX);
		intConXCalc = true;
		
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConX[i]-paddedIntConX[i-whalf])/whalf
						- (paddedIntConX[i+whalf]-paddedIntConX[i])/whalf)/(-1.92);
		}
	}
	else if ((!intConYCalc) && (&sig[0] != adressX))
	{
		adressY = &sig[0];
		paddedSig.clear();
		padSignal(sig, paddedSig, maxPadSize);
		calculateIntegralContour(paddedSig, paddedIntConY);
		intConYCalc = true;
		
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConY[i]-paddedIntConY[i-whalf])/whalf
						- (paddedIntConY[i+whalf]-paddedIntConY[i])/whalf)/(-1.92);
		}
	}
	else if (intConXCalc && (&sig[0] == adressX))
	{
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConX[i]-paddedIntConX[i-whalf])/whalf
						- (paddedIntConX[i+whalf]-paddedIntConX[i])/whalf)/(-1.92);
		}
	}
	else if (intConYCalc && (&sig[0] == adressY))
	{
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConY[i]-paddedIntConY[i-whalf])/whalf
						- (paddedIntConY[i+whalf]-paddedIntConY[i])/whalf)/(-1.92);
		}
	}
	else
	{
		std::cout << "Whoooops. Some padding went wrong!" << std::endl;
			// and do the standard padding		
		padSignal(sig, paddedSig, w);
		int longSigSizeStd = paddedSig.size();
		
		// Filtering
		std::vector<float> intConTmp;
		calculateIntegralContour(paddedSig, intConTmp);
		
		paddedRes.resize(longSigSizeStd);
		for(int i=whalf; i<longSigSizeStd-whalf; i++)
		{
			paddedRes[i]=((intConTmp[i]-intConTmp[i-whalf])/whalf-(intConTmp[i+whalf]-intConTmp[i])/whalf)/(-1.92);
		}
	}
	
	for(int i=0; i<sigSize; i++)
			res[i] = paddedRes[i+offset+w];
}

void CSSApproximation::approximateConvGauss2nd(std::vector<int>& sig, float sigma, std::vector<float>& res)
{
	//NOTE: Implementation like SURF
	int w = std::round(sigma*5.5);
	if (w%2 == 0)
		w--;
	int whalf=(w-1)/2;
	
	int w3 = std::round(w/3);
	if (w3%2 == 0)
		w3--;
	
	int w3half = (w3-1)/2;
	
		
	// Padding (executed only once...)
	int sigSize = sig.size();
	int offset = maxPadSize-w;
	
	static std::vector<float> paddedSig;
	std::vector<float> paddedRes;
	static std::vector<float> intCon;
	
	const int longSigSize = sigSize+2*maxPadSize;
	paddedRes.resize(longSigSize);
	res.resize(sigSize);
		
	
	if(!intConXCalc)
	{
		adressX = &sig[0];
		paddedSig.clear();
		padSignal(sig, paddedSig, maxPadSize);
		calculateIntegralContour(paddedSig, paddedIntConX);
		intConXCalc = true;
		
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConX[i+whalf]-paddedIntConX[i+w3half])/(whalf-w3half)
					- 2*(paddedIntConX[i+w3half]-paddedIntConX[i-w3half])/(w3half+w3half)
						+ (paddedIntConX[i-w3half]-paddedIntConX[i-whalf])/(whalf-w3half))/3.6;
		}
	}
	else if ((!intConYCalc) && (&sig[0] != adressX))
	{
		adressY = &sig[0];
		paddedSig.clear();
		padSignal(sig, paddedSig, maxPadSize);
		calculateIntegralContour(paddedSig, paddedIntConY);
		intConYCalc = true;
		
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConY[i+whalf]-paddedIntConY[i+w3half])/(whalf-w3half)
					- 2*(paddedIntConY[i+w3half]-paddedIntConY[i-w3half])/(w3half+w3half)
						+ (paddedIntConY[i-w3half]-paddedIntConY[i-whalf])/(whalf-w3half))/3.6;
		}
	}
	else if (intConXCalc && (&sig[0] == adressX))
	{
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConX[i+whalf]-paddedIntConX[i+w3half])/(whalf-w3half)
					- 2*(paddedIntConX[i+w3half]-paddedIntConX[i-w3half])/(w3half+w3half)
						+ (paddedIntConX[i-w3half]-paddedIntConX[i-whalf])/(whalf-w3half))/3.6;
		}
	}
	else if (intConYCalc && (&sig[0] == adressY))
	{
		for(int i=offset+whalf; i<longSigSize-offset-whalf; i++)
		{
			paddedRes[i] = ((paddedIntConY[i+whalf]-paddedIntConY[i+w3half])/(whalf-w3half)
					- 2*(paddedIntConY[i+w3half]-paddedIntConY[i-w3half])/(w3half+w3half)
						+ (paddedIntConY[i-w3half]-paddedIntConY[i-whalf])/(whalf-w3half))/3.6;
		}
	}
	else
	{
		std::cout << "Whoooops. Some padding went wrong!" << std::endl;
			// and do the standard padding		
		padSignal(sig, paddedSig, w);
		int longSigSizeStd = paddedSig.size();
		
		// Filtering
		std::vector<float> intConTmp;
		calculateIntegralContour(paddedSig, intConTmp);
		
		paddedRes.resize(longSigSizeStd);
		for(int i=whalf; i<longSigSizeStd-whalf; i++)
		{
			paddedRes[i] = ((intConTmp[i+whalf]-intConTmp[i+w3half])/(whalf-w3half)
					- 2*(intConTmp[i+w3half]-intConTmp[i-w3half])/(w3half+w3half)
						+ (intConTmp[i-w3half]-intConTmp[i-whalf])/(whalf-w3half))/3.6;
		}
	}
	
	// Smoothing for smaller sigmas
	// Otherwise findLocalMinMax finds a lot of "min's" or "max's" which aren't really local minima or maxima
	// Happens probably due to approximation errors
	if(sigma<20)
	{
		std::vector<float> pad2ndTmp;
		pad2ndTmp.resize(longSigSize);
		
		for(int i=0; i<longSigSize; i++)
			pad2ndTmp[i] = paddedRes[i];
		movingAveraging(pad2ndTmp, 11, paddedRes);
	}
	
	for(int i=0; i<sigSize; i++)
			res[i] = paddedRes[i+offset+w];
}

void CSSApproximation::approximateConvGaussAll(std::vector<int>& sig, float sigma, int nFlt, bool norm,
							std::vector<float>& res, std::vector<float>& res1st, std::vector<float>& res2nd)
{
	int sigSize=sig.size();
	float wIdeal = std::sqrt(12*sigma*sigma/nFlt+1);
	int maxPadSize = std::round(nFlt*(wIdeal-1)/2);
	
	std::vector<int> paddedSig;
	std::vector<float> paddedRes;
	res.resize(sigSize);
	res1st.resize(sigSize);
	res2nd.resize(sigSize);
	
	paddedSig.clear();
	int index = 0;
	
	// pad beginning of signal
	for (int i = 0; i < maxPadSize; i++)
	{
		index = positive_modulo(sigSize-maxPadSize+i, sigSize);
		paddedSig.push_back(sig[index]);
		//std::cout << "index beginning = " << index << std::endl;
	}
	// write signal
	for (int i = 0; i < sigSize; i++)
	{
		paddedSig.push_back(sig[i]);
	}
	// pad end of signal
	for (int i = 0; i < maxPadSize; i++)
	{
		index = positive_modulo(sigSize+i, sigSize);
		paddedSig.push_back(sig[index]);
		//std::cout << "index end = " << index << std::endl;
	}
	
	rectFiltering(paddedSig, sigma, nFlt, paddedRes);
	
	int longSigSize = paddedSig.size();
	std::vector<float> padded1st;
	std::vector<float> padded2nd;
	padded1st.resize(longSigSize);
	padded2nd.resize(longSigSize);
	
	//Calculate first derivative
	for (int i=0; i<longSigSize-1; i++)
		padded1st[i]=paddedRes[i+1]-paddedRes[i];
	padded1st[longSigSize-1] = padded1st[longSigSize-2];
	
	//Calculate second derivative (w/0 smoothing)
	for (int i=0; i<longSigSize-1; i++)
		padded2nd[i]=padded1st[i+1]-padded1st[i];
	padded2nd[longSigSize-1] = padded2nd[longSigSize-2];
	
	//Calculate second derivative (w/ smoothing)
	std::vector<float>pad2ndTmp;
	pad2ndTmp.resize(longSigSize);
	for(int i=0; i<longSigSize; i++)
		pad2ndTmp[i] = padded2nd[i];
	movingAveraging(pad2ndTmp, 15, padded2nd);
	
	if(!norm) // No Normalization
	{
		for(int i=0; i<sigSize; i++)
		{
			res[i] = paddedRes[i+maxPadSize];
			res1st[i] = padded1st[i+maxPadSize];
			res2nd[i] = padded2nd[i+maxPadSize];
		}
	}
	
	if(norm) // Normalization
	{
		for(int i=0; i<sigSize; i++)
		{
			res[i] = paddedRes[i+maxPadSize];
			res1st[i] = 1.2534*sigma*padded1st[i+maxPadSize];
			res2nd[i] = 1.05*sigma*sigma*padded2nd[i+maxPadSize];
		}
	}
}

void CSSApproximation::calculateIntegralContour(std::vector<float>& sig, std::vector<float>& res)
{
	int sigSize = sig.size();
	res.resize(sigSize);
	res[0] = sig[0];
	for(int i=1; i<sigSize; i++)
	{
		res[i] = res[i-1]+sig[i];
	}
}

int CSSApproximation::positive_modulo(int i, int n)
{
    return (i % n + n) % n;
}

void CSSApproximation::movingAveraging(std::vector<float>& sig, int winSize, std::vector<float>& res)
{
	int sigSize = sig.size();
	if((winSize%2)==0)
		winSize = winSize-1;
	
	std::vector<float> intConTmp;
	calculateIntegralContour(sig, intConTmp);
	
	res.resize(sigSize);
	
	int up = (winSize-1)/2;
	int down = up+1;
	for (int i=down; i<sigSize-up; i++)
		res[i]=(intConTmp[i+up]-intConTmp[i-down])/winSize;
}

void CSSApproximation::padSignal(std::vector<int>& sig, std::vector<float>& paddedSig, int padSize)
{
	int sigSize = sig.size();
	paddedSig.clear();
	
	int index = 0;
	
	// pad beginning of signal
	for (int i = 0; i < padSize; i++)
	{
		index = positive_modulo(sigSize-padSize+i, sigSize);
		paddedSig.push_back(sig[index]);
		//std::cout << "index beginning = " << index << std::endl;
	}
	// write signal
	for (int i = 0; i < sigSize; i++)
	{
		paddedSig.push_back(sig[i]);
	}
	// pad end of signal
	for (int i = 0; i < padSize; i++)
	{
		index = positive_modulo(sigSize+i, sigSize);
		paddedSig.push_back(sig[index]);
		//std::cout << "index end = " << index << std::endl;
	}
}

void CSSApproximation::findLocalMinMax(std::vector<float>& kappa, int length, float sigma, bool isClosed)
{
	// empty the vectors for every new calculation
	uMinima.clear();
	uMaxima.clear();
	
	// find local maxima and minima
	// algorithm is inspired by http://www.billauer.co.il/peakdet.html (last checked: 2016-03-30)
// 	const float thres = 0.002;
	float thres = calcThreshold(sigma);
	const float delta = 0.001;
	//std::cout << "delta = " << delta << std::endl;
	
	// make threshold available outside this function (e. g. for plot)
	thresLocalMinMax = thres;
	
	int maxPos, minPos;
	float min = 1000000; float max = -1000000;

	// do two iterations in case of closed contours so that extrema at
	// both back and front are considered
	int nuIterations = 1;
	if (isClosed)
		nuIterations=2;

	// first iteration is for initializing min, max, and lookformax
	// values are only stored in second run
	bool lookForMax = true;
	for (int h = 0; h < nuIterations; h++)
	{
		for (int i = 0; i < length; i++)
		{
			if(kappa[i]>max)
			{
				max = kappa[i];
				maxPos = i;
			}
			if(kappa[i]<min)
			{
				min = kappa[i];
				minPos = i;
			}
			
			if(lookForMax)
			{
				if (kappa[i]<(max-delta))
				{
					bool atBorder = !isClosed && (maxPos == 0 || maxPos == length-1);
					if((max>thres) && (h>0 || !isClosed) && !atBorder)
						uMaxima.push_back(maxPos);
					min = kappa[i];
					minPos=i;
					lookForMax = false;
				}
			}
			else
			{
				if(kappa[i]>(min+delta))
				{
					bool atBorder = !isClosed && (minPos == 0 || minPos == length-1);
					if (min < (-1.0 * thres) && (h > 0 || !isClosed) && !atBorder)
						uMinima.push_back(minPos);
					max=kappa[i];
					maxPos=i;
					lookForMax = true;
				}
			}
		}
	}
// 	std::cout << "Found " << uMaxima.size() << " Maxima" << std::endl;
// 	std::cout << "Found " << uMinima.size() << " Minima" << std::endl;
}

void CSSApproximation::traceMinMaxSignatures(std::vector<float>& kappa, int length, float sigma, bool isClosed)
{
	// this function is called for every sigma
	// because local min and max are always calculated for one specific sigma
	// and not all kappa values are stored (only corresponding u-coordinates)
	
	// data structure
	/*  [cssPoint] = [u, sigma, kappa]
	 *  t  +------------------------------------------------------+
	 *  r  | maxSignature: [cssPoint] [cssPoint] [cssPoint] ...   |
	 *  a  | maxSignature: [cssPoint] [cssPoint]                  |
	 *  c  | maxSignature: [cssPoint] [cssPoint] [cssPoint] ...   |
	 *  e  | ...                                                  |
	 *  d  +------------------------------------------------------+
	 *  MaxSignatures
	 */
	CssPoint cssPoint;

	/*********************************************************************/
	// NOTE: member variable vector uMaxima only stores values for current sigma
	for (size_t i = 0; i < uMaxima.size(); i++)
	{
		cssPoint.u = uMaxima[i];
		cssPoint.sigma = sigma;
		cssPoint.kappa = kappa[cssPoint.u];
		
		// initialize if no previous values (just push back)
		if (sigma == sigmaSt)
		{
			// new vector for every point
			std::vector<CssPoint> maxSignature;
			maxSignature.push_back(cssPoint);
			tracedMaxSignatures.push_back(maxSignature);
			
			/*
			std::cout << "cssPoint.u = " << cssPoint.u << std::endl;
			std::cout << "cssPoint.sigma = " << cssPoint.sigma << std::endl;
			std::cout << "cssPoint.kappa = " << cssPoint.kappa << std::endl;
			*/
		}
		// otherwise find correspondences
		else
		{			
			// find next one in previous vectors, take last elements
			// note: not most effective because only signatures which are large enough must be checked
			int lastMinUDistance = 1000000;
			int jMin = -1;

			// go trough all last
			int uDistance = 0;
			float sigmaDistance = 0;

			for (int j = 0; j < tracedMaxSignatures.size(); j++)
			{
				// calculate distances
				sigmaDistance = sigma - tracedMaxSignatures[j].back().sigma;
				uDistance = std::abs(cssPoint.u - tracedMaxSignatures[j].back().u);
				int maxUDistance = 0.2 * sigma + 10;
				
				// in case of closed contours use circular distance
				if ((uDistance > (length / 2)) && isClosed)
				{
					uDistance = length - uDistance;
				}
				
				// sigmaDistance is used to avoid false assignments in case of
				// overlapping curvy signatures
				if ((uDistance < lastMinUDistance) && (sigmaDistance <= maxSigmaDis) && (uDistance <= maxUDistance))
				{
					lastMinUDistance = uDistance;
					jMin = j;
				}
			}

			// store in vector to which distance was minimal
			// only assign if valid index found, otherwise throw away (or create new signature?)
			if (jMin != -1)
			{
				tracedMaxSignatures[jMin].push_back(cssPoint);

				//if (sigma == 10.0)
				//{
					//std::cout << cssPoint.u << ", ";
					//std::cout << uMaxima.size() << std::endl;
				//}
			}
			else
			{
				std::cout << "Warning: Could not trace local maximum (traceMinMaxSignatures) at sigma=" << sigma << std::endl;
				// TODO: maybe create new signature here?
			}
		}
	}
	/*********************************************************************/

	// do the same for the minima (this is just a copy of the upper part)
	/*********************************************************************/
	// Note: member variable vector uMinima only stores values for current sigma
	for (size_t i = 0; i < uMinima.size(); i++)
	{
		cssPoint.u = uMinima[i];
		cssPoint.sigma = sigma;
		cssPoint.kappa = kappa[cssPoint.u];
		
		// initialize if no previous values (just push back)
		if (sigma == sigmaSt)
		{
			// new vector for every point
			std::vector<CssPoint> minSignature;
			minSignature.push_back(cssPoint);
			tracedMinSignatures.push_back(minSignature);
			
			/*
			std::cout << "cssPoint.u = " << cssPoint.u << std::endl;
			std::cout << "cssPoint.sigma = " << cssPoint.sigma << std::endl;
			std::cout << "cssPoint.kappa = " << cssPoint.kappa << std::endl;
			*/
		}
		// otherwise find correspondences
		else
		{			
			// find next one in previous vectors, take last elements
			// note: not very effective because only signatures which are large enough must be checked
			int lastMinUDistance = 1000000;
			int jMin = -1;

			// go trough all last
			int uDistance = 0;
			float sigmaDistance = 0;

			for (int j = 0; j < tracedMinSignatures.size(); j++)
			{
				// calculate distances
				sigmaDistance = sigma - tracedMinSignatures[j].back().sigma;
				uDistance = std::abs(cssPoint.u - tracedMinSignatures[j].back().u);
				int maxUDistance = 0.2 * sigma + 10;
				
				// in case of closed contours use circular distance
				if ((uDistance > (length / 2)) && isClosed)
				{
					uDistance = length - uDistance;
				}
				
				// sigmaDistance is used to avoid false assignments in case of
				// overlapping curvy signatures
				if ((uDistance < lastMinUDistance) && (sigmaDistance <= maxSigmaDis) && (uDistance <= maxUDistance))
				{
					lastMinUDistance = uDistance;
					jMin = j;
				}
			}

			// store in vector to which distance was minimal
			// only assign if valid index found, otherwise throw away (or create new signature?)
			if (jMin != -1)
			{
				tracedMinSignatures[jMin].push_back(cssPoint);
			}
			else
			{
				std::cout << "Warning: Could not trace local minimum (traceMinMaxSignatures) at sigma=" << sigma << std::endl;
				// TODO: maybe create new signature here?
			}
		}
	}
	/*********************************************************************/
	
	/*
	if (sigma == sigmaStart)
	{
		std::cout << "tracedMaxSignatures.size() = " << tracedMaxSignatures.size() << std::endl;
		std::cout << "tracedMinSignatures.size() = " << tracedMinSignatures.size() << std::endl;
		std::cout << "sigma = " << sigma << std::endl;
	}
	*/
}

void CSSApproximation::traceMinMaxSignaturesV2(std::vector<float>& kappa, int length, float sigma, bool isClosed)
{
	// this function is called for every sigma
	// because local min and max are always calculated for one specific sigma
	// and not all kappa values are stored (only corresponding u-coordinates)
	
	// data structure
	/*  [cssPoint] = [u, sigma, kappa]
	 *  t  +------------------------------------------------------+
	 *  r  | maxSignature: [cssPoint] [cssPoint] [cssPoint] ...   |
	 *  a  | maxSignature: [cssPoint] [cssPoint]                  |
	 *  c  | maxSignature: [cssPoint] [cssPoint] [cssPoint] ...   |
	 *  e  | ...                                                  |
	 *  d  +------------------------------------------------------+
	 *  MaxSignatures
	 */
	CssPoint cssPoint;

	/*********************************************************************/
	// NOTE: member variable vector uMaxima only stores values for current sigma
	for (size_t i = 0; i < uMaxima.size(); i++)
	{
		cssPoint.u = uMaxima[i];
		cssPoint.sigma = sigma;
		cssPoint.kappa = kappa[cssPoint.u];
		
		
		if (sigma == sigmaSt)
		{
			// Initiaization when starting
			std::vector<CssPoint> maxSignature;
			maxSignature.push_back(cssPoint);
			tracedMaxSignatures.push_back(maxSignature);
		}
		
		
		// otherwise find correspondences
		else
		{			
			// find next one in previous vectors, take last elements
			// NOTE: not most effective because only signatures which are large enough must be checked
			int lastMinUDistance = 1000000;
			int jMin = -1;

			// go trough all last
			int uDistance = 0;
			float sigmaDistance = 0;

			for (int j = 0; j < tracedMaxSignatures.size(); j++)
			{
				// calculate distances
				sigmaDistance = sigma - tracedMaxSignatures[j].back().sigma;
				uDistance = std::abs(cssPoint.u - tracedMaxSignatures[j].back().u);
				int maxUDistance = 0.2 * sigma + 10;
				
				// in case of closed contours use circular distance
				if ((uDistance > (length / 2)) && isClosed)
				{
					uDistance = length - uDistance;
				}
				
				// sigmaDistance is used to avoid false assignments in case of
				// overlapping curvy signatures
				if ((uDistance < lastMinUDistance) && (sigmaDistance <= maxSigmaDis) && (uDistance <= maxUDistance))
				{
					lastMinUDistance = uDistance;
					jMin = j;
				}
			}

			// store in vector to which distance was minimal
			// only assign if valid index found, otherwise throw away (or create new signature?)
			if (jMin != -1)
			{
				tracedMaxSignatures[jMin].push_back(cssPoint);
			}
			else
			{
				std::cout << "Warning: Could not trace local maximum (traceMinMaxSignatures) at sigma=" << sigma << std::endl;
				// TODO: maybe create new signature here?
				
				std::vector<CssPoint> maxSignature;
				maxSignature.push_back(cssPoint);
				tracedMaxSignatures.push_back(maxSignature);
				
			}
		}
	}
	/*********************************************************************/

	// do the same for the minima (this is just a copy of the upper part)
	/*********************************************************************/
	// Note: member variable vector uMinima only stores values for current sigma
	for (size_t i = 0; i < uMinima.size(); i++)
	{
		cssPoint.u = uMinima[i];
		cssPoint.sigma = sigma;
		cssPoint.kappa = kappa[cssPoint.u];
		
		// initialize if no previous values (just push back)
		if (sigma == sigmaSt)
		{
			// new vector for every point
			std::vector<CssPoint> minSignature;
			minSignature.push_back(cssPoint);
			tracedMinSignatures.push_back(minSignature);
			
			/*
			std::cout << "cssPoint.u = " << cssPoint.u << std::endl;
			std::cout << "cssPoint.sigma = " << cssPoint.sigma << std::endl;
			std::cout << "cssPoint.kappa = " << cssPoint.kappa << std::endl;
			*/
		}
		// otherwise find correspondences
		else
		{			
			// find next one in previous vectors, take last elements
			// note: not very effective because only signatures which are large enough must be checked
			int lastMinUDistance = 1000000;
			int jMin = -1;

			// go trough all last
			int uDistance = 0;
			float sigmaDistance = 0;

			for (int j = 0; j < tracedMinSignatures.size(); j++)
			{
				// calculate distances
				sigmaDistance = sigma - tracedMinSignatures[j].back().sigma;
				uDistance = std::abs(cssPoint.u - tracedMinSignatures[j].back().u);
				int maxUDistance = 0.2 * sigma + 10;
				
				// in case of closed contours use circular distance
				if ((uDistance > (length / 2)) && isClosed)
				{
					uDistance = length - uDistance;
				}
				
				// sigmaDistance is used to avoid false assignments in case of
				// overlapping curvy signatures
				if ((uDistance < lastMinUDistance) && (sigmaDistance <= maxSigmaDis) && (uDistance <= maxUDistance))
				{
					lastMinUDistance = uDistance;
					jMin = j;
				}
			}

			// store in vector to which distance was minimal
			// only assign if valid index found, otherwise throw away (or create new signature?)
			if (jMin != -1)
			{
				tracedMinSignatures[jMin].push_back(cssPoint);
			}
			else
			{
				std::cout << "Warning: Could not trace local minimum (traceMinMaxSignatures) at sigma=" << sigma << std::endl;
				// TODO: maybe create new signature here?
				std::vector<CssPoint> minSignature;
				minSignature.push_back(cssPoint);
				tracedMinSignatures.push_back(minSignature);
			}
		}
	}
	/*********************************************************************/
	
	/*
	if (sigma == sigmaStart)
	{
		std::cout << "tracedMaxSignatures.size() = " << tracedMaxSignatures.size() << std::endl;
		std::cout << "tracedMinSignatures.size() = " << tracedMinSignatures.size() << std::endl;
		std::cout << "sigma = " << sigma << std::endl;
	}
	*/
}

float CSSApproximation::calcThreshold(float sigma)
{
	float val = 0.2807 * std::pow(sigma, -0.867);
	return val;
}

void CSSApproximation::deleteShortSignatures(int minLength)
{
	int cntMin = 0;
	int cntMax = 0;
	
	// minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		if (tracedMinSignatures[i].size() < minLength)
		{
			// myvector.erase (myvector.begin()+5); = erase 6th element
			tracedMinSignatures.erase(tracedMinSignatures.begin()+i);
			i--;

			cntMin++;
		}
	}

	// maxima
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		if (tracedMaxSignatures[i].size() < minLength)
		{
			// myvector.erase (myvector.begin()+5); = erase 6th element
			tracedMaxSignatures.erase(tracedMaxSignatures.begin()+i);
			i--;

			cntMax++;
		}
	}
	
	std::cout << "Deleted: MinSignatures: " << cntMin++ << " | MaxSignatures: " << cntMax << std::endl;
}

void CSSApproximation::calculateShapeFeatures(Edge* e)
{
	
	// at this point, tracing is finished, extract contour fragments
	std::vector<cv::Point> fragmentPoints;
	std::vector<cv::Point> p = e->getEdgePoints();
	
	// process minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// delete old fragmentPoints
		fragmentPoints.clear();
		
		// get x, y coordinates of local extremum
		int uCoord = tracedMinSignatures[i].at(0).u;
		cv::Point XYCoordCenter(p[uCoord].x, p[uCoord].y);
		
		// get curvature of local extremum
		float kappa = tracedMinSignatures[i].at(0).kappa;
		
		// determine current radius
		float radius = sF * tracedMinSignatures[i].back().sigma;
				
		// now extract the fragment
		// check points left from and including local extremum, break when beginning reached
		for (int j = uCoord; j > -100000; j--)
		{
			// index correction in case of closed contours
			if (j < 0 && e->getIsLoop())
			{
				j = j + p.size();
			}
			
			// do not exceed limits (j >= 0) for open contours
			if ((e->getDistance(XYCoordCenter, p[j]) <= radius) && j >= 0)
			{
				fragmentPoints.push_back(p[j]);
			}
			else
			{
				break;
			}
		}
		
		// inverse elements of vector
		std::reverse(fragmentPoints.begin(), fragmentPoints.end());

		// determine u-Coordinate of extremum in local fragment
		int uCoordFrag = fragmentPoints.size()-1;
		
		// check points right from local extremum, cancel when radius reached
		for (int j = uCoord+1; j < 100000; j++)
		{
			// index correction in case of closed contours
			if (j >= p.size() && e->getIsLoop())
			{
				j = j - p.size();
			}
			
			// TODO: avoid case when radius captures complete closed contour
			// in this case, no criterion for cancleing is reached
			
			// do not exceed limits (j < p.size()) for open contours
			if ((e->getDistance(XYCoordCenter, p[j]) <= radius) && j < p.size())
			{
				fragmentPoints.push_back(p[j]);
			}
			else
			{
				break;
			}
		}

		//std::cout << "min " << i << ":" << std::endl;
		//std::cout << "kappa = " << kappa << std::endl;
		
		// determine and save rotation for the fragment
		float rotMin = shapeFeat_t.calcOrientation(fragmentPoints, uCoordFrag, radius, kappa);
		rotMinSignatures.push_back(rotMin);
	}

	// process maxima (copy of function from above)
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// delete old fragmentPoints
		fragmentPoints.clear();
		
		// get x, y coordinates of local extremum
		int uCoord = tracedMaxSignatures[i].at(0).u;
		cv::Point XYCoordCenter(p[uCoord].x, p[uCoord].y);
		
		// get curvature of local extremum
		float kappa = tracedMaxSignatures[i].at(0).kappa;
		
		// determine current radius
		float radius = sF * tracedMaxSignatures[i].back().sigma;

		// check points left from and including local extremum, break when beginning reached
		for (int j = uCoord; j > -100000; j--)
		{
			// index correction in case of closed contours
			if (j < 0 && e->getIsLoop())
			{
				j = j + p.size();
			}
			
			// do not exceed limits (j >= 0) for open contours
			if ((e->getDistance(XYCoordCenter, p[j]) <= radius) && j >= 0)
			{
				fragmentPoints.push_back(p[j]);
			}
			else
			{
				break;
			}
		}
		
		// inverse elements of vector
		std::reverse(fragmentPoints.begin(), fragmentPoints.end());
		
		// determine u-Coordinate of extremum in local fragment
		int uCoordFrag = fragmentPoints.size()-1;
		
		// check points right from local extremum, cancel when radius reached
		for (int j = uCoord+1; j < 100000; j++)
		{
			// index correction in case of closed contours
			if (j >= p.size() && e->getIsLoop())
			{
				j = j - p.size();
			}
			
			// do not exceed limits (j < p.size()) for open contours
			if ((e->getDistance(XYCoordCenter, p[j]) <= radius) && j < p.size())
			{
				fragmentPoints.push_back(p[j]);
			}
			else
			{
				break;
			}
		}

		//std::cout << "max " << i << ":" << std::endl;
		//std::cout << "kappa = " << kappa << std::endl;

		// determine and save rotation for the fragment
		float rotMax = shapeFeat_t.calcOrientation(fragmentPoints, uCoordFrag, radius, kappa);
		rotMaxSignatures.push_back(rotMax);
	}
}

CSSApproximation::~CSSApproximation()
{
	// delete all elements
	uMinima.clear();
	uMaxima.clear();
	
	tracedMaxSignatures.clear();
	tracedMinSignatures.clear();
	
	rotMinSignatures.clear();
	rotMaxSignatures.clear();
}



