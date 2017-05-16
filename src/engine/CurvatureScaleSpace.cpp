#include "CurvatureScaleSpace.h"
#include <math.h>
#include <iostream>
//for GUI stuff can be deleted if no imshow() is used
//#include <opencv2/highgui/highgui.hpp>
//using namespace cv;

#include "./../util/Visualization.h"
Visualization visualization;

#include "./../util/ShapeFeatures.h"
ShapeFeatures shapeFeatures;

float lengthFactor = 4.5; // lengthFactor * sigma samples at each side of max of gaussian
float sigmaMax = 200.0;
float sigmaStep = 2.0; // default 2.0
float sigmaStart = 8.0; // default 8.0

float maxSigmaDistance = 3.0 * sigmaStep; // tracing of extrema signatures (see traceMinMaxSignatures)
const float sF = 0.3; // scale factor for radius
int minSigmaSteps = 8; // delete signatures shorter than N sigma steps, default 8

float testScale = 100;

CurvatureScaleSpace::CurvatureScaleSpace()
{
	//std::cout << "CurvatureScaleSpace object created" << std::endl;
}

void CurvatureScaleSpace::calculateScaleSpaces(std::vector<Edge*>* v)
{
	// delete old tracing results
	tracedMinSignatures.clear();
	tracedMaxSignatures.clear();
	
	rotMinSignatures.clear();
	rotMaxSignatures.clear();
	
	// calculate and analyse scale spaces for all edges
	for (std::vector<Edge*>::iterator it = v->begin(); it != v->end(); ++it)
	{
		calculateCSSInterestPoints((*it));
	}
}

// the following function is called for every edge e
void CurvatureScaleSpace::calculateCSSInterestPoints(Edge* e)
{
	float sigma = sigmaStart;
	int eLength = int(e->length()); // length of edge
	//std::cout << "eLength " << eLength << std::endl;

	// test: invert edge
	//e->invertRecursive();

	// data structures for kappa and derivatives
	float _kappa[eLength];
	std::vector<float> _xCoorU;
	std::vector<float> _xCoorUU;
	std::vector<float> _yCoorU;
	std::vector<float> _yCoorUU;

	// Get X and Y values seperately
	std::vector<int> X;
	std::vector<int> Y;
	X = e->getEdgeXCoordinatesVector();
	Y = e->getEdgeYCoordinatesVector();

	int numberSigmaValues = 0; // just a counter
	while (sigma <= sigmaMax)
	{
		//std::cout << "sigma = " << sigma << std::endl;

		std::vector<float> _gauss1st;
		std::vector<float> _gauss2nd;

		// calculate 1st and 2n derivative of Gaussian
		gauss1st(sigma, &_gauss1st, true);
		gauss2nd(sigma, &_gauss2nd, true);

		//clock_t t1 = clock();
		if(!e->getIsLoop())
		{
			convolve(&X, &_gauss1st, &_xCoorU);
			convolve(&Y, &_gauss1st, &_yCoorU);
			convolve(&X, &_gauss2nd, &_xCoorUU);
			convolve(&Y, &_gauss2nd, &_yCoorUU);
		}
		else
		{	// do this of loop is closed
			convolveClosedContour4(&X, &_gauss1st, &_xCoorU);
			convolveClosedContour4(&Y, &_gauss1st, &_yCoorU);
			convolveClosedContour4(&X, &_gauss2nd, &_xCoorUU);
			convolveClosedContour4(&Y, &_gauss2nd, &_yCoorUU);
		}
		//clock_t t2 = clock();
		//std::cout << double(t2 - t1) / CLOCKS_PER_SEC << std::endl;	
		
		// calculate curvature _kappa for current sigma
		calcKappa(&_xCoorU, &_yCoorU, &_xCoorUU, &_yCoorUU, eLength, _kappa);

		// search for local minima and maxima and store values in uMinima and uMaxima
		// uMinima and uMaxima are member variables for temporary storing the results
		findLocalMinMax(_kappa, eLength, sigma, e->getIsLoop());
		traceMinMaxSignatures(_kappa, eLength, sigma, e->getIsLoop());

		/*********************************************************************/
		// Note: this part is only for development purposes
		if (sigma == testScale)
		{
			// calculate corresponding smoothed curve for visualization
			std::vector<float> _gauss;
			std::vector<float> _xSmoothed;
			std::vector<float> _ySmoothed;
			
			gauss(sigma, &_gauss, true);

			if (!e->getIsLoop())
			{
				convolve(&X, &_gauss, &_xSmoothed);
				convolve(&Y, &_gauss, &_ySmoothed);
			}
			else
			{
				convolveClosedContour4(&X, &_gauss, &_xSmoothed);
				convolveClosedContour4(&Y, &_gauss, &_ySmoothed);
			}
	
			// convert to Edge object for drawing gauss filtered contour
			std::vector<cv::Point> edgePoints;
			cv::Point p;
			for (int i = 0; i < eLength; i++)
			{
				p.x = std::round(_xSmoothed[i]);
				p.y = std::round(_ySmoothed[i]);
				//std::cout << "p = (" << _xSmoothed[i] << ", " << _ySmoothed[i] << ")" << std::endl;
				
				edgePoints.push_back(p);
			}
			Edge edgeGauss(edgePoints);

			// write SVG file of colored smoothed curve
			visualization.saveSVGcolorWheel(&edgeGauss, sigma);

			// write Matlab file		
			visualization.saveMatlabKappa(_kappa, eLength, sigma, uMinima, uMaxima, thresLocalMinMax);
			std::cout << "sigmaPlotMatlab = " << sigma << std::endl;
			std::cout << "thresLocalMinMax = " << thresLocalMinMax << std::endl;
		}
		/*********************************************************************/

		// the following two functions are called for every sigma
		// new values are therefore attached to the files
		//visualization.saveCSVKappa(_kappa, eLength, sigma);
		visualization.saveCSVMinMax(&uMinima, &uMaxima, sigma);

		// display current status of convolution
		printConvolutionStatus(sigma);

		sigma += sigmaStep;
		numberSigmaValues++;
	}

	std::cout << "tracedMinSignatures.size() = " << tracedMinSignatures.size() << std::endl;
	std::cout << "tracedMaxSignatures.size() = " << tracedMaxSignatures.size() << std::endl;
	
	// delete signatures shorter than N sigma steps
	deleteShortSignatures(minSigmaSteps);
	
	// Tracing is finished here - tracedMin/MaxSignatures are cleared for every edge
	calculateShapeFeatures(e);

	visualization.generateColorWheel(eLength);
	visualization.saveCSVTracedMinMax(tracedMinSignatures, tracedMaxSignatures);
	visualization.saveSvgCssResult(e, tracedMinSignatures, tracedMaxSignatures, rotMinSignatures, rotMaxSignatures, sF);
	
	// output for calculation threshold function
	int no = 3;
	int sz = tracedMaxSignatures[no].size();
	
	std::cout << "kappa = [";
	for (int i = 0; i < sz; i++)
	{
		std::cout << std::fabs(tracedMaxSignatures[no].at(i).kappa);
		if (i != (sz-1)) { std::cout << ", "; }
	}
	std::cout << "];" << std::endl;
	
	std::cout << "scale = [";
	for (int i = 0; i < sz; i++)
	{
		std::cout << tracedMaxSignatures[no].at(i).sigma;
		if (i != (sz-1)) { std::cout << ", "; }
	}
	std::cout << "];" << std::endl;
}

inline void CurvatureScaleSpace::calcKappa(
		std::vector<float>* xCoorU,
		std::vector<float>* yCoorU,
		std::vector<float>* xCoorUU,
		std::vector<float>* yCoorUU,
		int length, float* output)
{
	for (int i = 0; i < length; i++)
	{
		// old version with checks (at commands)
		output[i] = (xCoorU->at(i) * yCoorUU->at(i) - xCoorUU->at(i) * yCoorU->at(i))
			/ pow(pow(xCoorU->at(i), 2.0) + pow(yCoorU->at(i), 2.0), 1.5); // 3/2 = 1.5

		// new version without checks
		/*
		output[i] = ((*xCoorU)[i] * (*yCoorUU)[i] - (*xCoorUU)[i] * (*yCoorU)[i])
			/ pow(pow((*xCoorU)[i], 2.0) + pow((*yCoorU)[i], 2.0), 1.5); // 3/2 = 1.5
		*/
	}
}

void CurvatureScaleSpace::gauss(float sigma, std::vector<float>* output, bool normalized)
{
	// noPts = number of samples at a side of max
	int noPts = ceil(lengthFactor*sigma);

	output->clear();
	output->reserve(2*noPts+1);

	float sum = 0;
	float tmp;
	
	for (int i = 0; i <= 2*noPts; i++)
	{
		tmp = exp(-pow((i-noPts) / sigma, 2) / 2) / SQRT_2_PI / sigma;
		output->push_back(tmp);

		//sum += fabsf(output->at(i));
		sum += fabsf((*output)[i]);
	}

	if (normalized)
	{
		for (int i = 0; i <= 2*noPts; i++)
		{
			(*output)[i] /= sum;
		}
	}
}

void CurvatureScaleSpace::gauss1st(float sigma, std::vector<float>* output, bool normalized)
{
	// noPts = number of samples at each side of max
	int noPts = ceil(lengthFactor*sigma);

	output->clear();
	output->reserve(2*noPts+1);

	float sum = 0;
	float tmp;

	for (int i = 0; i <= 2*noPts; i++)
	{
		tmp = exp(-pow((i-noPts) / sigma, 2) / 2) * (i-noPts) / SQRT_2_PI / pow(sigma, 3);
		output->push_back(tmp);
		
		//sum += fabsf(output->at(i));
		sum += fabsf((*output)[i]);
	}

	if (normalized)
	{
		for (int i = 0; i <= 2*noPts; i++)
		{
			//output->at(i) /= sum;
			(*output)[i] /= sum;
		}
	}
	
	/*
	// print data for plot
	std::cout << "i = [";
	for (int i = 0; i <= 2*noPts; i++)
	{
		std::cout << i << " "; 
	}
	std::cout << "];" << std::endl;
	
	std::cout << "g = [";
	for (int i = 0; i <= 2*noPts; i++)
	{
		std::cout << (*output)[i] << " "; 
	}
	std::cout << "];" << std::endl;
	*/
}

void CurvatureScaleSpace::gauss2nd(float sigma, std::vector<float>* output, bool normalized)
{
	// noPts = number of samples at each side of max
	int noPts = ceil(lengthFactor*sigma);

	output->clear();
	output->reserve(2*noPts+1);

	float sum = 0;
	float tmp;

	for (int i = 0; i <= 2*noPts; i++)
	{
		tmp = exp(-pow((i-noPts) / sigma, 2) / 2) * (pow((i-noPts), 2) - pow(sigma, 2)) / SQRT_2_PI / pow(sigma, 5);
		output->push_back(tmp);
		
		//sum += fabsf(output->at(i));
		sum += fabsf((*output)[i]);
	}

	if (normalized)
	{
		for (int i = 0; i <= 2*noPts; i++)
		{
			//output->at(i) /= sum;
			(*output)[i] /= sum;
		}
	}
}

void CurvatureScaleSpace::convolve(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res)
{
	const int sigSize = sig->size();
	const int kernSize = kern->size();
	const int padSize = (kernSize-1) / 2;
		
	// input and output should have same size
	res->resize(sigSize);

	// padding
	std::vector<float> paddedSig;
	paddedSig.clear();
	//paddedSig.reserve(2*padSize+sigSize);

	// pad beginning of signal
	for (int i = 0; i < padSize; i++)
	{
		paddedSig.push_back((*sig)[0]);
		//paddedSig.push_back(0); 
	}
	// write signal
	for (int i = 0; i < sigSize; i++)
	{
		paddedSig.push_back((*sig)[i]); 
	}
	// pad end of signal
	for (int i = 0; i < padSize; i++)
	{
		paddedSig.push_back((*sig)[sigSize-1]);
		//paddedSig.push_back(0);
	}

	// convolution
	for (int i = 0; i < sigSize; i++)
	{
		(*res)[i] = 0;
		for (int k = 0; k < kernSize; k++)
		{
			//res->at(i) = res->at(i) + (paddedSig.at(i+k) * kern->at(k));
			(*res)[i] = (*res)[i] + (paddedSig[i+k] * (*kern)[k]);
		}
	}
	
	//std::cout << "size padded sig = " << paddedSig.size() << std::endl;
}

void CurvatureScaleSpace::convolveClosedContour1(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res)
{
	//TODO: check commutativity for discrete convolution

	const int sigSize = sig->size();
	const int kernSize = kern->size();
	const int padSize = (kernSize-1) / 2;

	// input and output should have same size
	res->resize(sigSize);

	// convolution
	int idx = 0;
	for (int i = 0; i < sigSize; i++)
	{
		(*res)[i] = 0;
		for (int k = 0; k < kernSize; k++)
		{
			// Version Daniel
			//(*res)[i] = (*res)[i] + ((*sig)[(i+k) % sigSize] * (*kern)[k]);
			
			// Version Markus
			// beginning => jump to end
			if ((i < padSize) && (k < padSize))
			{
				idx = sigSize-1-k;
			}
			// end => jump to beginning
			else if ((i > sigSize-padSize-1) && (k > padSize))
			{
				idx = k; 
			}
			else
			{
				idx = i+k;
			}
			
			(*res)[i] = (*res)[i] + ((*sig)[idx] * (*kern)[k]);
		}
	}
}

void CurvatureScaleSpace::convolveClosedContour2(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res)
{
	const int sigSize = sig->size();
	const int kernSize = kern->size();
	const int padSize = (kernSize-1) / 2;

	// check for failures
	if (padSize > sigSize)
	{
		std::cout << "Error: padSize > sigSize in CurvatureScaleSpace::convolveClosedContour2" << std::endl;
	}

	// input and output should have same size
	res->resize(sigSize);

	// padding
	std::vector<float> paddedSig;
	paddedSig.clear();
	//paddedSig.reserve(2*padSize+sigSize);

	// pad beginning of signal
	for (int i = 0; i < padSize; i++)
	{
		paddedSig.push_back((*sig)[sigSize-padSize+i]);
		//paddedSig.push_back(0);
	}
	// write signal
	for (int i = 0; i < sigSize; i++)
	{
		paddedSig.push_back((*sig)[i]); 
	}
	// pad end of signal
	for (int i = 0; i < padSize; i++)
	{
		paddedSig.push_back((*sig)[i]);
		//paddedSig.push_back(0);
	}

	// convolution
	for (int i = 0; i < sigSize; i++)
	{
		(*res)[i] = 0;
		for (int k = 0; k < kernSize; k++)
		{
			//res->at(i) = res->at(i) + (paddedSig.at(i+k) * kern->at(k));
			(*res)[i] = (*res)[i] + (paddedSig[i+k] * (*kern)[k]);
		}
	}
	
	//std::cout << "size padded sig = " << paddedSig.size() << std::endl;
}

inline int positive_modulo(int i, int n)
{
    return (i % n + n) % n;
}

void CurvatureScaleSpace::convolveClosedContour3(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res)
{
	const int sigSize = sig->size();
	const int kernSize = kern->size();
	const int padSize = (kernSize-1) / 2;
	
	/*
	std::cout << "padSize = " << padSize << std::endl;
	std::cout << "kernSize = " << kernSize << std::endl;
	std::cout << "sigSize = " << sigSize << std::endl;
	*/

	// input and output should have same size
	res->resize(sigSize);

	// convolution
	int index = 0;
	for (int i = 0; i < sigSize; i++)
	{
		(*res)[i] = 0;

		for (int k = 0; k < kernSize; k++)
		{
			index = positive_modulo(i+k-padSize, sigSize);
			(*res)[i] = (*res)[i] + ((*sig)[index] * (*kern)[k]);
			//(*res)[i] = (*res)[i] + (sig->at(index) * (*kern)[k]);
		}
	}
}

void CurvatureScaleSpace::convolveClosedContour4(std::vector<int>* sig, std::vector<float>* kern, std::vector<float>* res)
{
	const int sigSize = sig->size();
	const int kernSize = kern->size();
    const int padSize = (kernSize-1) / 2;
    const int maxPadSize = ceil(lengthFactor*sigmaMax);

	static std::vector<float> paddedSig; // static duration via static keyword. This line is only executed once.
	res->resize(sigSize); // input and output should have same size
	
	// padded signal must only be recalculated for a new signal (= new edge = new adress)
	static int *lastAdress = NULL; // this line is only executed once

	if (lastAdress != &(*sig)[0]) // entered in case of new signal
	{
		// clean up paddedSig
		paddedSig.clear();
		int index = 0;
		
		// pad beginning of signal
		for (int i = 0; i < maxPadSize; i++)
		{
			index = positive_modulo(sigSize-maxPadSize+i, sigSize);
			paddedSig.push_back((*sig)[index]);
			//std::cout << "index beginning = " << index << std::endl;
		}
		// write signal
		for (int i = 0; i < sigSize; i++)
		{
			paddedSig.push_back((*sig)[i]);
		}
		// pad end of signal
		for (int i = 0; i < maxPadSize; i++)
		{
			index = positive_modulo(sigSize+i, sigSize);
			paddedSig.push_back((*sig)[index]);
			//std::cout << "index end = " << index << std::endl;
		}
		
		// save adress of signal (reason: see above)
		lastAdress = &(*sig)[0];
	}
	
	// convolution
	int offset = maxPadSize-padSize;
	for (int i = 0; i < sigSize; i++)
	{
		(*res)[i] = 0;
		for (int k = 0; k < kernSize; k++)
		{
			//res->at(i) = res->at(i) + (paddedSig.at(i+k) * kern->at(k));
			(*res)[i] = (*res)[i] + (paddedSig[i+k+offset] * (*kern)[k]);
			//std::cout << "i+k+offset = " << (i+k+offset) << std::endl;
		}
		//std::cout << "==========================" << std::endl;
	}
}

void CurvatureScaleSpace::findLocalMinMax(float* kappa, int length, float sigma, bool isClosed)
{
	// empty the vectors for every new calculation
	uMinima.clear();
	uMaxima.clear();

	// mean value of kappa
	float kappaAbsMean = 0;
	float kappaAbsVariance = 0;

	// mean
	/*
	for (int i = 0; i < length; i++)
	{
		kappaAbsMean = kappaAbsMean + std::abs(kappa[i]);
	}
	kappaAbsMean = kappaAbsMean / length;

	// variance
	for (int i = 0; i < length; i++)
	{
		kappaAbsVariance = kappaAbsVariance + ((kappa[i] - kappaAbsMean) * (kappa[i] - kappaAbsMean));
	}
	kappaAbsVariance = kappaAbsVariance / (length-1);
	*/

	//std::cout << "kappaMean = " << kappaAbsMean << std::endl;
	//std::cout << "kappaVariance = " << kappaVariance << std::endl;
	
	// find local maxima and minima
	// algorithm is inspired by http://www.billauer.co.il/peakdet.html (last checked: 2016-03-30)
	
	const float thres = 0.002;
	//float thres = calcThreshold(sigma);
	//const float thres = 0.008;
	//const float thres = 0.03 * (1 - (sigma * 1.0 / sigmaMax));
	//const float thres = 0.01 * (sigma * 1.0 / sigmaMax);
	//const float delta = 3.0 * kappaAbsVariance;
	const float delta = 0.001;
	//std::cout << "delta = " << delta << std::endl;
	
	// make threshold available outside this function (e. g. for plot)
	thresLocalMinMax = thres;
	
	float val;
	int maxPos, minPos;
	float min = 1000000; float max = -1000000;

	// do two iterations in case of closed contours so that extrema at
	// both back and front are considered
	int nuIterations = 1;
	if (isClosed) { nuIterations = 2; };

	// first iteration is for initializing min, max, and lookformax
	// values are only stored in second run
	bool lookForMax = true;
	for (int h = 0; h < nuIterations; h++)
	{
		for (int i = 0; i < length; i++)
		{
			val = kappa[i];

			// current value larger then last max?
			// this part saves last highest value
			if (val > max)
			{
				max = val;
				maxPos = i;
			}
			// current value smaller then last min?
			// this part saves last smallest value
			if (val < min)
			{
				min = val;
				minPos = i;
			}
			
			if (lookForMax) // true for i = 0, j = 0
			{
				if (val < (max - delta))
				{
					// only take maxima higher then certain threshold	
					// if closed: only entered in second run (h > 0 || !isClosed)
					bool atBorder = !isClosed && (maxPos == 0 || maxPos == length-1);
					if ((max > thres) && (h > 0 || !isClosed) && !atBorder)
					{
						uMaxima.push_back(maxPos);
					}
					
					// save current value for minimum search
					//uMaxima.push_back(maxPos);
					min = val; minPos = i;
					lookForMax = false;
				}
			}
			else
			{
				if (val > (min + delta))
				{
					// only take minima lower then certain threshold
					// if closed: only entered in second run (h > 0 || !isClosed)
					bool atBorder = !isClosed && (minPos == 0 || minPos == length-1);
					if (min < (-1.0 * thres) && (h > 0 || !isClosed) && !atBorder)
					{
						uMinima.push_back(minPos);
					}
					// save current value for maximum search
					//uMinima.push_back(minPos);
					max = val; maxPos = i;
					lookForMax = true;
				}
			}
		}
	}
}

float CurvatureScaleSpace::calcThreshold(float sigma)
{
	float val = 0.2807 * std::pow(sigma, -0.867);  // first variant

	//std::cout << "val = " << val << std::endl;
	return val;
}

// static int value = 1; // static duration via static keyword. This line is only executed once.
// TODO: check if this also works for multiple edges
void CurvatureScaleSpace::printConvolutionStatus(float currentSigma)
{
	// calculate number of sigmas (approximately)
	static int numberSigmas = (sigmaMax-sigmaStart) / sigmaStep;
	
	// save current number of sigma
	static int numberCurrentSigma = 0;
	numberCurrentSigma++;
	
	// calculate current percentage value
	int currentPercent = 100.0 * numberCurrentSigma / (1.0 * numberSigmas);
	float correctionFactor = currentPercent / 100.0; // assume linear function
	
	float tmp = currentPercent;
	currentPercent = tmp * correctionFactor;

	// only print new percent values
	static int lastPercent = currentPercent;
	if (currentPercent > lastPercent || numberCurrentSigma == 1)
	{ 
		std::cout << "Status: " << currentPercent << "%\r" << std::flush;
		lastPercent = currentPercent;
	}
}

void CurvatureScaleSpace::traceMinMaxSignatures(float* kappa, int length, float sigma, bool isClosed)
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
	// Note: member variable vector uMaxima only stores values for current sigma
	for (size_t i = 0; i < uMaxima.size(); i++)
	{
		cssPoint.u = uMaxima[i];
		cssPoint.sigma = sigma;
		cssPoint.kappa = kappa[cssPoint.u];
		
		// initialize if no previous values (just push back)
		if (sigma == sigmaStart)
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
				if ((uDistance < lastMinUDistance) && (sigmaDistance <= maxSigmaDistance) && (uDistance <= maxUDistance))
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
				std::cout << "Warning: Could not trace local maximum (traceMinMaxSignatures)" << std::endl;
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
		if (sigma == sigmaStart)
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
				if ((uDistance < lastMinUDistance) && (sigmaDistance <= maxSigmaDistance) && (uDistance <= maxUDistance))
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
				std::cout << "Warning: Could not trace local minimum (traceMinMaxSignatures)" << std::endl;
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

void CurvatureScaleSpace::deleteShortSignatures(int minLength)
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

void CurvatureScaleSpace::calculateShapeFeatures(Edge* e)
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
		float rotMin = shapeFeatures.calcOrientation(fragmentPoints, uCoordFrag, radius, kappa);
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
		float rotMax = shapeFeatures.calcOrientation(fragmentPoints, uCoordFrag, radius, kappa);
		rotMaxSignatures.push_back(rotMax);
	}
}

CurvatureScaleSpace::~CurvatureScaleSpace()
{
	// delete all elements
	uMinima.clear();
	uMaxima.clear();
	
	tracedMaxSignatures.clear();
	tracedMinSignatures.clear();
	
	rotMinSignatures.clear();
	rotMaxSignatures.clear();
}
