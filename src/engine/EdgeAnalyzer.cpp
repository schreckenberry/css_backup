#include "EdgeAnalyzer.h"
#include <iostream>

float sigmaMax = 120.0;
float sigmaStep = 0.2;
float threshMaxMin = 0.2;
float threshPtsMax = 0.3;
float threshPtsMin = 0.3;
float threshPtsZero = 0.15;

void EdgeAnalyzer::calculateScaleSpaces(std::vector<Edge*>* v,
		std::vector<InterestPoint*>* maximaVector,
		std::vector<InterestPoint*>* minimaVector,
		std::vector<InterestPoint*>* zeroVector) {
	//TODO set max parameter here
	// Ablauf:
	// 1. Unterabtastung (resizeEdge)
	// 2. calcInterestPoints

	std::vector<Edge*> resizedEdges;
	for (std::vector<Edge*>::iterator it = v->begin(); it != v->end(); ++it) {
		if ((*it)->length() > 200) {
			resizedEdges.push_back(
					resizeEdge((*it), 200)
					); // resampling wenn kante zu viele punkte hat (it ist eine Kante)

		} else if ((*it)->length() > 20) {
			resizedEdges.push_back(
					resizeEdge((*it), int((*it)->length()/2))
					);
		} else {
			resizedEdges.push_back(
					(*it)
					);
		}
		//std::cout << "Length before " << (*it)->length() << " Length after " << resizedEdges[resizedEdges.size()-1]->length() << std::endl;
	}

	// berechne für alle kanten die skalenräume
	for (std::vector<Edge*>::iterator it = resizedEdges.begin(); it != resizedEdges.end(); ++it) {
		calculateInterestPoints((*it), maximaVector, minimaVector, zeroVector);
	}
}

void EdgeAnalyzer::calculateInterestPoints(Edge* e,
		std::vector<InterestPoint*>* maximaVector,
		std::vector<InterestPoint*>* minimaVector,
		std::vector<InterestPoint*>* zeroVector){
	float sigma = 1.0;
	int eLength = int(e->length()); // length of edge
	int tail, lengthGauss;
	//std::cout << "eLength " << eLength << std::endl;
	// CSS
	float _kappa[eLength];
	float _xCoorU[eLength];
	float _xCoorUU[eLength];
	float _yCoorU[eLength];
	float _yCoorUU[eLength];

	// Maxima
	std::vector< std::vector<int> > maxima; // scale space for local maxima of curvature
	std::vector<int> maximaPts; // indices of interesting points of maxima scale space
	std::vector<float> maximaSigma; // sigma of interesting points of maxima scale space

	// Minima
	std::vector< std::vector<int> > minima;
	std::vector<int> minimaPts;
	std::vector<float> minimaSigma;

	// Zero Crossings
	std::vector< std::vector<int> > zero;
	std::vector<int> zeroPts;
	std::vector<float> zeroSigma;

	int numberConvs = 0;

	while (sigma < sigmaMax) {
		tail = ceil(4.5 * sigma); // calculate sampled size of kernel
		lengthGauss = 2*ceil(4.5 * sigma)+1; // length of gauss function

		float _gauss1st[lengthGauss];
		float _gauss2nd[lengthGauss];

		// calculate curvature of edge
		gauss1st(tail, sigma, _gauss1st, 1);
		gauss2nd(tail, sigma, _gauss2nd, 1);

		int X[eLength];
		int Y[eLength];
		e->getEdgeXCoordinates(X);
		e->getEdgeYCoordinates(Y);

		//clock_t t1 = clock();
		//
		convolution(X, _xCoorU, eLength, _gauss1st, lengthGauss);
		convolution(Y, _yCoorU, eLength, _gauss1st, lengthGauss);
		convolution(X, _xCoorUU, eLength, _gauss2nd, lengthGauss);
		convolution(Y, _yCoorUU, eLength, _gauss2nd, lengthGauss);
		//
		//clock_t t2 = clock();
		//std::cout << double(t2 - t1) / CLOCKS_PER_SEC << std::endl;
		
		numberConvs++;
		
/*
 *  def aus header:
	inline void calcKappa(xCoorU, yCoorU, xCoorUU, yCoorUU, int length, float* output);
*/
		calcKappa(_xCoorU, _yCoorU, _xCoorUU, _yCoorUU, eLength, _kappa);

		// find maxima, minima and zero crossings
		maxima.push_back(findLocalMaxima(_kappa, eLength, threshMaxMin));
		minima.push_back(findLocalMinima(_kappa, eLength, threshMaxMin));
		zero.push_back(findZeroCrossings(_kappa, eLength));

		// exit condition
		if (zero.back().size() < 2) {
			break;
		}
		sigma += sigmaStep;
	}
	std::cout << "numberConvs = " << numberConvs++ << std::endl;

	getMaximaMax(&maximaPts, &maximaSigma, &maxima, threshPtsMax);
	getMinimaMax(&minimaPts, &minimaSigma, &minima, threshPtsMin);
	getZeroMax(&zeroPts, &zeroSigma, &zero, threshPtsZero);

	//std::cout << "maximaPts.size " << maximaPts.size() << " minimaPts.size " << minimaPts.size() << " zeroPts.size " << zeroPts.size() << std::endl;

	std::vector<cv::Point> points = e->getEdgePoints();
	for(int i = 0; i < maximaPts.size(); i++){
		maximaVector->push_back(new InterestPoint(points[ maximaPts[i] ]));
	}
	for(int i = 0; i < minimaPts.size(); i++){
		minimaVector->push_back(new InterestPoint(points[ minimaPts[i] ]));
	}
	for(int i = 0; i < zeroPts.size(); i++){
		zeroVector->push_back(new InterestPoint(points[ zeroPts[i] ]));
	}
}

std::vector<int> EdgeAnalyzer::findLocalMaxima(float * input, int length, float threshold) {
    std::vector<int> maximaIdx;
    maximaIdx.reserve(50);
    int counter = 0;
    float currMax = 0;

    for (int i = 3; i < length-2; i++) {
        currMax = max(input+i-3, 7);
        if (input[i] == currMax && fabs(input[i]) > threshold) {
            maximaIdx.push_back(i);
            counter++;
        }
    }
    maximaIdx.resize(counter); //shrink size if allocated to much
    return maximaIdx;
}

std::vector<int> EdgeAnalyzer::findLocalMinima(float * input, int length, float threshold){
    std::vector<int> minimaIdx;
    minimaIdx.reserve(50);
    float currMin = 0;
    int counter = 0;

    for (int i = 3; i < length-2; i++) {
        currMin = min(input+i-3, 7);
        if (input[i] == currMin && fabs(input[i]) > threshold) {
            minimaIdx.push_back(i);
            counter++;
        }
    }
    minimaIdx.resize(counter); //shrink size if allocated to much
    return minimaIdx;
}

std::vector<int> EdgeAnalyzer::findZeroCrossings(float * input, int length)
{
    std::vector<int> zeroIdx;
    zeroIdx.reserve(50);
    int counter = 0;

    for (int i = 1; i < length; i++) {
        if ((input[i] > 0 && input[i-1] < 0) || (input[i] < 0 && input[i-1] > 0)) {
            zeroIdx.push_back(i);
            counter++;
        }
    }
    zeroIdx.resize(counter);
    return zeroIdx;
}

void EdgeAnalyzer::getMaximaMax(std::vector<int>* maximaPts,
		std::vector<float>* maximaSigma,
		std::vector< std::vector<int> >* maxima,
		float threshold)
{
    int len = maxima->size();
    int i = len-1;

    while (i > 0) {
        if ((*maxima)[i].size() > 0) {
            break;
        }
        i--;
    }

    float maxSigma = i * sigmaStep;
    if (maxSigma < 2)
        return;

    float minSigma = round(maxSigma * threshold);
    if (minSigma < 1)
        minSigma = 1;

    int currLen = 0;
    int prevLen = 0;
    while (i > 3) {
        i--;
        if (i * sigmaStep < minSigma)
            break;

        std::vector<int> indicesStart = (*maxima)[i];

        if (currLen > prevLen)
            prevLen = currLen;

        currLen = (int)indicesStart.size();

        if (currLen > prevLen) {
            int j = i;
            while (j > 2) {
                std::vector<int> indicesLow = (*maxima)[i];
                findOriginMax(maximaPts, maximaSigma, indicesLow, indicesStart, i);

                j = j - 10;
            }
        }
    }
}

void EdgeAnalyzer::getMinimaMax(std::vector<int>* minimaPts,
		std::vector<float>* minimaSigma,
		std::vector< std::vector<int> >* minima,
		float threshold)
{
    int len = minima->size();
    int i = len-1;

    while (i > 0) {
        if ((*minima)[i].size() > 0) {
            break;
        }
        i--;
    }

    float maxSigma = i * sigmaStep;
    if (maxSigma < 2)
        return;

    float minSigma = round(maxSigma * threshold);
    if (minSigma < 1)
        minSigma = 1;

    int currLen = 0;
    int prevLen = 0;
    while (i > 3) {
        i--;
        if (i * sigmaStep < minSigma)
            break;

        std::vector<int> indicesStart = (*minima)[i];

        if (currLen > prevLen)
            prevLen = currLen;

        currLen = (int)indicesStart.size();

        if (currLen > prevLen) {
            int j = i;
            while (j > 2) {
                std::vector<int> indicesLow = (*minima)[i];
                findOriginMin(minimaPts, minimaSigma, indicesLow, indicesStart, i);
                j = j - 10;
            }
        }
    }
}

void EdgeAnalyzer::getZeroMax(std::vector<int>* zeroPts,
		std::vector<float>* zeroSigma,
		std::vector< std::vector<int> >* zero,
		float threshold)
{
    int len = zero->size();
    float maxSigma = len * sigmaStep;   // biggest sigma

    if (maxSigma <= 3) {
        return;
    }

    float minSigma = roundf(maxSigma * threshold);
    if (minSigma < 1)
        minSigma = 1;

    for (int i = len-2; i > 0; i--) {
        if (i*sigmaStep < minSigma || maxSigma < 2)
            break;
        else {
            if (i == len-2) {
                findSimilar(zeroPts, zeroSigma, (*zero)[i], floorf((*zero)[i].size()/2), i);
            }

            int diff = (int)((*zero)[i].size() - (*zero)[i+1].size());
            if (diff > 1) {
                findSimilar(zeroPts, zeroSigma, (*zero)[i], floor(diff/2), i);
            }
        }
    }
}

Edge* EdgeAnalyzer::resizeEdge(Edge* e, int maxEdgePoints)
{
	//TODO better move to EdgeProcessor, because its not about analyzing
	int eLength = e->length();
	Edge* resizedEdge = new Edge();

	float dEdgeX[eLength];
	float dEdgeY[eLength];
	float arcLength[eLength];
	float sampledArcLength[maxEdgePoints];

	int X[eLength];
	int Y[eLength];

	e->getEdgeXCoordinates(X);
	e->getEdgeXCoordinates(Y);

	float kernel[3] = {-1, 0, 1};
	convolution(X, dEdgeX, int(eLength), kernel, 3);
	convolution(Y, dEdgeY, int(eLength), kernel, 3);


	// calculate arc length of the edge
	arcLength[0] = 0;
	for (int i = 0; i < eLength-1; i++) {
		arcLength[i+1] = arcLength[i] + sqrt(pow(dEdgeX[i], 2) + pow(dEdgeY[i], 2));
	}

	// approximate new edge
	sampledArcLength[0] = 0;
	for (int i = 0; i < maxEdgePoints-1; i++) {
		sampledArcLength[i+1] = (i+1) * arcLength[eLength-1] / (maxEdgePoints-1);
	}

	// match closest value of arcLength to sampledArcLength
	int j = 2;
	std::vector<cv::Point> ePoints = e->getEdgePoints();
	resizedEdge->addPoint(ePoints[0]);// first index

	for (int i = 1; i < maxEdgePoints-1; i++) {
		while (j < eLength-1) {
			if (float curr = arcLength[j]-sampledArcLength[i] >= 0) {
				if (fabs(curr) < fabs(arcLength[j-1]-sampledArcLength[i])) {
					resizedEdge->addPoint(ePoints[j]);
					j++;
					break;
				} else {
					resizedEdge->addPoint(ePoints[j-1]);
					j++;
					break;
				}
			}
			j++;
		}
	}
	resizedEdge->addPoint(ePoints[eLength-1]); // last index

	return resizedEdge;
}

void EdgeAnalyzer::findOriginMax(std::vector<int>* maximaPts, std::vector<float>* maximaSigma, std::vector<int> low, std::vector<int> start, int index){
    int i = 0;
    while (i < start.size()) {
        int ind = findClosest(low, start[i]);
        if (ind != -1 && start[i]-15 < low[ind] && start[i]+15 > low[ind]) {
            if (checkIfExists(*maximaPts, low[ind])) {
                maximaPts->push_back(low[ind]);
                maximaSigma->push_back(index * sigmaStep);
            }
        }
        i++;
    }
}

void EdgeAnalyzer::findOriginMin(std::vector<int>* minimaPts, std::vector<float>* minimaSigma, std::vector<int> low, std::vector<int> start, int index) {
    int i = 0;
    while (i < start.size()) {
        int ind = findClosest(low, start[i]);
        if (ind != -1 && start[i]-15 < low[ind] && start[i]+15 > low[ind]) {
            if (checkIfExists(*minimaPts, low[ind])) {
                minimaPts->push_back(low[ind]);
                minimaSigma->push_back(index * sigmaStep);
            }
        }
        i++;
    }
}

void EdgeAnalyzer::findSimilar(std::vector<int>* zeroPts, std::vector<float>* zeroSigma, std::vector<int> indices, int diff, int index) {
    std::vector<int> diffs;
    std::vector<int> pos;
    for (int i = 0; i < indices.size()-1; i++) {
        if (indices[i+1] - indices[i] < 10) {
            diffs.push_back(indices[i+1] - indices[i]);
            pos.push_back(round((indices[i+1] + indices[i]) / 2));
        }
    }
    for (int i = 0; i < diff; i++) {
        int ind = indexOfSmallestElement(diffs);
        if (diffs.size() > 0) {
            diffs[ind] = 999;
            if (checkIfExists(*zeroPts, pos[ind])) {
                zeroPts->push_back(pos[ind]);
                zeroSigma->push_back(index * sigmaStep);
            }
        }
    }
}

/*
 * (c) Moritz
 */
inline void EdgeAnalyzer::calcKappa(float * xCoorU, float * yCoorU, float * xCoorUU, float * yCoorUU, int length, float * output)
{
    for (int i = 0; i < length; i++) {
        output[i] = (xCoorU[i] * yCoorUU[i] - xCoorUU[i] * yCoorU[i]) / pow(pow(xCoorU[i], 2) + pow(yCoorU[i], 2), 3/2);
    }
}

void EdgeAnalyzer::convolution(int * input, float * output, int iLength, float * kernel, int kLength)
{
    float in[(2*(kLength-1)+iLength)];

    int i = 0;
    while (i < kLength-1) { // pad beginning of signal
        in[i] = float(input[0]);
        ++i;
    }
    i = iLength+kLength-1;
    while (i < 2*kLength+iLength) { // pad end of signal
        in[i] = float(input[iLength-1]);
        ++i;
    }
    i = 0;
    while (i < iLength) { // insert signal
        in[i+kLength-1] = float(input[i]);
        ++i;
    }

    for (int i = 0; i < iLength; i++) {
        output[i] = 0;
        for (int k=0; k<kLength; k++){
            output[i] = output[i] + in[i+kLength/2+k] * kernel[k];
        }
    }

	//std::cout << "size padded sig = " << (2*(kLength-1)+iLength) << std::endl;
}

void EdgeAnalyzer::gauss(int x, float sigma, float * output, int normalized) {
    float sum = 0;

    for (int i = 0; i <= 2*x; i++) {
        output[i] = exp(-pow((i-x) / sigma, 2) / 2) / SQRT_2_PI / sigma;
        sum += fabsf(output[i]);
    }

    if (normalized) {
        for (int i = 0; i <= 2*x; i++) {
            output[i] /= sum;
        }
    }
}

void EdgeAnalyzer::gauss1st(int x, float sigma, float * output, int normalized) {
    float sum = 0;

    for (int i = 0; i <= 2*x; i++) {
        output[i] = exp(-pow((i-x) / sigma, 2) / 2) * (i-x) / SQRT_2_PI / pow(sigma, 3);
        sum += fabsf(output[i]);
    }

    if (normalized) {
        for (int i = 0; i <= 2*x; i++) {
            output[i] /= sum;
        }
    }
}

void EdgeAnalyzer::gauss2nd(int x, float sigma, float * output, int normalized) {
    float sum = 0;

    for (int i = 0; i <= 2*x; i++) {
        output[i] = exp(-pow((i-x) / sigma, 2) / 2) * (pow((i-x), 2) - pow(sigma, 2)) / SQRT_2_PI / pow(sigma, 5);
        sum += fabsf(output[i]);
    }

    if (normalized) {
        for (int i = 0; i <= 2*x; i++) {
            output[i] /= sum;
        }
    }
}

float EdgeAnalyzer::min(float * input, int length) {
    float minimum = 0;
    for (int i = 0; i < length; i++) {
        if (input[i] <= minimum)
            minimum = input[i];
    }

    return minimum;
}

float EdgeAnalyzer::max(float * input, int length) {
    float maximum = 0;
    for (int i = 0; i < length; i++) {
        if (input[i] >= maximum)
            maximum = input[i];
    }

    return maximum;
}



int EdgeAnalyzer::checkIfExists(std::vector<int> arr, int x) {
	//TODO pointer here for arr
    for (int i = 0; i < arr.size(); i++) {
        if (abs(x - arr[i]) < 10)
            return 0;
    }
    return 1;
}

int EdgeAnalyzer::findClosest(std::vector<int> low, int x) {
    std::vector<int> arr;
    for (int i = 0; i < low.size(); i++)
        arr.push_back(abs(low[i] - x));

    return indexOfSmallestElement(arr);
}

int EdgeAnalyzer::indexOfSmallestElement(std::vector<int> array)
{
    int index = -1;
    if (array.size() == 0)
        return index;
    int currMin = array[0];

    for (int i = 0; i < array.size(); i++)
    {
        if (array[i] <= currMin)
            index = i;
            currMin = array[i];
    }

    return index;
}
