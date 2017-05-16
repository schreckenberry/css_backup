#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <vector>
#include <iostream>
#include <stdio.h>

// Util
#include "../util/Edge.h"
#include "../util/InterestPoint.h"
#include "../util/CssPoint.h"
// OpenCV
#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc/imgproc.hpp>

class Visualization
{
private:
	std::vector<cv::Scalar> rgbValues;				//!< vector to store RGB values
	std::vector<cv::Scalar> colorWheel;				//!< vector to store RGB values along colorwheel				
	cv::RNG rng;									//!< OpenCV random number generator

	int canvasWidth = 800;							//!< Width of SVG canvas
	int canvasHeight = 800;							//!< Height of SVG canvas
	/*
	 * fills the vector rgbValues
	 * @param: numberOfValues of RGB Values to be generated
	 */
	void generateRgbValues(int numberOfValues);

public:
	/*
	 * constructor
	 */
	Visualization();

	/*
	 * save analysis result as a bitmap in PNG format
	 */
	void savePNG();
	
	/*
	 * save result of CSS analysis in SVG format
	 * @param edge: vector which contains the traced edges
	 * @param tracedMinSignatures: vector of signatures of local minima in the CSS
	 * @param tracedMaxSignatures: vector of signatures of local maxima in the CSS
	 * @param rotMinSignatures: vector with char. rotations of minima
	 * @param rotMaxSignatures: vector with char. rotations of maxima
	 * @param sF: scale factor for size of radius of circles around extrema
	 * 
	 * TODO: modify so that it works with more than one edge
	 */
	void saveSvgCssResult(Edge* e,
			std::vector< std::vector<CssPoint> > &tracedMinSignatures,
			std::vector< std::vector<CssPoint> > &tracedMaxSignatures,
			std::vector<float> &rotMinSignatures,
			std::vector<float> &rotMaxSignatures,
			float sF);

	/*
	 * save analysis result of ambiguous pixels in SVG format
	 * @param image: binary edge input image (greyscale -> 1 channel)
	 * @param v: vector with edges to be drawn
	 * @param scaleFactor: drawing size of single points in px, default is 4px
	 */
	void saveSVG(cv::Mat &image, std::vector<Edge*>* v, int scaleFactor = 4);
	
	/*
	 * save SVG from OpenCV data, where fragPoints are also in img
	 * (if fragPts is empty only normal image is drawn)
	 */
	void saveCvDataSVG(cv::Mat &img, std::vector<cv::Point> &fragPts);

	/*
	 * save analysis result colored along color wheel in SVG format
	 * @param e: vector with edge to be drawn
	 * @param displayValue: Value to be displayed (e.g. sigma)
	 * @param scaleFactor: drawing size of single points in px, default is 4px
	 */
	void saveSVGcolorWheel(Edge* e, float displayValue = -1, int scaleFactor = 4);

	/*
	 * fills the vector rgbValues with RGB values
	 * @param numberOfValues: number of RGB Values to be generated
	 * @param rgbValues: vector to store color values
	 */
	void generateRgbValues(int numberOfValues, std::vector<cv::Scalar*>* rgbValues);	

	/*
	 * fills the vector colorWheel with RGB values along
	 * the closed color wheel
	 * @param length: number of needed values (must be > 1)
	 */
	void generateColorWheel(int length);

	/*
	 * get value from vector colorWheel
	 * @param i: index of RGB value
	 */
	cv::Scalar getColorWheelValueRGB(int i);

	/*
	 * get value from vector colorWheel
	 * @param i: index of BGR value
	 */
	cv::Scalar getColorWheelValueBGR(int i);

	/*
	 * help function for checking correctness of edge vectors, can be used before
	 * drawing or writing a file and during development
	 * @param image: input image
	 * @param vE: vector with edges
	 * @param vP: vector with detected interest points
	 * @returns: true if no problems detected
	 */
	int validateEdges(cv::Mat& image, std::vector<Edge*>* vE, std::vector<InterestPoint*>* vP);
	
	/*
	 * Set the size of the canvas for drawing SVG images
	 */
	void setCanvasSize(int width, int height);

	/*
	 * write Matlab file which plots kappa over arc length
	 * @param input: array with kappa values
	 * @param length: length of array
	 * @param sigma: sigma value for which kappa was calculated (for display in legend)
	 * @param uMinima: vector which stores the detected local minima
	 * @param uMaxima: vector which stores the detected local maxima
	 * @param thres: threshold which was calculated for selection of local minima and maxima
	 */
	void saveMatlabKappa(float* kappa, int length, float sigma, std::vector<int> uMinima, std::vector<int> uMaxima, float thres);
	void saveMatlabKappaVec(std::vector<float>& kappa, int length, float sigma, std::vector<int> uMinima, std::vector<int> uMaxima, float thres);
	/*
	 * write Matlab file which displays the CSS
	 * @param minima: vector with vectors of local minima
	 * @param maxima: vector with vectors of local maxima
	 * @param sigmaVec: vector with corresponding sigma values
	 */
	void saveMatlabCSS(std::vector< std::vector<int> > minima, std::vector< std::vector<int> > maxima, std::vector<float> sigmaVec);

	/*
	 * save all kappa values in CSV file
	 * @param input: array with kappa values
	 * @param length: length of array
	 * @param sigma: sigma value for which kappa was calculated (for display in legend)
	 */
	void saveCSVKappa(std::vector<float>& kappa, int length, float sigma);
	void saveCSVKappaAppr(std::vector<float>& kappa, int length, float sigma);

	bool firstRunSaveCSVKappaConv;
	bool firstRunSaveCSVKappaAppr;
	
	/*
	 * save signal and corresponding integral contour
	 */
	void saveCSVConvolution(std::vector<float>& sig, int derivative);
	bool firstRunSaveCSVConv0;
	bool firstRunSaveCSVConv1;
	bool firstRunSaveCSVConv2;
	
	void saveCSVApproximation(std::vector<float>& sig, int derivative);
	bool firstRunSaveCSVAppr0;
	bool firstRunSaveCSVAppr1;
	bool firstRunSaveCSVAppr2;
	
	
	/*
	 * save all min and max for every sigma in CSV file
	 * @param input: array with kappa values
	 * @param length: length of array
	 * @param sigma: sigma value for which kappa and corresponding min and max were calculated
	 */
	void saveCSVMinMax(std::vector<int>* uMinima, std::vector<int>* uMaxima, float sigma);
	
	/*
	 * save traced min and max in CSV file together with corresponding color values 
	 * @param tracedMaxSignatures: traced maxima signatures
	 * @param tracedMaxSignatures: traced minima signatures
	 */
	void saveCSVTracedMinMax(std::vector< std::vector<CssPoint> > &tracedMinSignatures, std::vector< std::vector<CssPoint> > &tracedMaxSignatures);
	
	bool firstRunSaveCSVMinMax;							//!< Control flag for function saveCSVMinMax
	/*
	 * method destructor
	 */
	~Visualization();
};

#endif // VISUALIZATION_H
