#include "Visualization.h"

bool drawOnlyFragmentArea = 0;
#define COLOR_SCHEME 1 // 0 = black, 1 = white

#if (COLOR_SCHEME == 0)
char bgColor[] = {"rgb(0,0,0)"};
char edgePxColor[] = {"rgb(255,255,255)"};
char triMaxColor[] = {"rgb(0,128,0)"};
char triMinColor[] = {"rgb(128,0,0)"};
char circFill[] = {"rgb(45,45,45)"};
char circStroke[] = {"rgb(38,38,38)"};
#endif
#if (COLOR_SCHEME == 1)
char bgColor[] = {"rgb(255,255,255)"};
char edgePxColor[] = {"rgb(0,0,0)"};
char triMaxColor[] = {"rgb(0,200,0)"};
char triMinColor[] = {"rgb(200,0,0)"};
char circFill[] = {"rgb(220,220,220)"};
char circStroke[] = {"rgb(200,200,200)"};
#endif

// p stores coordinates, indices stores indices of edges sharing the point
struct VisInterestPoint
{
	cv::Point p;
	std::vector<int> indices;
};

Visualization::Visualization()
{
	rng(12345); // OpenCV random number generator state init
	firstRunSaveCSVKappaConv = true;
	firstRunSaveCSVKappaAppr = true;
	firstRunSaveCSVMinMax = true;
	firstRunSaveCSVConv0 = true;
	firstRunSaveCSVConv1 = true;
	firstRunSaveCSVConv2 = true;
	firstRunSaveCSVAppr0 = true;
	firstRunSaveCSVAppr1 = true;
	firstRunSaveCSVAppr2 = true;
}

void Visualization::savePNG()
{
	//TODO: implement this function
}

void Visualization::generateRgbValues(int numberOfValues)
{
	rgbValues.clear();
	
	cv::Scalar color;
	for (int i = 0; i < numberOfValues; i++)
	{
		color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		rgbValues.push_back(color);
		//std::cout << rgbValues.at(i)->val[0] << std::endl;
	}
	std::cout << "Generated " << rgbValues.size() << " random RGB values" << std::endl;
}

void Visualization::generateColorWheel(int length)
{
	// Intervals: H[0,179], S[0,255], V[0,255]
	
	colorWheel.clear();

	// 1x1 px images for conversion
	cv::Mat colorRGB(1, 1, CV_8UC3, cv::Scalar(0, 0, 0));
	cv::Mat colorHSV(1, 1, CV_8UC3, cv::Scalar(0, 255, 255));

	/*
	std::cout << "R = " << (int) colorRGB.at<cv::Vec3b>(0,0)[0] << std::endl;
	std::cout << "G = " << (int) colorRGB.at<cv::Vec3b>(0,0)[1] << std::endl;
	std::cout << "B = " << (int) colorRGB.at<cv::Vec3b>(0,0)[2] << std::endl;
	*/
	
	int hueValue;
	cv::Scalar color;
	for (int i = 0; i < length; i++)
	{
		hueValue = 179.0 * i / (length-1);
		
		// write interval-scaled HSV value
		colorHSV.at<cv::Vec3b>(0,0)[0] = hueValue;
		colorHSV.at<cv::Vec3b>(0,0)[1] = 255;
		colorHSV.at<cv::Vec3b>(0,0)[2] = 255;
		
		// calculate corresponding RGB value
		cv::cvtColor(colorHSV, colorRGB, CV_HSV2RGB);
		
		color = cv::Scalar(colorRGB.at<cv::Vec3b>(0,0)[0],
							colorRGB.at<cv::Vec3b>(0,0)[1],
							colorRGB.at<cv::Vec3b>(0,0)[2]);
		
		// save value
		colorWheel.push_back(color);
	}
	
	std::cout << "Generated " << colorWheel.size() << " ordered color values" << std::endl;
}

cv::Scalar Visualization::getColorWheelValueRGB(int i)
{
	return colorWheel.at(i);
}

cv::Scalar Visualization::getColorWheelValueBGR(int i)
{
	return cv::Scalar(colorWheel[i].val[2], colorWheel[i].val[1], colorWheel[i].val[0]);
}

void Visualization::generateRgbValues(int numberOfValues, std::vector<cv::Scalar*>* rgbValues)
{
	//TODO: implement this function for external generation of RGB values
}

int Visualization::validateEdges(cv::Mat& image, std::vector<Edge*>* vE, std::vector<InterestPoint*>* vP)
{
	cv::Mat sharedPixels(image.size(), CV_8UC1, 0.0);
	std::vector<cv::Point> valIPVector;

	// go through all edges and count number of edges at a coordinate
	for (int i = 0; i < vE->size(); i++)
	{
		// go trough points of edge i
		for (int j = 0; j < vE->at(i)->getEdgePoints().size(); j++)
		{
			sharedPixels.at<uchar>(vE->at(i)->getEdgePoints()[j])++;
		}
	}
	
	// now collect the coordinates with more than one edge (called interest points)
	//VisInterestPoint visInterestPointTmp;
	int numberOfEdgePointsSharedPixels = 0; // which are not an interest point
	int numberOfInterestPoints = 0; // which are not an interest point
	int numberOfContributions = 0;
	
	for (int y = 0; y < sharedPixels.rows; y++)
	{
		for (int x = 0; x < sharedPixels.cols; x++)
		{
			if (sharedPixels.at<uchar>(y, x) > 1.0)
			{
				//visInterestPointTmp.p.x = x;
				//visInterestPointTmp.p.y = y;
				valIPVector.push_back(cv::Point(x, y));
				
				numberOfInterestPoints++;
				numberOfContributions = numberOfContributions + sharedPixels.at<uchar>(y, x);
			}
			else if (sharedPixels.at<uchar>(y, x) == 1.0)
			{
				numberOfEdgePointsSharedPixels++;
			}
		}
	}
	
	// count pixels in edges which are not an interest point
	int numberOfEdgePointsInEdges = 0;
	bool isIP = false;

	for (int i = 0; i < vE->size(); i++)
	{
		for (int j = 0; j < vE->at(i)->getEdgePoints().size(); j++)
		{
			// check if point is one of the interest points
			isIP = false;
			for (int k = 0; k < valIPVector.size(); k++)
			{
				// note that case can only be entered once because only one edge point is
				// analyzed and compared with all interest points
				if ((vE->at(i)->getEdgePoints()[j].x == valIPVector.at(k).x) &&
					(vE->at(i)->getEdgePoints()[j].y == valIPVector.at(k).y))
				{
					isIP = true;
				}
			}
			
			// it is an edge point if it is not an interest point
			if (!isIP)
			{
				numberOfEdgePointsInEdges++;
			}
		}
	}
	
	/*
	// print interest points vectors in matlab compatible format
	std::cout << "vP = [";
	for (int i = 0; i < vP->size(); i++)
	{
		std::cout << vP->at(i)->getX() << " " << vP->at(i)->getY() << "; ";
	}
	std::cout << "]" << std::endl;

	std::cout << "valIPVector = [";
	for (int i = 0; i < valIPVector.size(); i++)
	{
		std::cout << valIPVector.at(i).x << " " << valIPVector.at(i).y << "; ";
	}
	std::cout << "]" << std::endl;
	*/
	
	// check if there is a difference between external and interal interest points
	// take i-th interest point from internal vector ...
	// TODO: implement check vice versa (if internal vector is larger then external)
	for (int i = 0; i < valIPVector.size(); i++)
	{
		// ... and go trough all interest points of external vector valIPVector.size(); vP->size()
		int pointFound = 0;
		int jTmp = 0;

		for (int j = 0; j < vP->size(); j++)
		{
			if (valIPVector.at(i) == vP->at(j)->getPoint())
			{
				pointFound++;
			}
		}
		
		if (pointFound == 0)
		{
			std::cout << "Inconsistent Interest Point at (x=" << valIPVector.at(i).x
				<< ",y=" << valIPVector.at(i).y << ")" << std::endl;
		}
		else if (pointFound > 1)
		{
			std::cout << "Inconsistent Interest Point at (x=" << valIPVector.at(i).x
				<< ",y=" << valIPVector.at(i).y << ") - Found" << pointFound << " times" << std::endl;
		}
	}
	
	// calculate length of all edges together
	int accumulatedLength = 0;
	for (int i = 0; i < vE->size(); i++)
	{
		accumulatedLength = accumulatedLength + vE->at(i)->getEdgePoints().size();
	}
	
	std::cout << "numberOfInterestPoints = " << numberOfInterestPoints << std::endl;
	std::cout << "accumulatedLength = " << accumulatedLength << std::endl;
	std::cout << "numberOfContributions = " << numberOfContributions << std::endl;
	std::cout << "accumulatedLength - numberOfContributions = "
		<< accumulatedLength - numberOfContributions << std::endl;
	std::cout << "numberOfEdgePointsSharedPixels = " << numberOfEdgePointsSharedPixels << std::endl;
	std::cout << "numberOfEdgePointsInEdges = " << numberOfEdgePointsInEdges << std::endl;

	//TODO: Add further comparisons of the values above to increase security
	if (numberOfEdgePointsSharedPixels == numberOfEdgePointsInEdges)
	{
		return 1;
	}
	else
	{
		std::cout << "Error: Incorrect Edge Points (Visualization::validateEdges)" << std::endl;
		return 0;
	}
}

void Visualization::saveSVGcolorWheel(Edge* e, float displayValue, int scaleFactor)
{
	int eLength = e->length();
	std::vector<cv::Point> p = e->getEdgePoints();

	FILE *file;
	file = fopen("./gauss.svg", "w");

	//fprintf (file, "<svg>\n");
	// SVG canvas with black background
	fprintf (file, "<svg width=\"%d\" height=\"%d\">",
		canvasWidth * scaleFactor, canvasHeight * scaleFactor);
	fprintf (file, "<rect width=\"%d\" height=\"%d\" style=\"fill:rgb(0,0,0);\" />",
		canvasWidth * scaleFactor, canvasHeight * scaleFactor);
	
	generateColorWheel(eLength);
	
	for (int i = 0; i < eLength; i++)
	{
		// format: <rect x="0" y="0" width="400" height="180" style="fill:rgb(0,0,0);" />
			fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(%d,%d,%d);\" />",
				(double) scaleFactor * p[i].x,
				(double) scaleFactor * p[i].y,
				(double) scaleFactor, (double) scaleFactor,
				(int) colorWheel[i].val[0], (int) colorWheel[i].val[1], (int) colorWheel[i].val[2]);
	}
	
	// mark first element
	fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(255,255,255);\" />",
		(double) scaleFactor * p[0].x,
		(double) scaleFactor * p[0].y,
		(double) scaleFactor, (double) scaleFactor);
	
	// display value at top right if passed
	if (displayValue != -1)
	{
		int xPos = scaleFactor * canvasWidth - scaleFactor * 150;
		int yPos = 25 * scaleFactor;
		int factor = 20 * scaleFactor;
		fprintf (file, "<text x=\"%d\" y=\"%d\" style=\"fill:white;font-size: %dpx;\">&#963; = %f</text>",
					xPos, yPos, factor, displayValue);
	}
	
	fprintf (file, "</svg>");
	fclose (file);
	std::cout << "Writing file gauss.svg with " << eLength << " points finished." << std::endl;
}

void Visualization::saveSvgCssResult(Edge* e,
		std::vector< std::vector<CssPoint> > &tracedMinSignatures,
		std::vector< std::vector<CssPoint> > &tracedMaxSignatures,
		std::vector<float> &rotMinSignatures,
		std::vector<float> &rotMaxSignatures,
		float sF)
{
	// drawing determines layer, therefore single functions for every object

	/*** parameters ******************************************************/
	float os = 3.0; // offset of edge points of triangles (=> size control)
	/*********************************************************************/
	
	// variables for reuse
	int uCoord = 0;
	float xCoord = 0.0; float yCoord = 0.0;
	
	// new SVG file
	FILE *file;
	file = fopen("./cssresult.svg", "w");
	
	int eLength = e->length();
	std::vector<cv::Point> p = e->getEdgePoints();
	
	// SVG canvas with black or white background
	fprintf (file, "<svg width=\"%d\" height=\"%d\">\n", canvasWidth, canvasHeight);
	fprintf (file, "<rect width=\"%d\" height=\"%d\" style=\"fill:%s;\" />\n", canvasWidth, canvasHeight, bgColor);

	
	// draw circles for maxima
	// format: <circle cx="40" cy="40" r="24" />
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// determine origin
		uCoord = tracedMaxSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;
		
		float cScale = tracedMaxSignatures[i].back().sigma; // characteristic scale
		
		fprintf (file, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"fill:%s;\" fill-opacity=\"1.0\" />\n",
			xCoord+0.5, yCoord+0.5 , cScale*sF, circFill);
	}

	// draw circles for minima
	// format: <circle cx="40" cy="40" r="24" />
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// determine origin
		uCoord = tracedMinSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;
		
		float cScale = tracedMinSignatures[i].back().sigma; // characteristic scale

		fprintf (file, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"fill:%s;\" fill-opacity=\"1.0\" />\n",
			xCoord+0.5, yCoord+0.5, cScale*sF, circFill);
	}
	
	//draw rotations of maxima
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// determine origin
		uCoord = tracedMaxSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;
		
		float cScale = tracedMaxSignatures[i].back().sigma; // characteristic scale

		// format: <line x1="0" y1="0" x2="200" y2="200" style="stroke:rgb(255,0,0);stroke-width:2" />
		fprintf (file, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:%s;stroke-width:0.5\" />\n",
			xCoord+0.5, yCoord+0.5, xCoord+0.5+(cScale*sF*cos(rotMaxSignatures[i])), yCoord+0.5+(cScale*sF*sin(rotMaxSignatures[i])), triMaxColor);
	}

	// draw rotations of minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// determine origin
		uCoord = tracedMinSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;
		
		float cScale = tracedMinSignatures[i].back().sigma; // characteristic scale

		// format: <line x1="0" y1="0" x2="200" y2="200" style="stroke:rgb(255,0,0);stroke-width:2" />
		fprintf (file, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:%s;stroke-width:0.5\" />\n",
			xCoord+0.5, yCoord+0.5, xCoord+0.5+(cScale*sF*cos(rotMinSignatures[i])), yCoord+0.5+(cScale*sF*sin(rotMinSignatures[i])), triMinColor);
	}
	
	// draw circle strokes of maxima
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// determine origin
		uCoord = tracedMaxSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;
		
		float cScale = tracedMaxSignatures[i].back().sigma; // characteristic scale

		// stroke-dasharray: 2 1
		fprintf (file, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke:%s; stroke-width:0.5; fill:none;\" />\n",
			xCoord+0.5, yCoord+0.5 , cScale*sF, circStroke);
	}

	// draw circle strokes of minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// determine origin
		uCoord = tracedMinSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;
		
		float cScale = tracedMinSignatures[i].back().sigma; // characteristic scale

		// stroke-dasharray: 2 1
		fprintf (file, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke:%s; stroke-width:0.5; fill:none;\" />\n",
			xCoord+0.5, yCoord+0.5 , cScale*sF, circStroke);
	}

	// draw triangles for maxima
	// format: <path d="M100 100 L150 200 L50 200 L100 100" .../> 1st point (tip), 2nd point, 3rd point, last point
	// subpixel also start at top left, therefore correction terms
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// take first element (position on lowest scale)
		uCoord = tracedMaxSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;

		fprintf (file, "<path d=\"M%f %f L%f %f L%f %f L%f %f\" style=\"fill:%s;\" fill-opacity=\"1.0\" />\n",
			xCoord+0.5, yCoord-os-0.5, xCoord+os+0.5, yCoord+os-0.5, xCoord-os+0.5, yCoord+os-0.5, xCoord+0.5, yCoord-os-0.5, triMaxColor);
	}
	
	// draw triangles for minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// take first element (position on lowest scale)
		uCoord = tracedMinSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;

		fprintf (file, "<path d=\"M%f %f L%f %f L%f %f L%f %f\" style=\"fill:%s;\" fill-opacity=\"1.0\" />\n",
			xCoord+0.5, yCoord+os+1.5, xCoord-os+0.5, yCoord-os+1.5, xCoord+os+0.5, yCoord-os+1.5, xCoord+0.5, yCoord+os+1.5, triMinColor);
	}
	
	// draw pixels of edge
	for (int i = 0; i < eLength; i++)
	{
		// format: <rect x="0" y="0" width="400" height="180" style="fill:rgb(0,0,0);" />
		fprintf (file, "<rect x=\"%d\" y=\"%d\" width=\"1\" height=\"1\" style=\"fill:%s;\" />\n",
			p[i].x, p[i].y, edgePxColor);
	}
	
	// draw center pixels of maxima
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// take first element (position on lowest scale)
		uCoord = tracedMaxSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;

		fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"1\" height=\"1\" style=\"fill:rgb(0,255,0);\" />\n", xCoord, yCoord);
	}

	// draw center pixels of minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// take first element (position on lowest scale)
		uCoord = tracedMinSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;

		fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"1\" height=\"1\" style=\"fill:rgb(255,0,0);\" />\n", xCoord, yCoord);
	}
	
	// draw numbers of maxima
	for (int i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// take first element (position on lowest scale)
		uCoord = tracedMaxSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;

		fprintf (file, "<text x=\"%f\" y=\"%f\" style=\"fill:rgb(0,255,0); font-size:6px; font-family:arial;\">%d</text>\n",
			xCoord+5, yCoord+1.5, i);
	}
	
	// draw numbers of minima
	for (int i = 0; i < tracedMinSignatures.size(); i++)
	{
		// take first element (position on lowest scale)
		uCoord = tracedMinSignatures[i].at(0).u;
		xCoord = p[uCoord].x; yCoord = p[uCoord].y;

		fprintf (file, "<text x=\"%f\" y=\"%f\" style=\"fill:rgb(255,0,0); font-size:6px; font-family:arial;\">%d</text>\n",
			xCoord+5, yCoord+2.5, i);
	}
	
	fprintf (file, "</svg>");
	fclose (file);
	std::cout << "Writing file cssresult.svg finished." << std::endl;
}

void Visualization::setCanvasSize(int width, int height)
{
	canvasWidth = width;
	canvasHeight = height;
}

void Visualization::saveSVG(cv::Mat& image, std::vector<Edge*>* v, int scaleFactor)
{
	generateRgbValues(v->size());
	std::vector<VisInterestPoint> visInterestPointVector; // see global def. of struct interestPoint

	// go through all edges and count number of edges at a coordinate
	cv::Mat sharedPixels(image.size(), CV_8UC1, 0.0);
	for (int i = 0; i < v->size(); i++)
	{
		// go trough points of edge i
		for (int j = 0; j < v->at(i)->getEdgePoints().size(); j++)
		{
			sharedPixels.at<uchar>(v->at(i)->getEdgePoints()[j])++;
		}
	}
	
	// now collect the coordinates with more than one edge (called interest points)
	VisInterestPoint visInterestPointTmp;
	for (int y = 0; y < sharedPixels.rows; y++)
	{
		for (int x = 0; x < sharedPixels.cols; x++)
		{
			if (sharedPixels.at<uchar>(y, x) > 1.0)
			{
				visInterestPointTmp.p.x = x;
				visInterestPointTmp.p.y = y;
				visInterestPointVector.push_back(visInterestPointTmp);
			}
		}
	}
	
	//std::cout << "visInterestPointVector.size() = " << visInterestPointVector.size() << std::endl;
	
	// again go trough all edges and save indices of edges sharing a coordinate
	for (int i = 0; i < v->size(); i++)
	{
		// go trough points j of edge i
		for (int j = 0; j < v->at(i)->getEdgePoints().size(); j++)
		{
			// go trough the vector with the k shared coordinates
			for (int k = 0; k < visInterestPointVector.size(); k++)
			{
				// check if point j of edge i shares coordinate and save index
				if (visInterestPointVector.at(k).p == v->at(i)->getEdgePoints()[j])
				{
					visInterestPointVector.at(k).indices.push_back(i);
				}
			}
		}
	}

	/*********** Write SVG file ******************************************/
	FILE *file;
	file = fopen("./edges.svg", "w");
	
	// SVG canvas with black background
	fprintf (file, "<svg width=\"%d\" height=\"%d\">",
		image.cols * scaleFactor, image.rows * scaleFactor);
	fprintf (file, "<rect width=\"%d\" height=\"%d\" style=\"fill:rgb(0,0,0);\" />",
		image.cols * scaleFactor, image.rows * scaleFactor);

	// draw all edges except for interest points
	int isIP = false;
	for (int i = 0; i < v->size(); i++)
	{
		for (int j = 0; j < v->at(i)->getEdgePoints().size(); j++)
		{
			// j should be counted up without drawing the pixel if in vector visInterestPointVector
			isIP = false;

			// check if point is one of the interest points
			for (int k = 0; k < visInterestPointVector.size(); k++)
			{
				// note that this case can only be entered once because only one edge point is
				// analyzed and compared with all interest points
				if ((v->at(i)->getEdgePoints()[j].x == visInterestPointVector.at(k).p.x) &&
					(v->at(i)->getEdgePoints()[j].y == visInterestPointVector.at(k).p.y))
				{
					isIP = true;
				}
			}
			
			if (isIP)
			{
				continue; // next j (next point of edge i)
			}
			
			// format: <rect x="0" y="0" width="400" height="180" style="fill:rgb(0,0,0);" />
			fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(%d,%d,%d);\" />",
				(double) scaleFactor * v->at(i)->getEdgePoints()[j].x,
				(double) scaleFactor * v->at(i)->getEdgePoints()[j].y,
				(double) scaleFactor, (double) scaleFactor,
				(int) rgbValues.at(i).val[0], (int) rgbValues.at(i).val[1], (int) rgbValues.at(i).val[2]);
				//128, 25, 135);
				
			// draw number of edge
			fprintf (file, "<text x=\"%f\" y=\"%f\" style=\"fill:rgb(255,255,255); font-size:0.5px; font-family:arial;\">%d</text>\n",
				(double) scaleFactor * v->at(i)->getEdgePoints()[j].x + scaleFactor * 0.25,
				(double) scaleFactor * v->at(i)->getEdgePoints()[j].y + scaleFactor * 0.75,
				(int) i);
		}
	}

	// draw interest points
	for (int k = 0; k < visInterestPointVector.size(); k++)
	{

		if (visInterestPointVector.at(k).indices.size() == 1 ||
			visInterestPointVector.at(k).indices.size() > 8)
		{
			std::cout << "Error: Incorrect Interest Point (Visualization::saveSVG)." << std::endl;
		}
		else
		{
			for (int l = visInterestPointVector.at(k).indices.size(); l > 0; l--)
			{
				// decrease height for overlay
				double hFactor = ((double) l) / visInterestPointVector.at(k).indices.size();

				fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(%d,%d,%d);\" />",
					(double) scaleFactor * visInterestPointVector.at(k).p.x,
					(double) scaleFactor * visInterestPointVector.at(k).p.y,
					(double) scaleFactor, (double) hFactor * scaleFactor,
					(int) rgbValues.at(visInterestPointVector.at(k).indices.at(l-1)).val[0],
					(int) rgbValues.at(visInterestPointVector.at(k).indices.at(l-1)).val[1],
					(int) rgbValues.at(visInterestPointVector.at(k).indices.at(l-1)).val[2]);
			}
		}
		
		// draw/overlay all interest points with defined color
		/*
		fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(%d,%d,%d);\" />",
			(double) scaleFactor * visInterestPointVector.at(k).p.x,
			(double) scaleFactor * visInterestPointVector.at(k).p.y,
			(double) scaleFactor, (double) scaleFactor,
			255, 255, 255);
		*/

		// draw/overlay border around interest points with defined color
		fprintf (file, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\""
			" style=\"stroke-width:0.1;stroke:rgb(%d,%d,%d);fill:none;\" />",
			(double) scaleFactor * visInterestPointVector.at(k).p.x,
			(double) scaleFactor * visInterestPointVector.at(k).p.y,
			(double) scaleFactor, (double) scaleFactor,
			255, 255, 255);
	}
	
	fprintf (file, "</svg>");
	fclose (file);
	std::cout << "Writing file edges.svg finished." << std::endl;
	/*********************************************************************/
}

void Visualization::saveCSVKappa(std::vector<float>& kappa, int length, float sigma)
{
	FILE *file;
	// new file if new run
	if (firstRunSaveCSVKappaConv)
	{
		std::cout << "New file kappaConv.csv started." << std::endl;

		file = fopen("./Data/kappaConv.csv", "w");
		firstRunSaveCSVKappaConv = false;
	}
	// otherwise just add new lines
	else
	{
		file = fopen("./Data/kappaConv.csv", "a");
	}
	
	/*********** Write CSV file ******************************************/
	/* Format for every line in CSV file is:
	 * sigma, kappa_1st_px, kappa_2nd_px, ..., kappa_nth_px
	 */
	fprintf (file, "%f, ", sigma);
	// all kappa values
	for (int i = 0; i < length; i++)
	{
		// no comma behind last element, but new line
		if (i == (length-1))
		{
			fprintf (file, "%f\n", kappa[i]);		
		}
		else
		{
			fprintf (file, "%f, ", kappa[i]);
		}
	}
	/*********************************************************************/

	fclose (file);
}

void Visualization::saveCSVKappaAppr(std::vector<float>& kappa, int length, float sigma)
{
	FILE *file;
	// new file if new run
	if (firstRunSaveCSVKappaAppr)
	{
		std::cout << "New file kappaAppr.csv started." << std::endl;

		file = fopen("./Data/kappaAppr.csv", "w");
		firstRunSaveCSVKappaAppr = false;
	}
	// otherwise just add new lines
	else
	{
		file = fopen("./Data/kappaAppr.csv", "a");
	}
	
	/*********** Write CSV file ******************************************/
	/* Format for every line in CSV file is:
	 * sigma, kappa_1st_px, kappa_2nd_px, ..., kappa_nth_px
	 */
	fprintf (file, "%f, ", sigma);
	// all kappa values
	for (int i = 0; i < length; i++)
	{
		// no comma behind last element, but new line
		if (i == (length-1))
		{
			fprintf (file, "%f\n", kappa[i]);		
		}
		else
		{
			fprintf (file, "%f, ", kappa[i]);
		}
	}
	/*********************************************************************/

	fclose (file);
}

void Visualization::saveCSVConvolution(std::vector<float>& sig, int derivative)
{
	FILE *fileIntCon;
	
	if (derivative == 0)
	{
		if (firstRunSaveCSVConv0)
		{
			fileIntCon = fopen("./Data/conv0.csv","w");
			std::cout << "New file conv0.csv started." << std::endl;
			firstRunSaveCSVConv0 = false;
		}
		else
			fileIntCon = fopen("./Data/conv0.csv","a"); //just add new lines
	}
	else if (derivative == 1)
	{
		if (firstRunSaveCSVConv1)
		{
			fileIntCon = fopen("./Data/conv1.csv","w");
			std::cout << "New file conv1.csv started." << std::endl;
			firstRunSaveCSVConv1 = false;
		}
		else
			fileIntCon = fopen("./Data/conv1.csv","a"); //just add new lines
	}
	else if (derivative == 2)
	{
		if (firstRunSaveCSVConv2)
		{
			fileIntCon = fopen("./Data/conv2.csv","w");
			std::cout << "New file conv2.csv started." << std::endl;
			firstRunSaveCSVConv2 = false;
		}
		else
			fileIntCon = fopen("./Data/conv2.csv","a"); //just add new lines
	}
	
		
	/*********** Write CSV file ******************************************/
	/*
	 * One Line per convolution result
	 */
	int sigSize = sig.size();

	for (int i = 0; i < sigSize; i++)
	{		
		// no comma behind last element, but new line
		if (i == (sigSize-1))
		{
			fprintf (fileIntCon, "%f\n", sig[i]);		
		}
		else
		{
			fprintf (fileIntCon, "%f, ", sig[i]);
		}
	}	
	fclose (fileIntCon);
}

void Visualization::saveCSVApproximation(std::vector<float>& sig, int derivative)
{
	FILE *fileIntCon;
	
	if (derivative == 0)
	{
		if (firstRunSaveCSVAppr0)
		{
			fileIntCon = fopen("./Data/appr0.csv","w");
			std::cout << "New file appr0.csv started." << std::endl;
			firstRunSaveCSVAppr0 = false;
		}
		else
			fileIntCon = fopen("./Data/appr0.csv","a"); //just add new lines
	}
	else if (derivative == 1)
	{
		if (firstRunSaveCSVAppr1)
		{
			fileIntCon = fopen("./Data/appr1.csv","w");
			std::cout << "New file appr1.csv started." << std::endl;
			firstRunSaveCSVAppr1 = false;
		}
		else
			fileIntCon = fopen("./Data/appr1.csv","a"); //just add new lines
	}
	else if (derivative == 2)
	{
		if (firstRunSaveCSVAppr2)
		{
			fileIntCon = fopen("./Data/appr2.csv","w");
			std::cout << "New file appr2.csv started." << std::endl;
			firstRunSaveCSVAppr2 = false;
		}
		else
			fileIntCon = fopen("./Data/appr2.csv","a"); //just add new lines
	}
		
	/*********** Write CSV file ******************************************/
	/*
	 * One Line per approximation result
	 */
	int sigSize = sig.size();

	for (int i = 0; i < sigSize; i++)
	{		
		// no comma behind last element, but new line
		if (i == (sigSize-1))
		{
			fprintf (fileIntCon, "%f\n", sig[i]);		
		}
		else
		{
			fprintf (fileIntCon, "%f, ", sig[i]);
		}
	}	
	fclose (fileIntCon);
}

void Visualization::saveCSVMinMax(std::vector<int>* uMinima, std::vector<int>* uMaxima, float sigma)
{
	FILE *fileMin;
	FILE *fileMax;
	// new files if new run
	if (firstRunSaveCSVMinMax)
	{
		fileMin = fopen("./Data/min.csv", "w");
		fileMax = fopen("./Data/max.csv", "w");
		
		std::cout << "New files min.csv and max.csv started." << std::endl;
		firstRunSaveCSVMinMax = false;
	}
	// otherwise just add new lines
	else
	{
		fileMin = fopen("./Data/min.csv", "a");
		fileMax = fopen("./Data/max.csv", "a");
	}

	/*********** Write CSV file ******************************************/
	/* Format for every line in CSV file is:
	 * min.csv: sigma, #min, minU_1, minU_2, ..., minU_n
	 * max.csv: sigma, #max, maxU_1, maxU_2, ..., maxU_n
	 */
	/*********** write min ***********/
	int nuMin = (*uMinima).size();
	if (nuMin > 0)
	{
		// print with comma
		fprintf (fileMin, "%f, ", sigma);
		fprintf (fileMin, "%d, ", nuMin);
	}
	else
	{
		// print without comma and enter new line
		fprintf (fileMin, "%f, ", sigma);
		fprintf (fileMin, "0\n");
	}

	for (int i = 0; i < nuMin; i++)
	{		
		// no comma behind last element, but new line
		if (i == (nuMin-1))
		{
			fprintf (fileMin, "%d\n", (*uMinima)[i]);		
		}
		else
		{
			fprintf (fileMin, "%d, ", (*uMinima)[i]);
		}
	}
	/*********************************/
	
	/*********** write max ***********/
	int nuMax = (*uMaxima).size();
	if (nuMax > 0)
	{
		// print with comma
		fprintf (fileMax, "%f, ", sigma);
		fprintf (fileMax, "%d, ", nuMax);
	}
	else
	{
		// print without comma and enter new line
		fprintf (fileMax, "%f, ", sigma);
		fprintf (fileMax, "0\n");
	}

	for (int i = 0; i < nuMax; i++)
	{		
		// no comma behind last element, but new line
		if (i == (nuMax-1))
		{
			fprintf (fileMax, "%d\n", (*uMaxima)[i]);	
		}
		else
		{
			fprintf (fileMax, "%d, ", (*uMaxima)[i]);
		}
	}
	/*********************************/

	/*********************************************************************/

	fclose (fileMin);
	fclose (fileMax);
}

void Visualization::saveMatlabKappa(float* kappa, int length, float sigma, std::vector<int> uMinima, std::vector<int> uMaxima, float thres)
{
	std::cout << "uMinima.size() = " << uMinima.size() << std::endl;
	std::cout << "uMaxima.size() = " << uMaxima.size() << std::endl;

	/*********** Write Matlab file ***************************************/	
	FILE *file;
	file = fopen("./kappa.m", "w");
	fprintf (file, "figure;\n");
	
	// I - plot graph of points
	// save sigma and threshold as variable
	fprintf (file, "s=%f*ones(1,%d);\n", sigma, length);
	fprintf (file, "T=%f;\n", thres);
	// create u vector with appropriate length
	fprintf (file, "u=0:1:%d;\n", length-1);

	fprintf (file, "k=[");
	for (int i = 0; i < length; i++)
	{
		fprintf (file, "%f ", kappa[i]);
	}
	fprintf (file, "];\n");
	fprintf (file, "plot(u, k, \'Color\', [0.7 0.7 0.7]);\n");
	fprintf (file, "hold on;\n");
	
	// II - plot single points with color
	generateColorWheel(length);
	float r, g, b;
	for (int i = 0; i < length; i++)
	{
		r = colorWheel.at(i).val[0] / 255.0;
		g = colorWheel.at(i).val[1] / 255.0;
		b = colorWheel.at(i).val[2] / 255.0;
		
		fprintf (file, "plot(%d, %f, \'.\', \'Color\', [%f %f %f]);\n", i, kappa[i], r, g, b);
	}

	// legend and labels
	fprintf (file, "legend(\'\\sigma = %f\');\n", sigma);
	fprintf (file, "xlabel(\'arc length (pixel)\');\n");
	fprintf (file, "ylabel(\'kappa\');\n");
	
	// plot values
	for (int i = 0; i < uMaxima.size(); i++)
	{
		fprintf (file, "plot(%d,k(%d),\'g*\');\n", uMaxima[i], uMaxima[i]+1);
	}
	for (int i = 0; i < uMinima.size(); i++)
	{
		fprintf (file, "plot(%d,k(%d),\'r*\');\n", uMinima[i], uMinima[i]+1);
	}
	
	// plot threshold
	fprintf (file, "line([%d %d], [%f %f]);\n", 0, length, thres, thres);
	fprintf (file, "line([%d %d], [%f %f]);\n", 0, length, -thres, -thres);
	/*********************************************************************/
	
	fclose (file);
	std::cout << "Writing file kappa.m finished." << std::endl;
}

void Visualization::saveMatlabKappaVec(std::vector<float>& kappa, int length, float sigma, std::vector<int> uMinima, std::vector<int> uMaxima, float thres)
{
	std::cout << "uMinima.size() = " << uMinima.size() << std::endl;
	std::cout << "uMaxima.size() = " << uMaxima.size() << std::endl;

	/*********** Write Matlab file ***************************************/	
	FILE *file;
	//file = fopen("./kappaAppr.m", "w");
	file=fopen("./../matlab/kappaApprDev.m", "w");
	fprintf (file, "figure;\n");
	
	// I - plot graph of points
	// save sigma and threshold as variable
	fprintf (file, "s=%f*ones(1,%d);\n", sigma, length);
	fprintf (file, "T=%f;\n", thres);
	// create u vector with appropriate length
	fprintf (file, "u=0:1:%d;\n", length-1);

	fprintf (file, "k=[");
	for (int i = 0; i < length; i++)
	{
		fprintf (file, "%f ", kappa[i]);
	}
	fprintf (file, "];\n");
	fprintf (file, "plot(u, k, \'Color\', [0.7 0.7 0.7]);\n");
	fprintf (file, "hold on;\n");
	
	// II - plot single points with color
	generateColorWheel(length);
	float r, g, b;
	for (int i = 0; i < length; i++)
	{
		r = colorWheel.at(i).val[0] / 255.0;
		g = colorWheel.at(i).val[1] / 255.0;
		b = colorWheel.at(i).val[2] / 255.0;
		
		fprintf (file, "plot(%d, %f, \'.\', \'Color\', [%f %f %f]);\n", i, kappa[i], r, g, b);
	}

	// legend and labels
	fprintf (file, "legend(\'\\sigma = %f\');\n", sigma);
	fprintf (file, "xlabel(\'arc length (pixel)\');\n");
	fprintf (file, "ylabel(\'kappa\');\n");
	
	// plot values
	for (int i = 0; i < uMaxima.size(); i++)
	{
		fprintf (file, "plot(%d,k(%d),\'g*\');\n", uMaxima[i], uMaxima[i]+1);
	}
	for (int i = 0; i < uMinima.size(); i++)
	{
		fprintf (file, "plot(%d,k(%d),\'r*\');\n", uMinima[i], uMinima[i]+1);
	}
	
	// plot threshold
	fprintf (file, "line([%d %d], [%f %f]);\n", 0, length, thres, thres);
	fprintf (file, "line([%d %d], [%f %f]);\n", 0, length, -thres, -thres);
	/*********************************************************************/
	
	fclose (file);
	std::cout << "Writing file kappa.m finished." << std::endl;
}

void Visualization::saveMatlabCSS(std::vector< std::vector<int> > minima, std::vector< std::vector<int> > maxima, std::vector<float> sigmaVec)
{
	if (minima.size() != sigmaVec.size())
	{
		std::cout << "Error: Visualization::saveMatlabCSS: Sizes do not match." << std::endl;
		return;
	}
	
	/*********** Write Matlab file ***************************************/
	FILE *file;
	file = fopen("./css.m", "w");

	fprintf(file, "figure;\n");
	fprintf(file, "hold on;\n");

	// go trough all sigma values
	for (int i = 0; i < sigmaVec.size(); i++)
	{
		float sigma = sigmaVec[i];
		
		// minima coordinates (plot every value pair)
		for (int j = 0; j < minima[i].size(); j++)
		{
			fprintf (file, "plot(%d,%f,\'r.\');\n", minima[i].at(j), sigma);
		}

		// mmaxima coordinates (plot every value pair)
		for (int j = 0; j < maxima[i].size(); j++)
		{
			fprintf (file, "plot(%d,%f,\'g.\');\n", maxima[i].at(j), sigma);
		}
	}
	
	// legend and labels
	fprintf (file, "legend(\'\\sigma_{max} = %f\');\n", sigmaVec.back());
	fprintf (file, "xlabel(\'arc length (pixel)\');\n");
	fprintf (file, "ylabel(\'\\sigma (pixel)\');\n");
	
	fclose (file);
	std::cout << "Writing file css.m finished." << std::endl;
	/*********************************************************************/
}

void Visualization::saveCSVTracedMinMax(std::vector< std::vector<CssPoint> > &tracedMinSignatures, std::vector< std::vector<CssPoint> > &tracedMaxSignatures)
{
	FILE *fileTracedMin;
	FILE *fileTracedMax;

	fileTracedMin = fopen("./Data/tracedMin.csv", "w");
	fileTracedMax = fopen("./Data/tracedMax.csv", "w");
	
	std::cout << "New files tracedMin.csv and tracedMax.csv started." << std::endl;

	/*********** Write CSV files *****************************************/
	/* Format for every line in CSV file (every line = one signature):
	 * tracedMin.csv: #min, R, G, B, u1, s1, k1, u2, s2, k2, ...
	 * tracedMax.csv: #max, R, G, B, u1, s1, k1, u2, s2, k2, ...
	 */

	/*********** write min ***********/
	// write one 0 if no minima
	if (tracedMinSignatures.size() == 0) { fprintf (fileTracedMin, "0, 0, 0, 0\n"); };
	
	// one traced signature for every i
	for (size_t i = 0; i < tracedMinSignatures.size(); i++)
	{
		// one j contains one [ui, si, ki]
		for (size_t j = 0; j < tracedMinSignatures[i].size(); j++)
		{
			// save corresponding color value and size
			if (j == 0)
			{
				// size
				fprintf (fileTracedMin, "%d, ", (int) tracedMinSignatures[i].size());
				
				// color value
				int u1 = tracedMinSignatures[i].at(0).u;
				fprintf (fileTracedMin, "%d, %d, %d, ",
					(int) colorWheel[u1].val[0],
					(int) colorWheel[u1].val[1],
					(int) colorWheel[u1].val[2]);
			}
			
			// write ui, si, ki
			fprintf (fileTracedMin, "%d, ", tracedMinSignatures[i].at(j).u);
			fprintf (fileTracedMin, "%f, ", tracedMinSignatures[i].at(j).sigma);
			// no comma behind last value
			if ((j+1) == tracedMinSignatures[i].size())
			{
				fprintf (fileTracedMin, "%f", tracedMinSignatures[i].at(j).kappa);
			}
			else
			{
				fprintf (fileTracedMin, "%f, ", tracedMinSignatures[i].at(j).kappa);
			}
		}

		fprintf (fileTracedMin, "\n");
	}
	/*********************************/

	/*********** write max ***********/
	// write one 0 if no minima
	if (tracedMaxSignatures.size() == 0) { fprintf (fileTracedMax, "0, 0, 0, 0\n"); };

	// one traced signature for every i
	for (size_t i = 0; i < tracedMaxSignatures.size(); i++)
	{
		// one j contains one [ui, si, ki]
		for (size_t j = 0; j < tracedMaxSignatures[i].size(); j++)
		{
			// save corresponding color value and size
			if (j == 0)
			{
				// size
				fprintf (fileTracedMax, "%d, ", (int) tracedMaxSignatures[i].size());
				
				// color value
				int u1 = tracedMaxSignatures[i].at(0).u;
				fprintf (fileTracedMax, "%d, %d, %d, ",
					(int) colorWheel[u1].val[0],
					(int) colorWheel[u1].val[1],
					(int) colorWheel[u1].val[2]);
			}

			// write ui, si, ki
			fprintf (fileTracedMax, "%d, ", tracedMaxSignatures[i].at(j).u);
			fprintf (fileTracedMax, "%f, ", tracedMaxSignatures[i].at(j).sigma);
			// no comma behind last value
			if ((j+1) == tracedMaxSignatures[i].size())
			{
				fprintf (fileTracedMax, "%f", tracedMaxSignatures[i].at(j).kappa);
			}
			else
			{
				fprintf (fileTracedMax, "%f, ", tracedMaxSignatures[i].at(j).kappa);
			}
		}

		fprintf (fileTracedMax, "\n");
	}
	/*********************************/
	
	fclose(fileTracedMin);
	fclose(fileTracedMax);
	/*********************************************************************/
}

void Visualization::saveCvDataSVG(cv::Mat &img, std::vector<cv::Point> &fragPts)
{
	// new SVG file
	FILE *file;
	file = fopen("./fragment.svg", "w");
	
	// SVG canvas with black background
	fprintf (file, "<svg width=\"%d\" height=\"%d\">\n", img.cols, img.rows);
	fprintf (file, "<rect width=\"%d\" height=\"%d\" style=\"fill:rgb(45,45,45);\" />\n", img.cols, img.rows);
	
	/*
	std::cout << "img.cols = " << img.cols << std::endl;
	std::cout << "img.rows = " << img.rows << std::endl;
	*/
	
	// draw area
	for (size_t x = 0; x < img.cols; x++)
	{
		for (size_t y = 0; y < img.rows; y++)
		{
			if (img.at<uchar>(y, x) > 0)
			{
				// format: <rect x="0" y="0" width="400" height="180" style="fill:rgb(0,0,0);" />
				fprintf (file, "<rect x=\"%d\" y=\"%d\" width=\"1\" height=\"1\" style=\"fill:rgb(255,255,255);\" />\n",
					(int) x, (int) y);
			}
		}
	}
	
	if (!drawOnlyFragmentArea)
	{
		// draw fragment
		for (size_t i = 0; i < fragPts.size(); i++)
		{
			// format: <rect x="0" y="0" width="400" height="180" style="fill:rgb(0,0,0);" />
			fprintf (file, "<rect x=\"%d\" y=\"%d\" width=\"1\" height=\"1\" style=\"fill:rgb(255,0,0);\" />\n",
				fragPts[i].x, fragPts[i].y);
		}

		// numbering of fragment points
		for (size_t i = 0; i < fragPts.size(); i++)
		{		
			fprintf (file, "<text x=\"%f\" y=\"%f\" style=\"fill:rgb(0,0,255); font-size:0.5px; font-family:arial;\">%d</text>\n",
				fragPts[i].x+0.25, fragPts[i].y+0.75, (int) i);
		}
	}
	
	fprintf (file, "</svg>");
	fclose (file);
	std::cout << "Writing file fragment.svg finished." << std::endl;
}

Visualization::~Visualization()
{
	// delete all elements
	rgbValues.clear();
	colorWheel.clear();
}
