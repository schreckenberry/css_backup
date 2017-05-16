#ifndef EDGE_PROCESSOR_H
#define EDGE_PROCESSOR_H

#include <vector>
// Util
#include "../util/Edge.h"
#include "../util/InterestPoint.h"
// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

class EdgeProcessor
{
public:
	/*
	 * method finds all edges in sorted order with OpenCV funtion findContours()
	 * source: http://docs.opencv.org/doc/tutorials/imgproc/shapedescriptors/find_contours/find_contours.html
	 * @param image: binary edge input image
	 */
	std::vector<Edge> getContours(cv::Mat& image);
	/*
	 * method finds all edges in sorted manner recursively by method of Daniel Gaspers
	 * @param image: binary edge input image (greyscale -> 1 channel)
	 * @param edgeVector: vector traced Edges get pushed to
	 */
	void getEdges(cv::Mat& image, std::vector<Edge*>* edgeVector);
	/*
	 * method finds all edges in sorted manner recursively by method of Daniel Gaspers
	 * @param image: binary edge input image (greyscale -> 1 channel)
	 * @param edgeVector: vector traced Edges get pushed to
	 * @param pointVecto: vector traced InterstPoints get pushed to
	 */
	void getEdgesAndInterestPoints(cv::Mat& image,
			std::vector<Edge*>* edgeVector,
			std::vector<InterestPoint*>* pointVector);
	/*
	 * method delete Edges shorter then T / append Edges based on following Edges
	 * @param edgeVector: vector containing references to Edges to be postprocessed
	 * @param T: threshold for deleting short Edges
	 * @param lengthForDirection: number of pixels that will be included in calculation the start/end direction of an Edge
	 * @param angleDifference: max. difference between direction of 2 edges for appending
	 */
	void postProcessEdges(std::vector<Edge*>* edgeVector, int T,
			int lengthForDirection, float angleDifference);
			
	/*
	 * method closes loops
	 * @param edgeVector: vector containing references to Edges to be postprocessed
	 * TODO: Make private again
	 */
	static void closeLoops(std::vector<Edge*>* edgeVector);

private:
	/*
	 * method returns distance between the given points
	 * @param p1: first point
	 * @param o2: second point
	 */
	static float getDistance(cv::Point p1, cv::Point p2);
	/*
	 * method recursive edge tracing by Daniel Gaspers
	 * @param image: binary edge input image (greyscale -> 1 channel)
	 * @param visited: image to store visited edges (1: visited 0: not visited)
	 * @param startPoint: Point to start each recursion step from
	 * @param e: Edge which will be extended
	 * @param vE: the vector traced Edges get pushed to when recursion ends
	 * @param vP: the vector traced InterestPoints get pushed to when recursion ends
	 */
	static void traceEdgeAndInterestPoints(cv::Mat& img, cv::Mat& visited,
				cv::Point startPoint, Edge* e, std::vector<Edge*>* vE,
				std::vector<InterestPoint*>* vP);
	/*
	 * method generates 8 neighbours for the given point clockwise, starting top left
	 * @param p: point who's neighbours will be generated
	 */
	static std::vector<cv::Point> get8Neighbours(cv::Point p);
	/*
	 * method generates 8 neighbours for the given point, including previous point in edge
	 * @param p1: point n in edge, who's neighbours will be generated ()
	 * @param p2: point n-1 in edge
	 */
	static std::vector<cv::Point> getNeighbours(cv::Point lastPoint,
			cv::Point previousPoint);
	/*
	 * method delete adges that are shorter then T
	 * @param edgeVector: vector containing references to Edges to be postprocessed
	 * @param T: threshold
	 */
	static void deleteShortEdges(std::vector<Edge*>* edgeVector, int T);
	/*
	 * method appends Edges
	 * @param edgeVector: vector containing references to Edges to be postprocessed
	 * @param lengthForDirection: length that is used for calculation of direction
	 * @param angleDifference: difference that to edges can have
	 */
	static void appendEdges(std::vector<Edge*>* edgeVector,
			int lengthForDirection, float angleDifference);
			
	static void selectionOfNeighbours(cv::Point startPoint, std::vector<cv::Point> &neighbours, cv::Mat &img, cv::Mat &visited);
};

#endif // EDGE_PROCESSOR_H
