#ifndef EDGE_H
#define EDGE_H

// OpenCV
#include <opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#include "InterestPoint.h"

class Edge
{
private:
	std::vector<cv::Point> edgePoints;
	std::vector<Edge*> followingEdges;
	std::vector<Edge*> previousEdges;
	InterestPoint* startInterestPoint;
	InterestPoint* endInterestPoint;
	bool isLoop;

public:
	/*
	 * method constructor
	 */
	Edge();
	Edge(std::vector<cv::Point> edgePointVector);
	Edge(cv::Point p);
	Edge(cv::Point p, InterestPoint* iP);
	/*
	 * makes a deep copy of the referenced edge
	 * @param e: referenced edge
	 */
	Edge(Edge *e);
	/*
	 * method destructor
	 */
	~Edge();
	void setStartInterestPoint(InterestPoint* p);
	void setEndInterestPoint(InterestPoint* p);
	InterestPoint* getStartInterestPoint();
	InterestPoint* getEndInterestPoint();
	/*
	 * method returns all cv::Point on Edge in std::vector
	 */
	std::vector<cv::Point> getEdgePoints();
	/*
	 * method returns all x coordinates of all point on Edge in an array
	 */
	void getEdgeXCoordinates(int* a);
	/*
	 * method returns all y coordinates of all point on Edge in an array
	 */
	void getEdgeYCoordinates(int* a);
	/*
	 * method returns all x coordinates of all point on Edge in an array
	 */
	std::vector<int> getEdgeXCoordinatesVector();
	/*
	 * method returns all y coordinates of all point on Edge in an array
	 */
	std::vector<int> getEdgeYCoordinatesVector();
	/*
	 * method adds a cv::Point to the edgePoints std::vector
	 * @param p: cv::Point to be added
	 */
	void addPoint(cv::Point p);
	/*
	 * method returns length of edgePoints std::vector
	 */
	int length();
	/*
	 * method prints all cv::Point in edgePoints std::vector to console/terminal
	 */
	void printPoints();
	/*
	 * method returns first cv::Point in Edge
	 */
	cv::Point getStartPoint();
	/*
	 * method returns last cv::Point in Edge
	 */
	cv::Point getEndPoint();
	/*
	 * method adds a following edge to the vector
	 */
	void addFollowingEdge(Edge* e);
	void addFollowingEdge(std::vector<Edge*> v);
	/*
	 * method adds a previous edge to the vector
	 */
	void addPreviousEdge(Edge* e);
	void addPreviousEdge(std::vector<Edge*> v);
	/*
	 * method deleted Edge with reference e out of the vector followingEdges
	 * @param e: reference to edge that will be deleted
	 */
	void deleteFollowingEdge(Edge* e);
	/*
	 * method deleted Edge with reference e out of the vector previousEdges
	 * @param e: reference to edge that will be deleted
	 */
	void deletePreviousEdge(Edge* e);
	/*
	 * method deletes first point in EdgePoints
	 */
	void deleteStartPoint();
	/*
	 * method deletes last point in EdgePoints
	 */
	void deleteEndPoint();
	/*
	 * method appends given edge to this edge
	 */
	void appendEdge(Edge* e);
	/*
	 * method appends given edge from back to front to this edge
	 */
	void appendInverseEdge(Edge* e);
	/*
	 * method return number of following edges
	 */
	int getNumberFollowingEdges();
	/*
	 * method return number of previous edges
	 */
	int getNumberPreviousEdges();
	/*
	 * method return vector of following edges
	 */
	std::vector<Edge*> getFollowingEdges();
	/*
	 * method return vector of following edges
	 */
	std::vector<Edge*> getPreviousEdges();
	/*
	 * method returns angle of edge start with length in rad
	 */
	float getStartDirection(int length);
	/*
	 * method returns angle of edge end with length in rad
	 */
	float getEndDirection(int length);
	/*
	 * method sets private variable isLoop
	 * @param b: value isLoop is set to
	 */
	void setIsLoop(bool b);
	/*
	 * method returns private variable isLoop
	 */
	bool getIsLoop();
	/*
	 * method inverts the edge and all previous edges
	 */
	void invertRecursive();
	/*
	 * returns a vector of moments
	 * @param index: position of fragment
	 * @param radius: radius of fragment
	 */
	cv::Moments getMoments(int index, float radius);
	double* getHuMoments(int index, float radius);

	float getDistance(cv::Point p1, cv::Point p2);
};

#endif // EDGE_H
