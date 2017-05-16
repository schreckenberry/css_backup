#include "Edge.h"
#include <iostream>

cv::Moments Edge::getMoments(int index, float radius){
	cv::Moments mom;

	//get vector of points of fragment
	std::vector<cv::Point> fragmentPoints;

	for(int i = 0; i < this->length(); i++)
	{
		if(getDistance(edgePoints[index], edgePoints[i]) <= radius)
		{
			fragmentPoints.push_back(edgePoints[i]);
		}
	}

	//get moments
	mom = cv::moments(fragmentPoints, false );

	return mom;
}

double* Edge::getHuMoments(int index, float radius){
	static double hu[7];

	// 7 Hu-Moments
	cv::HuMoments(getMoments(index, radius), hu);

	return hu;
}

void Edge::setStartInterestPoint(InterestPoint* p)
{
	startInterestPoint = p;
}

void Edge::setEndInterestPoint(InterestPoint* p)
{
	endInterestPoint = p;
}

InterestPoint* Edge::getStartInterestPoint()
{
	return startInterestPoint;
}

InterestPoint* Edge::getEndInterestPoint()
{
	return endInterestPoint;
}


std::vector<cv::Point> Edge::getEdgePoints()
{
	return edgePoints;
}

void Edge::getEdgeXCoordinates(int* a){
	for(int i = 0; i < this->length(); i++){
		a[i] = edgePoints[i].x;
	}
}

void Edge::getEdgeYCoordinates(int* a){
	for(int i = 0; i < this->length(); i++){
		a[i] = edgePoints[i].y;
	}
}

std::vector<int> Edge::getEdgeXCoordinatesVector(){
	std::vector<int> v;
	for(int i = 0; i < this->length(); i++){
			v.push_back(edgePoints[i].x);
		}
	return v;
}
std::vector<int> Edge::getEdgeYCoordinatesVector(){
	std::vector<int> v;
		for(int i = 0; i < this->length(); i++)
		{
			v.push_back(edgePoints[i].y);
		}
		return v;
}

void Edge::invertRecursive(){
	// invert previous edges
	for (std::vector<Edge*>::iterator it = followingEdges.begin();
				it != followingEdges.end(); ++it)
		{
			(*it)->invertRecursive();
		}

	// invert the points
	std::vector<cv::Point> newEdgePoints;
	for (std::vector<cv::Point>::reverse_iterator it = edgePoints.rbegin();
				it != edgePoints.rend(); ++it)
		{
			newEdgePoints.push_back((*it));
		}
	edgePoints = newEdgePoints;

	// switch previous and following edges
	std::vector<Edge*> tmp = previousEdges;
	previousEdges = followingEdges;
	followingEdges = tmp;

	return;
}

void Edge::addPoint(cv::Point p)
{
	edgePoints.push_back(p);
}

int Edge::length()
{
	return edgePoints.size();
}

void Edge::printPoints()
{
	for (std::vector<cv::Point>::iterator it = edgePoints.begin();
			it != edgePoints.end(); ++it)
	{
		std::cout << "[" << it->x << ", " << it->y << "]" << std::endl;
	}
}

cv::Point Edge::getStartPoint()
{
	return edgePoints.front();
}

cv::Point Edge::getEndPoint()
{
	return edgePoints.back();
}

void Edge::addFollowingEdge(Edge* e)
{
	followingEdges.push_back(e);
	return;
}

void Edge::addFollowingEdge(std::vector<Edge*> v)
{
	for (std::vector<Edge*>::iterator it = v.begin(); it != v.end(); it++)
	{
		followingEdges.push_back((*it));
	}
	//why no reference? because getFollowingEdges() returns a vector of objects
	//TODO later change to reference for perfonmance increase
	return;
}

void Edge::addPreviousEdge(Edge* e)
{
	previousEdges.push_back(e);
	return;
}

void Edge::addPreviousEdge(std::vector<Edge*> v)
{
	for (std::vector<Edge*>::iterator it = v.begin(); it != v.end(); it++)
	{
		previousEdges.push_back((*it));
	}
	//why no reference? because getPreviousEdges() returns a vector of objects
	//TODO later change to reference for perfonmance increase
	return;
}

int Edge::getNumberFollowingEdges()
{
	return followingEdges.size();
}

int Edge::getNumberPreviousEdges()
{
	return previousEdges.size();
}

std::vector<Edge*> Edge::getFollowingEdges()
{
	// maybe return reference here? no, because someone could manipulate the object
	return followingEdges;
}

std::vector<Edge*> Edge::getPreviousEdges()
{
	// maybe return reference here? no, because someone could manipulate the object
	return previousEdges;
}

void Edge::deleteFollowingEdge(Edge* e)
{
	std::vector<Edge*>::iterator it = std::find(followingEdges.begin(), followingEdges.end(), e);
	if(it != followingEdges.end()){ //if element is not found return value is .end()
		followingEdges.erase(it);
	}

	return;
}

void Edge::deletePreviousEdge(Edge* e)
{
	std::vector<Edge*>::iterator it = std::find(previousEdges.begin(), previousEdges.end(), e);
	if(it != previousEdges.end()){ //if element is not found return value is .end()
		previousEdges.erase(it);
	}

	return;
}

void Edge::deleteStartPoint()
{
	edgePoints.erase(edgePoints.begin());
	return;
}

void Edge::deleteEndPoint()
{
	edgePoints.pop_back();
	return;
}

void Edge::appendEdge(Edge* e)
{
	std::vector<cv::Point> newPoints = e->getEdgePoints();

	// do not appent the fist point, because it overlaps
	for(int i = 0; i < newPoints.size(); i++){
		this->edgePoints.push_back(newPoints[i]);
	}

	//get following edges and append
	followingEdges = e->getFollowingEdges();

	return;
}

void Edge::appendInverseEdge(Edge* e)
{
	// TODO what happens to following/previous?
	// this method is only used for loop closing -> maybe rename it to make the name more clear
	std::vector<cv::Point> newPoints = e->getEdgePoints();
	for (std::vector<cv::Point>::reverse_iterator it = newPoints.rbegin();
			it != newPoints.rend(); ++it)
	{
		this->edgePoints.push_back(*it);
	}
}

float Edge::getStartDirection(int length)
{
	if (this->length() >= length)
	{
		float dx = this->edgePoints.at(length - 1).x - this->edgePoints.at(0).x;
		float dy = this->edgePoints.at(length - 1).y - this->edgePoints.at(0).y;
		return atan2(dy, dx);
	}
	else
	{
		float dx = this->edgePoints.at(this->edgePoints.size() - 1).x
				- this->edgePoints.at(0).x;
		float dy = this->edgePoints.at(this->edgePoints.size() - 1).y
				- this->edgePoints.at(0).y;
		return atan2(dy, dx);
	}
}

float Edge::getEndDirection(int length)
{
	if (this->length() >= length)
	{
		float dx = this->edgePoints.at(this->edgePoints.size() - 1).x
				- this->edgePoints.at(this->edgePoints.size() - length).x;
		float dy = this->edgePoints.at(this->edgePoints.size() - 1).y
				- this->edgePoints.at(this->edgePoints.size() - length).y;
		return atan2(dy, dx);
	}
	else
	{
		float dx = this->edgePoints.at(this->edgePoints.size() - 1).x
				- this->edgePoints.at(0).x;
		float dy = this->edgePoints.at(this->edgePoints.size() - 1).y
				- this->edgePoints.at(0).y;
		return atan2(dy, dx);
	}
}

void Edge::setIsLoop(bool b)
{
	isLoop = b;
	return;
}

bool Edge::getIsLoop()
{
	return isLoop;
}

float Edge::getDistance(cv::Point p1, cv::Point p2)
{
	cv::Point diff = p1 - p2;
	return cv::sqrt(diff.x*diff.x + diff.y*diff.y);
}

Edge::Edge()
{
	isLoop = false;
}

Edge::Edge(std::vector<cv::Point> edgePointVector)
{
	edgePoints = edgePointVector;
	isLoop = false;
}

Edge::Edge(cv::Point p)
{
	edgePoints.push_back(p);
	isLoop = false;
}

Edge::Edge(cv::Point p, InterestPoint* iP)
{
	startInterestPoint = iP;
	edgePoints.push_back(p);
	isLoop = false;
}

Edge::Edge(Edge* e)
{
	std::vector<cv::Point> v = e->getEdgePoints();
	for (std::vector<cv::Point>::iterator it = v.begin(); it != v.end(); ++it)
	{
		this->addPoint(*it);
	}
	isLoop = false;
}

Edge::~Edge()
{
	//delete all elements
	edgePoints.clear();
	followingEdges.clear();
	previousEdges.clear();
}
