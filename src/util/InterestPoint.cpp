#include "InterestPoint.h"

cv::Point InterestPoint::getPoint()
{
	return position;
}

void InterestPoint::setPosition(cv::Point p)
{
	position = p;
	return;
}

int InterestPoint::getX()
{
	return position.x;
}

void InterestPoint::setX(int x)
{
	position.x = x;
	return;
}

int InterestPoint::getY()
{
	return position.y;
}

void InterestPoint::setY(int y)
{
	position.y = y;
	return;
}

void InterestPoint::addStartingEdge(Edge* e)
{
	startingEdges.push_back(e);
}

void InterestPoint::addEndingEdge(Edge* e)
{
	endingEdges.push_back(e);
}

void InterestPoint::switchStartingAndEndingEdges()
{
	std::vector<Edge*> tmp = startingEdges;
	startingEdges = endingEdges;
	endingEdges = tmp;
}

InterestPoint::InterestPoint(cv::Point p)
{
	position = p;
}

InterestPoint::~InterestPoint()
{

}
