#ifndef INTERESTPOINT_H
#define INTERESTPOINT_H

// OpenCV
#include <opencv2/core/core.hpp>

class Edge; //resolve circular dependency of header files

class InterestPoint
{
	private:
		cv::Point position;
		std::vector<Edge*> startingEdges;
		std::vector<Edge*> endingEdges;

	public:
		/*
		 * method constructor
		 */
		InterestPoint(cv::Point p);
		/*
		 * method destructor
		 */
		~InterestPoint();

		cv::Point getPoint();
		void setPosition(cv::Point p);
		int getX();
		void setX(int x);
		int getY();
		void setY(int y);
		void addStartingEdge(Edge* e);
		void addEndingEdge(Edge* e);
		void switchStartingAndEndingEdges();

};

#endif // INTERESTPOINT_H
