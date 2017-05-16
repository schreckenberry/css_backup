#include "EdgeProcessor.h"
#include <math.h>
#include <iostream>
//for GUI stuff can be deleted if no imshow() is used
//#include <opencv2/highgui/highgui.hpp>

using namespace cv;

void EdgeProcessor::postProcessEdges(std::vector<Edge*>* edgeVector, int T,
		int lengthForDirection, float angleDifference)
{
	// just loop closing of two edges
	std::cout << ">>>LoopClosing started" << std::endl;
	closeLoops(edgeVector);
	std::cout << ">>>LoopClosing finished" << std::endl;

	//TODO also append loops?
	std::cout << ">>>EdgeAppending started" << std::endl;
	appendEdges(edgeVector, lengthForDirection, angleDifference);
	std::cout << ">>>EdgeAppending finished" << std::endl;

	std::cout << ">>>EdgeDeleting started" << std::endl;
	deleteShortEdges(edgeVector, T);
	std::cout << ">>>EdgeDeleting finished" << std::endl;
	return;
}

void EdgeProcessor::closeLoops(std::vector<Edge*>* edgeVector)
{
	bool loopClosed = false;
	//iterate trough all edges
	for (std::vector<Edge*>::iterator baseEdge = edgeVector->begin();
			baseEdge != edgeVector->end(); baseEdge++)
	{
		loopClosed = false;

		// check if following edges exist
		if ((*baseEdge)->getNumberFollowingEdges() >= 2) //need of min 2 edges to close a loop
		{
			// if it is the case -> get vector of following edges
			std::vector<Edge*> followingEdges = (*baseEdge)->getFollowingEdges();

			//iterate through following edges of baseEdge
			for (std::vector<Edge*>::iterator followingEdge1 =
					followingEdges.begin();
					followingEdge1 != followingEdges.end(); followingEdge1++)
			{
				//iterate through following edges of baseEdge a second time
				for (std::vector<Edge*>::iterator followingEdge2 = followingEdges.begin();
						followingEdge2 != followingEdges.end();
						followingEdge2++)
				{
					//check if endpoints are neighboured -> distance <= sqrt(2)
					if (getDistance((*followingEdge1)->getEndPoint(), (*followingEdge2)->getEndPoint()) <= sqrt(2)
							&& (*followingEdge1) != (*followingEdge2))
					{// delete one point if overlapping
						if((*followingEdge1)->getStartPoint() == (*followingEdge2)->getStartPoint()){
							(*followingEdge2)->deleteStartPoint();
						}

						// append the reversed edge
						(*followingEdge1)->appendInverseEdge(*followingEdge2);
						// mark as loop
						(*followingEdge1)->setIsLoop(true);
						loopClosed = true;
						std::cout << "   closed a loop! resulting edge: " << (*followingEdge1)->getStartPoint() << (*followingEdge1)->getEndPoint() << std::endl;

						// delete appended Edge from vector of all Followers
						(*baseEdge)->deleteFollowingEdge((*followingEdge2));

						// delete appended Edge from vector of all Edges
						std::vector<Edge*>::iterator it = std::find(edgeVector->begin(), edgeVector->end(), (*followingEdge2));
						if(it != edgeVector->end()){ //if element is not found return value is .end()
							edgeVector->erase(it);
						}

						//delete edge from all previous edges
						std::vector<Edge*> previousEdges = (*followingEdge2)->getPreviousEdges();
						for (std::vector<Edge*>::iterator previousEdge =previousEdges.begin();
								previousEdge != previousEdges.end();
								previousEdge++)
						{
							(*previousEdge)->deleteFollowingEdge(*followingEdge2);
						}

						//delete(*followingEdge2);

						//break followingEdge2
						break;
					}
				}

				if(loopClosed)
				{
					//break followingEdge1
					break;
				}

			}
		}

		if(loopClosed)
		{
			//when loop closed start iterator from beginning,
			//because we manipulated data in front of the iterators position
			//-> this invalidates the iterator and can result in segfaults
			baseEdge = edgeVector->begin();
		}

	}
}

void EdgeProcessor::deleteShortEdges(std::vector<Edge*>* edgeVector, int T)
{
	for (std::vector<Edge*>::iterator it = edgeVector->begin();
			it != edgeVector->end(); it++)
	{
		if ((*it)->length() < T)
		{
			// delete short ends
			if ((*it)->getNumberFollowingEdges() == 0)
			{
				std::cout << "   deleted a short end edge: " << (*it)->getStartPoint() << (*it)->getEndPoint() << std::endl;
				// delete Edge from vector of all Edges
				edgeVector->erase(it);

				//delete edge from all previous edges
				std::vector<Edge*> previousEdges = (*it)->getPreviousEdges();
				for (std::vector<Edge*>::iterator previousEdge = previousEdges.begin();
						previousEdge != previousEdges.end();
						previousEdge++)
				{
					(*previousEdge)->deleteFollowingEdge(*it);
				}
				//TODO
				//delete(*it);

				it = edgeVector->begin();
			}
			// delete short starts
			else if((*it)->getNumberPreviousEdges() == 0)
			{
				std::cout << "   deleted a short start edge: " << (*it)->getStartPoint() << (*it)->getEndPoint() << std::endl;
				// delete Edge from vector of all Edges
				edgeVector->erase(it);

				//delete edge from all following edges
				std::vector<Edge*> followingEdges = (*it)->getFollowingEdges();
				for (std::vector<Edge*>::iterator followingEdge = followingEdges.begin();
						followingEdge != followingEdges.end();
						followingEdge++)
				{
					(*followingEdge)->deletePreviousEdge(*it);
				}
				//TODO
				//delete(*it);

				it = edgeVector->begin();
			}
		}
	}

	return;
}

void EdgeProcessor::appendEdges(std::vector<Edge*>* edgeVector, int lengthForDirection, float angleDifference)
{
	bool edgeAppended = false;

	for (std::vector<Edge*>::iterator baseEdge = edgeVector->begin();
			baseEdge != edgeVector->end(); baseEdge++)
	{
		edgeAppended = false;
		// check if following edges exist
		if ((*baseEdge)->getNumberFollowingEdges() > 0)
		{
			// get vector of following edges
			std::vector<Edge*> followingEdges = (*baseEdge)->getFollowingEdges();

			for (std::vector<Edge*>::iterator followingEdge =
					followingEdges.begin();
					followingEdge != followingEdges.end(); followingEdge++)
			{
				// check the difference in angle and the start and endpoint between the edges
				if (fabs( (*baseEdge)->getEndDirection(lengthForDirection)
							- (*followingEdge)->getStartDirection(lengthForDirection) ) <= angleDifference
							&& (*baseEdge)->getEndPoint() == (*followingEdge)->getStartPoint())
				{

					// delete overlapping point
					if((*followingEdge)->getStartPoint() == (*baseEdge)->getEndPoint()){
						(*followingEdge)->deleteStartPoint();
					}

					// append the edge
					(*baseEdge)->appendEdge(*followingEdge);

					// delete appended Edge from vector of all Edges
					std::vector<Edge*>::iterator it = std::find(edgeVector->begin(), edgeVector->end(), (*followingEdge));
					if(it != edgeVector->end()){ //if element is not found return value is .end()
						edgeVector->erase(it);
					}

					std::cout << "   appended an edge! resulting edge: " << (*followingEdge)->getStartPoint() << (*followingEdge)->getEndPoint() << std::endl;
					edgeAppended = true;

					// delete edge as follower of baseEdge
					(*baseEdge)->deleteFollowingEdge(*followingEdge);

					// set new follower of baseEdge
					(*baseEdge)->addFollowingEdge( (*followingEdge)->getFollowingEdges() );

					// delete and set previous edges of followers follower
					std::vector<Edge*> followersFollowingEdges = (*followingEdge)->getFollowingEdges();
					for (std::vector<Edge*>::iterator followersFollowingEdge = followersFollowingEdges.begin();
							followersFollowingEdge != followersFollowingEdges.end(); followersFollowingEdge++)
					{
						(*followersFollowingEdge)->deletePreviousEdge(*followingEdge);
						(*followersFollowingEdge)->addPreviousEdge(*baseEdge);
					}
					//TODO
					//delete(*followingEdge);

					// end the iteration over followingEdges, because match was found
					break;
				}
			}
		}
		if(edgeAppended)
		{
			baseEdge = edgeVector->begin();
		}
	}
}

std::vector<Edge> EdgeProcessor::getContours(Mat& image)
{
	std::vector<Edge> edgeVector;

	std::vector<Point> contours;
	std::vector<Vec4i> hierarchy;

	// Find contours
	findContours(image, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE,
			Point(0, 0));

	// generate edgeVector
	for (int i = 0; i < contours.size(); i++)
	{
		edgeVector.push_back(Edge(contours.at(i)));
	}

	return edgeVector;
}

void EdgeProcessor::getEdges(Mat& image, std::vector<Edge*>* edgeVector)
{
	cv::Mat visitedPoints(image.size(), image.type(), 0.0);
	std::vector<InterestPoint*>* pointVector;

	for (int y = 1; y < image.rows - 1; y++) // start at 1 because neighbours must exist for tracing
	{
		for (int x = 1; x < image.cols - 1; x++) // start at 1 because neighbours must exist for tracing
		{
			if (image.at<uchar>(y, x) == 255
					&& visitedPoints.at<uchar>(y, x) == 0) // singe grayscale image [0,255]
			{
				traceEdgeAndInterestPoints(image, visitedPoints, Point(x, y),
										new Edge(), edgeVector, pointVector);
			}
		}
	}

	return;
}

void EdgeProcessor::getEdgesAndInterestPoints(Mat& image,
		std::vector<Edge*>* edgeVector,
		std::vector<InterestPoint*>* pointVector)
{
	cv::Mat visitedPoints(image.size(), image.type(), 0.0);

	for (int y = 1; y < image.rows - 1; y++) // start at 1 because neighbours must exist for tracing
	{
		for (int x = 1; x < image.cols - 1; x++) // start at 1 because neighbours must exist for tracing
		{
			if (image.at<uchar>(y, x) == 255
					&& visitedPoints.at<uchar>(y, x) == 0) // singe grayscale image [0,255]
			{
				traceEdgeAndInterestPoints(image, visitedPoints, Point(x, y),
						new Edge(), edgeVector, pointVector);
			}
		}
	}

	return;
}

void EdgeProcessor::traceEdgeAndInterestPoints(cv::Mat& img, cv::Mat& visited,
		cv::Point startPoint, Edge* e, std::vector<Edge*>* vE,
		std::vector<InterestPoint*>* vP)
{
	//TODO maybe use pointer for startPoint to improve performance
	e->addPoint(startPoint);
	visited.at<uchar>(startPoint) = 1; // mark as visited
	std::vector<cv::Point> neighbours = get8Neighbours(startPoint);
	selectionOfNeighbours(startPoint, neighbours, img, visited);

	if (e->length() != 1)
	/*
	 * common case2 -> n-1 point in edge exists
	 * therefore follow continuity
	 */
	{
		// compare if first and last in array are direct neighbours
		if (neighbours.size() > 1)
		{
			if (getDistance(neighbours[0],neighbours[neighbours.size()-1]) == 1)
			{
				visited.at<uchar>(neighbours[0].y, neighbours[0].x) = 0;
				neighbours.erase(neighbours.begin());
			}
		}

		int possibilities = neighbours.size();

		// one pixel with more than one possibility
		if (possibilities > 1)
		/*
		 * multiple possibilities
		 */
		{
			//end edge we came from here
			vE->push_back(e);
			InterestPoint* iP = new InterestPoint(startPoint);
			vP->push_back(iP);

			e->setEndInterestPoint(iP);
			iP->addEndingEdge(e);
			//std::cout << "case2: trace multiple possibilities" << std::endl;

			//and use the last point to start 1 new edge for every possibility
			for (std::vector<cv::Point>::iterator it = neighbours.begin();
					it != neighbours.end(); ++it)
			{
				if (img.at<uchar>(it->y, it->x) == 255)
				{
					//std::cout << "case2: trace one of multiple possibilities: " << *it << std::endl;
					// tracing new edge from startPoint & add startPoint to new empty edge for overlapping
					Edge *etmp = new Edge(startPoint);
					etmp->addPreviousEdge(e);
					etmp->setStartInterestPoint(iP);
					iP->addStartingEdge(etmp);
					e->addFollowingEdge(etmp);
					traceEdgeAndInterestPoints(img, visited, *it, etmp, vE, vP);
				}
			}
		}
		else if (possibilities == 1)
		/*
		 * only one possibility
		 */
		{
			if (img.at<uchar>(neighbours[0].y, neighbours[0].x) == 255)
				{
					// mark pixel, when possibility of running a trace over there
					visited.at<uchar>(neighbours[0].y, neighbours[0].x) = 1;

					traceEdgeAndInterestPoints(img, visited, neighbours[0], e, vE, vP);
				}

		}
		else if (possibilities == 0)
		/*
		 * no possibility
		 */
		{
			//std::cout << "case2: push back edge because of end" << std::endl;
			// end recursion
			vE->push_back(e);
		}

	}
	else
	/*
	 * init case1 -> point was just found / no n-1 point in edge exists
	 * therefore start one tracing to every white point which has not been visited
	 * 		-> edges have to be inverted and merged
	 */
	{
		int possibilities = neighbours.size();

		// one pixel with more than one possibility
		if (possibilities > 2)
		/*
		 * multiple possibilities
		 */
		{
			// end edge we came from here
			InterestPoint* iP = new InterestPoint(startPoint);
			vP->push_back(iP);
			//std::cout << "case2: trace multiple possibilities" << std::endl;

			// and use the last point to start 1 new edge for every possibility
			for (std::vector<cv::Point>::iterator it = neighbours.begin();
					it != neighbours.end(); ++it)
			{
				if (img.at<uchar>(it->y, it->x) == 255)
				{
					//std::cout << "case2: trace one of multiple possibilities: " << *it << std::endl;
					// tracing new edge from startPoint & add startPoint to new empty edge for overlapping
					Edge *etmp = new Edge(startPoint);
					etmp->setStartInterestPoint(iP);
					iP->addStartingEdge(etmp);
					traceEdgeAndInterestPoints(img, visited, *it, etmp, vE, vP);
				}
			}
		}
		else if (possibilities == 2)
		/*
		 * two possibilities -> edge belongs together and goes into two opposite directions
		 */
		{
			Edge *e1 = new Edge(startPoint);
			Edge *e2 = new Edge(startPoint);

			traceEdgeAndInterestPoints(img, visited, neighbours[0], e1, vE, vP);
			traceEdgeAndInterestPoints(img, visited, neighbours[1], e2, vE, vP);

			e2->invertRecursive();
			//TODO handle InterestPoint as well

			// delete start point if necessary
			if (e2->getEndPoint() == e1->getStartPoint())
			{
				e1->deleteStartPoint();
			}
			e2->appendEdge(e1);

			vE->erase(std::find(vE->begin(), vE->end(), e1));

			if(getDistance(e2->getEndPoint(),e2->getStartPoint()) <= cv::sqrt(2)){
				e2->setIsLoop(true);
			}
			//TODO
			//delete(e1);
		}
		else if (possibilities == 1)
		/*
		 * only one possibility
		 */
		{
			if (img.at<uchar>(neighbours[0].y, neighbours[0].x) == 255)
			{
				//std::cout << "case2: trace only possibility " << *it << std::endl;
				// mark pixel, when possibility of running a trace over there
				visited.at<uchar>(neighbours[0].y, neighbours[0].x) = 1;
				traceEdgeAndInterestPoints(img, visited, neighbours[0], e, vE, vP);
			}
		}
		else if (possibilities == 0)
		/*
		 * no possibility
		 */
		{
			//std::cout << "case2: push back edge because of end" << std::endl;
			// end recursion
			vE->push_back(e);
		}
	}

	neighbours.clear();
}

std::vector<cv::Point> EdgeProcessor::get8Neighbours(cv::Point p)
{
	std::vector<cv::Point> v;

	v.push_back(cv::Point(p.x - 1, p.y - 1));
	v.push_back(cv::Point(p.x, p.y - 1));
	v.push_back(cv::Point(p.x + 1, p.y - 1));
	v.push_back(cv::Point(p.x + 1, p.y));
	v.push_back(cv::Point(p.x + 1, p.y + 1));
	v.push_back(cv::Point(p.x, p.y + 1));
	v.push_back(cv::Point(p.x - 1, p.y + 1));
	v.push_back(cv::Point(p.x - 1, p.y));

	return v;
}

std::vector<cv::Point> EdgeProcessor::getNeighbours(cv::Point lastPoint,
		cv::Point previousPoint)
{
	std::vector<cv::Point> v;
	// add neighbour in direction of edge first
	v.push_back(lastPoint + (lastPoint - previousPoint));

	v.push_back(cv::Point(lastPoint.x - 1, lastPoint.y - 1));
	v.push_back(cv::Point(lastPoint.x, lastPoint.y - 1));
	v.push_back(cv::Point(lastPoint.x + 1, lastPoint.y - 1));
	v.push_back(cv::Point(lastPoint.x + 1, lastPoint.y));
	v.push_back(cv::Point(lastPoint.x + 1, lastPoint.y + 1));
	v.push_back(cv::Point(lastPoint.x, lastPoint.y + 1));
	v.push_back(cv::Point(lastPoint.x - 1, lastPoint.y + 1));
	v.push_back(cv::Point(lastPoint.x - 1, lastPoint.y));

	//delete the redundant point
	v.erase(std::find(v.begin()+1, v.end(), lastPoint + (lastPoint - previousPoint)));

	return v;
}

float EdgeProcessor::getDistance(cv::Point p1, cv::Point p2){
	cv::Point diff = p1 - p2;
	return cv::sqrt(diff.x*diff.x + diff.y*diff.y);
}

void EdgeProcessor::selectionOfNeighbours(cv::Point startPoint, std::vector<cv::Point> &neighbours, cv::Mat &img, cv::Mat &visited)
{
	// initialize lastWhiteNeighbour with position far outside the image because lastWhiteNeighbour
	// is used to check if two neighbours are direct neighbours (distance == 1)
	cv::Point lastWhiteNeighbour(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());

	// calculate number of possibilities to forward the recursion
	for (std::vector<cv::Point>::iterator it = neighbours.begin();
			it != neighbours.end();)
	{	
		if (img.at<uchar>(it->y, it->x) == 255
				&& visited.at<uchar>(it->y, it->x) == 0)
		/*
		 * check if neighbours are white and have not been visited before
		 */
		{
			/* check if the neighbours are direct neighbours (distance == 1.0)
			 * if so: take pixel which is 4-neighbour to startPoint (currently analyzed pixel)
			 */
			if (getDistance(*it, lastWhiteNeighbour) == 1.0)
			{
				if (getDistance(startPoint, (*it)) == 1.0)
				{
					// withdraw lastWhiteNeighbour
					visited.at<uchar>(lastWhiteNeighbour.y, lastWhiteNeighbour.x) = 0;
					neighbours.erase(std::find(neighbours.begin(), neighbours.end(), lastWhiteNeighbour));

					lastWhiteNeighbour = *(it-1);
					visited.at<uchar>(lastWhiteNeighbour.y, lastWhiteNeighbour.x) = 1;
					
					//std::cout << "Special case of pixel neighbourhood" << std::endl;
				}
				else
				{
					// erase from neighbours and set iterator to following element
					it = neighbours.erase(it);
				}
			}
			else
			{
				lastWhiteNeighbour = cv::Point(it->x, it->y);
				// mark as visited so that no following point in the recursion goes there
				visited.at<uchar>(it->y, it->x) = 1;

				++it;
			}
		}
		else
		{
			// remove neighbour from vector so that pixel is followed in next recursion step
			it = neighbours.erase(it);
		}
	}
}
