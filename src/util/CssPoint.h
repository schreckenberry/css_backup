#ifndef CSSPOINT_H
#define CSSPOINT_H

/*
 * Structure to save specific points from the 3D CSS,
 * especially signatures
 */ 
struct CssPoint {
	int u;			//!< u coordinate
	float sigma;	//!< corresponding sigma (scale) value
	float kappa;	//!< corresponding curvature value
};

#endif // CSSPOINT_H
