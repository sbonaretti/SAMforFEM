/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef STATISTICSDISTANCECALCULATOR_H
#define STATISTICSDISTANCECALCULATOR_H

#include <vnl/vnl_vector.h>


namespace statistics{

	/**
	* Class for the calculation of distances
	*/

	class StatisticsDistanceCalculator{

	public:

		// constructor/destructor
		StatisticsDistanceCalculator();
		~StatisticsDistanceCalculator();

		// accessors
		/**
		* Sets the first vector
		*/
		void SetVectorOne (vnl_vector<double> vector) {_vectorOne = vector;}
		/**
		* Sets the second vector
		*/
		void SetVectorTwo (vnl_vector<double> vector) {_vectorTwo = vector;}
		/**
		* Gets the distance vector
		*/
		vnl_vector<double> GetDistanceVector(){return _distancesVector;}
		/**
		* Gets the distance average
		*/
		double GetDistanceAverage();
		/**
		* Gets the distance standard deviation. Call GetDistanceAverage() before it.
		*/
		double GetDistanceStdDev();
		/**
		* Gets the distance standard error. Call GetDistanceAverage() and GetDistanceStdError() before it.
		*/
		double GetDistanceStdError();

		// methods
		/**
		* Calculates the Euclidean distance (square root of the square of the differences) for each element
		*/
		void CalculateEuclideanDistance();
		/**
		* Calculates the Euclidean distance for each point (x,y,z) of the vector
		*/
		void CalculateEuclideanDistancePoints();


	protected:
		vnl_vector<double> _vectorOne;
		vnl_vector<double> _vectorTwo;
		vnl_vector<double> _distancesVector;
		double _distanceAverage;
		double _distanceStdDev;
		double _distanceStdError;

	

	};
}
#endif //STATISTICSDISTANCECALCULATOR_H