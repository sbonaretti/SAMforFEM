/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef STATISTICSBASICS_H
#define STATISTICSBASICS_H

#include <vnl/vnl_vector.h>

#include <vector>


namespace statistics{

	/**
	* Class for the calculation of basic statistics
	*/

	class StatisticsBasics{

	public:

		// constructor/destructor
		StatisticsBasics();
		~StatisticsBasics();

		// accessors
		/**
		* Sets the vector
		*/
		void SetVector (vnl_vector<double> vector) {_vector = vector;}
		
		/**
		* Gets the standard deviation
		*/
		double GetStandardDeviation() {return _standardDeviation;}
		/**
		* Gets the boxplot
		*/
		vnl_vector<double> GetBoxPlot() {return _boxPlot;}

		
		// methods
		/**
		* Calculates the standard deviation
		*/
		void CalculateStandardDeviation();
		
		/**
		* Calculates the box plot values
		*/
		void CalculateBoxPlot();
		

	private:
		/**
		* Calculates the percentile using the linear interpolation
		*/
		void CalculatePercentile(std::vector<double> vector, double percentage, double & percentileIndex, double & percentileValue);


	protected:
		vnl_vector<double> _vector;
		double _standardDeviation;
		vnl_vector<double> _boxPlot;
		double _percentage;
		double _percentile;

	};

}
#endif //STATISTICSBASICS_H
		