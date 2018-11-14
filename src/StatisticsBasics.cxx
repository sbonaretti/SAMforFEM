/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <StatisticsBasics.h>
#include <iostream>
#include <algorithm>


namespace statistics{

	// constructor
	StatisticsBasics::StatisticsBasics(){

		_standardDeviation = 0.0;
	}
	// destructor
	StatisticsBasics::~StatisticsBasics(){
	}

	
	// member funcions
	void StatisticsBasics::CalculateStandardDeviation(){
		
		_standardDeviation = 0.0;

		for (int i=0; i<_vector.size(); i++)
			_standardDeviation += (_vector(i) - _vector.mean())*(_vector(i) - _vector.mean());
		_standardDeviation /= (_vector.size()-1);
		_standardDeviation = std::sqrt (_standardDeviation);

	}

	void StatisticsBasics::CalculateBoxPlot(){

		// vnl to std
		std::vector<double> vector;
		for (int i=0; i<_vector.size(); i++)
			vector.push_back(_vector(i)); 
		
		// sort
		sort (vector.begin(), vector.end()); 
		
		// first percentile (25%)
		double firstPercentileValue = 0.0; double firstPercentileIndex = 0.0; 
		CalculatePercentile(vector, 25, firstPercentileIndex, firstPercentileValue);
		std::cout << "firstPercentileIndex: " << firstPercentileIndex << std::endl;
		std::cout << "firstPercentileValue: " << firstPercentileValue << std::endl;

		// second percentile (50%)
		double secondPercentileValue = 0.0; double secondPercentileIndex = 0.0; 
		CalculatePercentile(vector, 50, secondPercentileIndex, secondPercentileValue);
		std::cout << "secondPercentileValue: " << secondPercentileValue << std::endl;

		// third percentile (75%)
		double thirdPercentileValue = 0.0; double thirdPercentileIndex = 0.0; 
		CalculatePercentile(vector, 75, thirdPercentileIndex, thirdPercentileValue);
		std::cout << "thirdPercentileIndex: " << thirdPercentileIndex << std::endl;
		std::cout << "thirdPercentileValue: " << thirdPercentileValue << std::endl;

		// from here down - not tested
		// lower value 
		double lowerIndex = 1.5 * (thirdPercentileIndex-firstPercentileIndex) - firstPercentileIndex;
		std::cout << "lowerIndex: " << lowerIndex << std::endl;
		double lowerValue = vector.at(lowerIndex);
		//double lowerValue = 1.5 * (thirdPercentileValue-firstPercentileValue) - firstPercentileValue;
		std::cout << "lowerValue: " << lowerValue << std::endl;

		// upper value
		double upperIndex = 1.5 * (thirdPercentileIndex-firstPercentileIndex) + thirdPercentileIndex;
		std::cout << "upperIndex: " << upperIndex << std::endl;
		double upperValue = vector.at(upperIndex);
		std::cout << "upperValue: " << upperValue << std::endl;

		// outliers
		std::vector<double> lowerOutliers;
		std::vector<double> upperOutliers;
		
		for (int i=0; i<vector.size(); i++){
			if (vector.at(i) < lowerValue)
				lowerOutliers.push_back(vector.at(i));
			else if (vector.at(i) > upperIndex)
				upperOutliers.push_back(vector.at(i));
		}
		
		std::cout << "lower outliers" << std::endl;
		for (int i=0; i<lowerOutliers.size(); i++)
			std::cout << lowerOutliers.at(i) << ' ';
		std::cout << std::endl;

		std::cout << "upper outliers" << std::endl;
		for (int i=0; i<upperOutliers.size(); i++)
			std::cout << upperOutliers.at(i) << ' ';
		std::cout << std::endl;

		
	}
	
	// private
	void StatisticsBasics::CalculatePercentile(std::vector<double> vector, double percentage, double & percentileIndex, double & percentileValue){
		
		int n = vector.size();
		
		// index
		percentileIndex = (n+1) * percentage/100;
		double lowerIndex = floor(percentileIndex) -1; // array starts at 0
		double upperIndex = floor(percentileIndex); // +1 -1
		
		// value
		percentileValue = vector.at(lowerIndex) + (vector.at(upperIndex)-vector.at(lowerIndex)) * (percentage/100);



		







	}
}