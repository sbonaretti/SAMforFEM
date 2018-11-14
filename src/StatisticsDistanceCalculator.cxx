/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <StatisticsDistanceCalculator.h>
#include <iostream>
#include <math.h>

#include <StatisticsBasics.h>

namespace statistics{

	// constructor
	StatisticsDistanceCalculator::StatisticsDistanceCalculator(){
	}
	// destructor
	StatisticsDistanceCalculator::~StatisticsDistanceCalculator(){
	}

	// "gets" accessors
	double StatisticsDistanceCalculator::GetDistanceAverage(){

		_distanceAverage = _distancesVector.mean();
		
		return _distanceAverage;
	}
	
	double StatisticsDistanceCalculator::GetDistanceStdDev(){

		
		StatisticsBasics* stddevCalculator = new StatisticsBasics;
		stddevCalculator->SetVector(_distancesVector);
		stddevCalculator->CalculateStandardDeviation();
		_distanceStdDev = stddevCalculator->GetStandardDeviation();
				
		return _distanceStdDev;
	
	}

	double StatisticsDistanceCalculator::GetDistanceStdError(){

		for (int i=0; i<_distancesVector.size(); i++)
			_distanceStdDev += (_distancesVector(i) - _distanceAverage)*(_distancesVector(i) - _distanceAverage);
		_distanceStdDev /= (_distancesVector.size()-1);
		_distanceStdDev = std::sqrt (_distanceStdDev);
		
		_distanceStdError = _distanceStdDev / std::sqrt(double(_distancesVector.size()));
				
		return _distanceStdError;
	}


	// methods
	void StatisticsDistanceCalculator::CalculateEuclideanDistance(){

		_distancesVector.set_size(_vectorOne.size());

		for (int i=0; i<_distancesVector.size(); i++)
			_distancesVector(i) = std::sqrt((_vectorOne(i) - _vectorTwo(i)) * (_vectorOne(i) - _vectorTwo(i)));
			
	}
	
	void StatisticsDistanceCalculator::CalculateEuclideanDistancePoints(){
		
		_distancesVector.set_size(_vectorOne.size()/3);

		for (int i=0; i<_distancesVector.size(); i++){
			_distancesVector(i) = std::sqrt (	pow(_vectorOne(3*i) - _vectorTwo(3*i), 2) +
												pow(_vectorOne(3*i +1) - _vectorTwo(3*i +1), 2) +
												pow(_vectorOne(3*i +2) - _vectorTwo(3*i +2), 2) );
			//std::cout << _vectorOne(3*i) << ' ' << _vectorTwo(3*i) << ' ' << pow(_vectorOne(3*i) - _vectorTwo(3*i), 2) << std::endl;
			//std::cout << _vectorOne(3*i +1) << ' ' << _vectorTwo(3*i +1) << ' ' << pow(_vectorOne(3*i +1) - _vectorTwo(3*i +1), 2) << std::endl;
			//std::cout << _vectorOne(3*i +2) << ' ' << _vectorTwo(3*i +2) << ' ' << pow(_vectorOne(3*i +2) - _vectorTwo(3*i +2), 2) << std::endl;

			//std::cout << _distancesVector(i) << std::endl;
		}
	
	}
		

}