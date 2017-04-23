#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if (estimations.size() == 0 || estimations.size() != ground_truth.size())
  {
  	std::cout << "Invalid estimation or ground truth " << std::endl;
  	return rmse;
  } 
  for (int i=0; i< estimations.size(); i++)
  {
  	VectorXd residual = estimations[i] - ground_truth[i];
  	
  	//Coefficient wise mutiplication
  	residual = residual.array() * residual.array();
  	rmse += residual;  	
  }
  
  //Calculate the mean
  rmse = rmse / estimations.size();
  
  //Calculte the square root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}
