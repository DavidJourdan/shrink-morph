#pragma once

// Function to update and save advanced parameters to a JSON file
void updateAndSaveAdvancedParams(double lambda1, double lambda2, double thickness, double deltaLambda, int n_layers, int n_iter, double lim);

// Function to load advanced parameters from a JSON file
void loadAdvancedParams(double& lambda1, double& lambda2, double& thickness, double& deltaLambda, int& n_layers, int& n_iter, double& lim);
