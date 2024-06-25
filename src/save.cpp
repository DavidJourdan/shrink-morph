#include "save.h"

#include <fstream>
#include <nlohmann/json.hpp>

// Function to update and save advanced parameters to a JSON file
void updateAndSaveAdvancedParams(double lambda1, double lambda2, double thickness, double deltaLambda, int n_layers, int n_iter, double lim) {
    // Create a JSON object to store advanced parameters
    nlohmann::json advancedParams;
    advancedParams["lambda1"] = lambda1;
    advancedParams["lambda2"] = lambda2;
    advancedParams["thickness"] = thickness;
    advancedParams["deltaLambda"] = deltaLambda;
    advancedParams["n_layers"] = n_layers;
    advancedParams["n_iter"] = n_iter;
    advancedParams["lim"] = lim;

    // Serialize JSON to file
    std::ofstream outputFile("advanced_params.json");
    outputFile << advancedParams.dump(4); // dump with indentation for readability
    outputFile.close();
}

// Function to load advanced parameters from a JSON file
void loadAdvancedParams(double& lambda1, double& lambda2, double& thickness, double& deltaLambda, int& n_layers, int& n_iter, double& lim) {
    // Read the JSON file
    std::ifstream inputFile("advanced_params.json");
    if (!inputFile.is_open()) {
        // Material data
        lambda1 = 0.58;
        lambda2 = 1.08;
        thickness = 1.218;
        deltaLambda = 0.0226764665509417;
        n_layers = 10;
        lim = 1e-6;
        n_iter = 1000;
        return;
    }

    // Parse JSON from file
    nlohmann::json advancedParams;
    inputFile >> advancedParams;

    // Extract parameters from JSON
    lambda1 = advancedParams["lambda1"];
    lambda2 = advancedParams["lambda2"];
    thickness = advancedParams["thickness"];
    deltaLambda = advancedParams["deltaLambda"];
    n_layers = advancedParams["n_layers"];
    n_iter = advancedParams["n_iter"];
    lim = advancedParams["lim"];

    // Close the file
    inputFile.close();
}

