/*
 * This file is part of TinyAD and released under the MIT license.
 * Author: Patrick Schmidt
 */
#include "LocalGlobalSolver.h"
#include "constants.h"
#include "functions.h"
#include "newton.h"
#include "parameterization.h"
#include "path_extraction.h"
#include "simulation_utils.h"
#include "stretch_angles.h"
#include "timer.h"

#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <igl/file_dialog_open.h>
#include <igl/loop.h>
#include <igl/readOBJ.h>
#include <polyscope/curve_network.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <thread>

int main(int argc, char* argv[])
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  // Initialize polyscope
  polyscope::init();
  polyscope::view::upDir = polyscope::UpDir::NegYUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Planar;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  polyscope::options::openImGuiWindowForUserCallback = false;
  polyscope::loadColorMap("twilight", DATA_PATH_STR "twilight_colormap.png");

  // Some state about imgui windows to stack them
  float imguiStackMargin = 10;
  float rightWindowsWidth = 370;

  // Material data
  double lambda1 = 0.58;
  double lambda2 = 1.08;
  double thickness = 1.218;
  double deltaLambda = 0.0226764665509417;

  // Load a mesh in OBJ format
  std::string filename;
  if(argc < 2)
    filename = igl::file_dialog_open();
  else
    filename = std::string(DATA_PATH_STR) + std::string(argv[1]) + ".obj";
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::readOBJ(filename, V, F))
  {
    std::cout << "File " << filename << " not found\n";
    return 0;
  }
  // resize mesh to a 100mm width
  V *= 100 / (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();
  V.col(1) *= -1;

  // subdivide until number of mesh faces is at least nFmin
  int nFmin = 1000;
  if(argc == 3)
    nFmin = std::stoi(argv[2]);
  while(F.rows() < nFmin)
  {
    Eigen::MatrixXd tempV = V;
    Eigen::MatrixXi tempF = F;
    igl::loop(tempV, tempF, V, F);
  }

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  VertexPositionGeometry geometry(mesh, V);

  // Run local-global parameterization algorithm
  Timer paramTimer("Parameterization");
  LocalGlobalSolver LGsolver(V, F);
  Eigen::MatrixXd P = localGlobal(V, F, lambda1, lambda2, LGsolver);
  paramTimer.stop();

  // Add mesh and initial param to polyscope viewer
  polyscope::registerSurfaceMesh2D("Param", P, F)->setEdgeWidth(0.0);
  polyscope::getSurfaceMesh("Param")->addFaceScalarQuantity("sigma1", LGsolver.s1);
  polyscope::getSurfaceMesh("Param")->addFaceScalarQuantity("sigma2", LGsolver.s2);
  polyscope::getSurfaceMesh("Param")
      ->addFaceScalarQuantity("stretch orientation", LGsolver.stretchAngles())
      ->setColorMap("twilight")
      ->setMapRange({-PI / 2, PI / 2})
      ->setEnabled(true);

  // callback function for ImGui interface
  polyscope::state::userCallback = [&]() {
    // parameters modifiable in the GUI
    static double wD = 0.1;
    static double lim = 1e-6;
    static int n_iter = 1000;
    static double rotAngle = 0;
    static double curvature = deltaLambda / thickness / lambda1;

    ImGui::PushID("user_callback");
    ImGui::SetNextWindowPos(
        ImVec2(polyscope::view::windowWidth - (rightWindowsWidth + imguiStackMargin), imguiStackMargin));
    ImGui::SetNextWindowSize(ImVec2(rightWindowsWidth, 0.));

    ImGui::Begin("Command UI", nullptr);

    ImGui::PushItemWidth(100);

    if(ImGui::Button("Smoothing"))
    {
      paramTimer.start();
      auto func = parameterizationFunction(geometry, wD, lambda1, lambda2);

      // optimize the parameterization function with Newton's method
      newton(geometry, P, func, n_iter, lim);
      paramTimer.stop();

      // Show resulting parametrization
      polyscope::getSurfaceMesh("Param")->updateVertexPositions2D(P);

      // compute data from SVDs and display them
      LGsolver.solveOneStep(P, 1 / lambda2, 1 / lambda1);
      // center P
      P.rowwise() -= P.colwise().sum() / P.rows();

      polyscope::getSurfaceMesh("Param")->addFaceScalarQuantity("sigma1", LGsolver.s1);
      polyscope::getSurfaceMesh("Param")->addFaceScalarQuantity("sigma2", LGsolver.s2);
      polyscope::getSurfaceMesh("Param")
          ->addFaceScalarQuantity("stretch orientation", LGsolver.stretchAngles())
          ->setColorMap("twilight")
          ->setMapRange({-PI / 2, PI / 2})
          ->setEnabled(true);
    }
    ImGui::SameLine();
    ImGui::PushItemWidth(40);
    ImGui::InputDouble("wD", &wD, 0, 0, "%.1f");
    ImGui::PopItemWidth();

    if(ImGui::TreeNode("Advanced"))
    {
      ImGui::InputDouble("lambda1", &lambda1, 0, 0, "%.1e");
      ImGui::InputDouble("lambda2", &lambda2, 0, 0, "%.1e");
      ImGui::InputDouble("Thickness", &thickness, 0, 0, "%.1e");
      if (ImGui::InputDouble("Curvature", &curvature, 0, 0, "%.1e"))
      {
        deltaLambda = thickness * lambda1 * curvature;
      }
      ImGui::InputInt("Iterations", &n_iter);
      ImGui::InputDouble("Limit", &lim, 0, 0, "%.1e");
      if(ImGui::Button("Rotate"))
      {
        // run ARAP
        LGsolver.solve(P, 1., 1.);

        // center and rotate
        P.rowwise() -= P.colwise().sum() / P.rows();
        P = P * Eigen::Rotation2D<double>(-rotAngle / 180 * PI).toRotationMatrix();

        // change aspect ratio
        P.col(0) /= lambda1;
        P.col(1) /= lambda2;

        // run Local-Global algorithm
        LGsolver.solve(P, 1 / lambda2, 1 / lambda1);
        // center again
        P.rowwise() -= P.colwise().sum() / P.rows();

        polyscope::getSurfaceMesh("Param")->updateVertexPositions2D(P);
      }
      ImGui::SameLine();
      ImGui::PushItemWidth(40);
      ImGui::InputDouble("Angle", &rotAngle, 0, 0, "%.0f");
      ImGui::PopItemWidth();
      ImGui::TreePop();
    }

    float lastWindowHeightUser = imguiStackMargin + ImGui::GetWindowHeight();
    ImGui::End();
    ImGui::PopID();

    ImGui::SetNextWindowPos(ImVec2(polyscope::view::windowWidth - (rightWindowsWidth + imguiStackMargin),
                                   2 * imguiStackMargin + lastWindowHeightUser));
    ImGui::SetNextWindowSize(ImVec2(rightWindowsWidth, 0.));

    ImGui::Begin("Parameterization", nullptr);
    ImGui::Text("Computation time: %.2f s.", paramTimer.seconds());
    if(ImGui::Button("Next"))
      polyscope::popContext();
    ImGui::End();
  };
  polyscope::show();
  polyscope::removeAllStructures();

  // Resize to 100mm width
  const double scaleFactor = 100 / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
  V *= scaleFactor;
  P *= scaleFactor;
  geometry.inputVertexPositions *= scaleFactor;
  geometry.refreshQuantities();

  // Init mesh and compute theta1
  polyscope::registerSurfaceMesh("Input mesh", V, F);
  FaceData<Eigen::Matrix2d> MrInv = precomputeSimData(mesh, P, F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, F, MrInv);
  polyscope::getSurfaceMesh("Input mesh")
      ->addFaceScalarQuantity("theta1", theta1)
      ->setColorMap("twilight")
      ->setMapRange({-PI / 2, PI / 2});

  // initialize theta2 (stored on vertices for convenience and speed)
  VertexData<double> theta2(mesh, 0);

  // Find face closest to mesh center and fix its vertices
  std::vector<int> fixedIdx = findCenterFaceIndices(P, F);

  // Save the input mesh DOFs
  Eigen::VectorXd xTarget(V.size());
  for(int i = 0; i < V.rows(); ++i)
    for(int j = 0; j < 3; ++j)
      xTarget(3 * i + j) = V(i, j);

  // Init timers for GUI
  Timer optimTimer("Directions optimization");
  optimTimer.stop();
  Timer trajTimer("Trajectory optimization");
  trajTimer.stop();

  polyscope::state::userCallback = [&]() {
    // App state
    static bool optimRun = false;
    static bool trajectoriesRun = false;
    // threads
    static std::unique_ptr<std::thread> thr;
    static bool threadFinished = false;
    // GUI parameters
    static double wM = 0.01;
    static double wL = 0.01;
    static double timeLimit = 5;
    static double E1 = 10;
    static double width = 100;
    static double lim = 1e-6;
    static int n_iter = 1000;
    static bool displayTravel = false;

    ImGui::PushID("user_callback");
    ImGui::SetNextWindowPos(
        ImVec2(polyscope::view::windowWidth - (rightWindowsWidth + imguiStackMargin), imguiStackMargin));
    ImGui::SetNextWindowSize(ImVec2(rightWindowsWidth, 0.));

    ImGui::Begin("Command UI", nullptr);

    ImGui::PushItemWidth(100);
    if(ImGui::InputDouble("Width", &width, 0, 0, "%.0f"))
    {
      if(width > 0)
      {
        // Resize with respect to width
        const double scaleFactor = width / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
        V *= scaleFactor;
        P *= scaleFactor;
        geometry.inputVertexPositions *= scaleFactor;
        xTarget *= scaleFactor;
        geometry.refreshQuantities();
        MrInv = precomputeSimData(mesh, P, F);
      }
    }

    if(ImGui::Button("Simulation"))
    {
      // Load input mesh to Polyscope (not displayed by default)
      polyscope::removeAllStructures();
      Eigen::MatrixXd VTarget = xTarget.reshaped<Eigen::RowMajor>(V.rows(), 3);
      polyscope::registerSurfaceMesh("Input mesh", VTarget, F)->setEnabled(false);

      // Run the optimization in a separate thread
      thr = std::make_unique<std::thread>([&]() {
        threadFinished = false;
        polyscope::options::alwaysRedraw = true;
        optimTimer.start();
        optimRun = true;
        // Define simulation function
        auto func = simulationFunction(geometry, MrInv, theta1, theta2, E1, lambda1, lambda2, deltaLambda, thickness);

        // (Projected) Newton optimization
        newton(geometry, V, func, n_iter, lim, true, fixedIdx, [&](const Eigen::VectorXd& X) {
          // convert X to V for visualizing individual iterations
          for(int i = 0; i < V.rows(); ++i)
            for(int j = 0; j < 3; ++j)
              V(i, j) = X(3 * i + j);
        });

        optimTimer.stop();
        threadFinished = true;
        polyscope::options::alwaysRedraw = false;
      });

      // Display simulated mesh
      polyscope::registerSurfaceMesh("Simulation", V, F);
      polyscope::getSurfaceMesh("Simulation")
          ->addFaceScalarQuantity("theta1", theta1)
          ->setColorMap("twilight")
          ->setMapRange({-PI / 2, PI / 2});
      polyscope::getSurfaceMesh("Simulation")->addVertexScalarQuantity("theta2", theta2);
      polyscope::view::flyToHomeView();
    }

    if(ImGui::Button("Directions optimization"))
    {
      // Display simulated mesh (will be updated at each iteration)
      polyscope::removeAllStructures();
      polyscope::registerSurfaceMesh("Simulation", V, F);
      polyscope::getSurfaceMesh("Simulation")
          ->addFaceScalarQuantity("theta1", theta1)
          ->setColorMap("twilight")
          ->setMapRange({-PI / 2, PI / 2});
      polyscope::view::flyToHomeView();

      // Load input mesh to Polyscope (not displayed by default)
      Eigen::MatrixXd VTarget = xTarget.reshaped<Eigen::RowMajor>(V.rows(), 3);
      polyscope::registerSurfaceMesh("Input mesh", VTarget, F)->setEnabled(false);

      // Run the optimization in a separate thread
      thr = std::make_unique<std::thread>([&]() {
        threadFinished = false;
        polyscope::options::alwaysRedraw = true;
        optimTimer.start();
        optimRun = true;

        // Define optimization function
        auto adjointFunc = adjointFunction(geometry, F, MrInv, theta1, E1, lambda1, lambda2, deltaLambda, thickness);

        // Optimize this energy function using SGN [Zehnder et al. 2021]
        sparse_gauss_newton(
            geometry, V, xTarget, MrInv, theta1, theta2, adjointFunc, fixedIdx, n_iter, lim, wM, wL, E1, lambda1,
            lambda2, deltaLambda, thickness, [&](const Eigen::VectorXd& X) {
                              // convert X to V for visualizing individual iterations
                              for(int i = 0; i < V.rows(); ++i)
                                for(int j = 0; j < 3; ++j)
                                  V(i, j) = X(3 * i + j);
                            });

        optimTimer.stop();

        threadFinished = true;
        polyscope::options::alwaysRedraw = false;
      });
    }

    if(polyscope::options::alwaysRedraw)
    {
      polyscope::getSurfaceMesh("Simulation")->updateVertexPositions(V);
    }

    if(threadFinished && thr->joinable())
    {
      thr->join();

      // Display theta2
      polyscope::getSurfaceMesh("Simulation")->addVertexScalarQuantity("theta2", theta2);

      // Display distance with target mesh
      Eigen::MatrixXd VTarget = xTarget.reshaped<Eigen::RowMajor>(V.rows(), 3);
      Eigen::VectorXd d = (V - VTarget).cwiseProduct((V - VTarget)).rowwise().sum();
      d = d.array().sqrt();
      std::cout << "Avg distance = "
                << 100 * d.sum() / d.size() / (VTarget.colwise().maxCoeff() - VTarget.colwise().minCoeff()).norm()
                << "\n";
      std::cout << "Max distance = "
                << 100 * d.lpNorm<Eigen::Infinity>() /
                       (VTarget.colwise().maxCoeff() - VTarget.colwise().minCoeff()).norm()
                << "\n";
      polyscope::getSurfaceMesh("Simulation")->addVertexScalarQuantity("Distance", d)->setEnabled(true);
    }

    if(ImGui::Button("Generate trajectories"))
    {
      // Center flat mesh at 0
      Eigen::RowVector2d center = P.colwise().sum() / P.rows();
      for(int i = 0; i < P.rows(); ++i)
        P.row(i) = P.row(i) - center;

      polyscope::removeAllStructures();
      polyscope::view::style = polyscope::view::NavigateStyle::Planar;
      polyscope::registerSurfaceMesh2D("Param", P, F)->setEnabled(false);
      polyscope::view::flyToHomeView();

      trajTimer.start();
      trajectoriesRun = true;
      stripePattern(geometry, V, P, F, theta2, filename, timeLimit);
      trajTimer.stop();
    }
    ImGui::SameLine();
    ImGui::PushItemWidth(40);
    ImGui::InputDouble("(s) Time limit per layer", &timeLimit, 0, 0, "%.1f");
    ImGui::PopItemWidth();

    if(ImGui::TreeNode("Advanced parameters"))
    {
      ImGui::InputInt("Iterations", &n_iter);
      ImGui::InputDouble("Limit", &lim, 0, 0, "%.1e");
      ImGui::InputDouble("E1 / E2 ratio", &E1, 0, 0, "%.0f");

      ImGui::PushItemWidth(40);
      ImGui::InputDouble("wM", &wM, 0, 0, "%.2f");
      ImGui::SameLine();
      ImGui::InputDouble("wL", &wL, 0, 0, "%.2f");

      ImGui::TreePop();
    }

    float lastWindowHeightUser = imguiStackMargin + ImGui::GetWindowHeight();
    ImGui::End();
    ImGui::PopID();

    if(optimRun)
    {
      ImGui::SetNextWindowPos(ImVec2(polyscope::view::windowWidth - (rightWindowsWidth + imguiStackMargin),
                                     2 * imguiStackMargin + lastWindowHeightUser));
      ImGui::SetNextWindowSize(ImVec2(rightWindowsWidth, 0.));

      ImGui::Begin("Optimization", nullptr);
      if(polyscope::options::alwaysRedraw)
        ImGui::TextUnformatted("Computing...");
      else
      {
        double s = optimTimer.seconds();
        if(s < 60)
          ImGui::Text("Computation time: %.2f s.", s);
        else
        {
          double r = std::fmod(s, 60);
          double q = (s - r) / 60;
          ImGui::Text("Computation time: %.0f min %.0f s.", q, r);
        }
      }
      lastWindowHeightUser += imguiStackMargin + ImGui::GetWindowHeight();
      ImGui::End();
    }

    static int lastLayer = 10;
    if(trajectoriesRun)
    {
      ImGui::SetNextWindowPos(ImVec2(polyscope::view::windowWidth - (rightWindowsWidth + imguiStackMargin),
                                     2 * imguiStackMargin + lastWindowHeightUser));
      ImGui::SetNextWindowSize(ImVec2(rightWindowsWidth, 0.));

      ImGui::Begin("Paths preview", nullptr);
      if(polyscope::options::alwaysRedraw)
        ImGui::TextUnformatted("Computing...");
      else
      {
        auto onChangeLayer = [&]() {
          for(int i = 1; i <= 10; ++i)
          {
            if(i == lastLayer)
            {
              polyscope::getCurveNetwork("Layer " + std::to_string(i))->setEnabled(true);
              polyscope::getCurveNetwork("Travel " + std::to_string(i))->setEnabled(displayTravel);
            }
            else
            {
              polyscope::getCurveNetwork("Layer " + std::to_string(i))->setEnabled(false);
              polyscope::getCurveNetwork("Travel " + std::to_string(i))->setEnabled(false);
            }
          }
        };
        if(ImGui::SliderInt("Layer", &lastLayer, 1, 10))
        {
          onChangeLayer();
        }
        if(ImGui::Checkbox("Travel", &displayTravel))
        {
          onChangeLayer();
        }
        double s = trajTimer.seconds();
        if(s < 60)
          ImGui::Text("Computation time: %.2f s.", s);
        else
        {
          double r = std::fmod(s, 60);
          double q = (s - r) / 60;
          ImGui::Text("Computation time: %.0f min %.0f s.", q, r);
        }
      }
      lastWindowHeightUser += imguiStackMargin + ImGui::GetWindowHeight();
      ImGui::End();
    }
  };

  polyscope::view::style = polyscope::view::NavigateStyle::Turntable;
  polyscope::view::flyToHomeView();
  polyscope::show();
}
