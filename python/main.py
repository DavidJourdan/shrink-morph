import shrink_morph_py
import numpy as np
import polyscope as ps
import polyscope.imgui as gui
import togcode
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
filename = filedialog.askopenfilename(defaultextension=".obj", filetypes=[("Wavefront OBJ", "*.obj")])
V, F = shrink_morph_py.read_from_OBJ(filename)

class ShrinkMorph:
  lambda1 = 0.58
  lambda2 = 1.08
  wD = 1e-5
  E1 = 10
  thickness = 1.218
  deltaLambda = 0.0226764665509417
  n_layers = 10
  lim = 1e-6
  n_iter = 1000
  width = 200
  wM = 0.01
  wL = 0.01
  num_rectangles = 5
  num_layers = 10
  layer_height = 0.08
  with_smoothing = False
  printer_profile = "Prusa_MK3S"
  printers_list = [
    'Anet_A8', 
    'Anycubic_Mega_Zero',
    'Artillery_SW_X1', 
    'Creality_K1_Max', 
    'CR10',
    'CR10S_Pro', 
    'CR10_S5', 
    'E3D_Toolchanger', 
    'Ender2', 
    'Ender3', 
    'Flsun_SR', 
    'Kingroon_P3', 
    'Kywoo3D_Tycoon'
    'Lutum_4.6', 
    'Micro_Delta_Rework', 
    'Monoprice_select_mini_v2', 
    'Pharaoh_xd', 
    'Prusa_MK2S', 
    'Prusa_MK3S', 
    'Replicator2', 
    'Replicator2x', 
    'Snapmaker', 
    'Strateo3D', 
    'Ultimaker2', 
    'Ultimaker3', 
    'UltimakerS3', 
    'UltimakerS5', 
    'UltimakerS7', 
    'Voron_V0', 
    'Volumic_Stream30_Pro_MK2', 
    'Wasp_2040_Pro', 
  ]
  printer = togcode.Printer(printer_profile)

  # GUI variables
  leave = True

  # Display printer buildplate
  def display_buildplate(self):
    build_vert = np.array([[-1,-1,-0.1], [1, -1, -0.1], [1, 1, -0.1], [-1, 1, -0.1]])
    build_vert[:, 0] *= self.printer.bed_size[0] / 2
    build_vert[:, 1] *= self.printer.bed_size[1] / 2
    build_face = np.array([[0, 1, 2, 3]])

    ps.register_surface_mesh("Buildplate", build_vert, build_face, color=(0.95, 0.95, 0.95), edge_width=5, edge_color=(0.5, 0.5, 0.5), material="flat")

  def param_screen(self):
    self.P = shrink_morph_py.parameterization(self.V, F, self.lambda1, self.lambda2, 0, self.n_iter, self.lim)
    scale = self.width / (np.max(self.P) - np.min(self.P))
    self.P *= scale
    self.V *= scale

    ps.register_surface_mesh("Parameterization", self.P, self.F, material="flat")

    sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
    ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
    ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
    ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')

    self.display_buildplate()

    ps.set_user_callback(self.callback_param)
    ps.show()

  def callback_param(self):
    # global lambda1, lambda2, wD, E1, thickness, deltaLambda, n_layers, lim, n_iter, width, P, self.V, self.F, printer_profile, with_smoothing
    gui.PushItemWidth(200)
    changed = gui.BeginCombo("Select printer", self.printer_profile)
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val, self.printer_profile==val)
        if selected:
          self.printer_profile = val
          self.printer = togcode.Printer(self.printer_profile)
          build_vert = np.array([[-1,-1,-0.1], [1, -1, -0.1], [1, 1, -0.1], [-1, 1, -0.1]])
          build_vert[:, 0] *= self.printer.bed_size[0] / 2
          build_vert[:, 1] *= self.printer.bed_size[1] / 2
          ps.get_surface_mesh("Buildplate").update_vertex_positions(build_vert)
      gui.EndCombo()
    gui.PopItemWidth()

    gui.PushItemWidth(100)
    changed, self.with_smoothing = gui.Checkbox("With smoothing", self.with_smoothing) 
    if changed:
      self.P = shrink_morph_py.reparameterization(self.V, self.P, self.F, self.lambda1, self.lambda2, self.wD if self.with_smoothing else 0, self.n_iter, self.lim)
      ps.get_surface_mesh("Parameterization").update_vertex_positions(self.P)

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')
    changed, self.width = gui.DragFloat("Width", self.width, 1, 1, 500, "%.0f")
    if changed and self.width > 0:
      scale = self.width / (np.max(self.P) - np.min(self.P))
      self.P *= scale
      self.V *= scale
      ps.get_surface_mesh("Parameterization").update_vertex_positions(self.P)

    if gui.Button("Increase mesh resolution"):
      self.V, self.P, self.F, _ = shrink_morph_py.subdivide(self.V, self.P, self.F, np.array([]))
      self.P = shrink_morph_py.reparameterization(self.V, self.P, self.F, self.lambda1, self.lambda2, self.wD if self.with_smoothing else 0, self.n_iter, self.lim)
      ps.register_surface_mesh("Parameterization", self.P, self.F, material="flat")

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')

    if gui.TreeNode("Advanced"):
      changed, self.lambda1 = gui.InputFloat("self.lambda1", self.lambda1, 0, 0, "%.1e")
      changed, self.lambda2 = gui.InputFloat("self.lambda2", self.lambda2, 0, 0, "%.1e")
      changed, self.thickness = gui.InputFloat("self.Thickness", self.thickness, 0, 0, "%.1e")
      changed, self.n_iter = gui.InputInt("Iterations", self.n_iter, step=1)
      changed, self.lim = gui.InputFloat("self.Limit", self.lim, 0, 0, "%.1e")

    if gui.Button("Next"):
      self.leave = False
      ps.unshow()
    
    # TODO pass the new variables: rect_width and rect_length should be passed to the generate_calibration function
    _, self.num_rectangles = gui.DragFloat("Number of Rectangles", self.num_rectangles, 1, 1, 10, "%.0f")
    _, self.num_layers = gui.DragFloat("Number of Layers", self.num_layers, 1, 1, 20, "%.0f")
    _, self.layer_height = gui.DragFloat("Layer Height", self.layer_height, 0.01, 0, 50, "%.2f")
    _, self.printer.bed_temp = gui.DragFloat("Bed temperature", self.printer.bed_temp, 1, 20, 60, "%.0f")
    _, self.printer.extruder_temp = gui.DragFloat("Nozzle temperature", self.printer.extruder_temp, 1, 180, 230, "%.0f")
    _, self.printer.print_speed = gui.DragFloat("Printing speed (mm/s)", self.printer.print_speed, 1, 10, 100, "%.0f")
    _, rect_width = gui.DragFloat("Rectangle width (mm)", 20, 1, 1, 20, "%.0f")
    _, rect_length = gui.DragFloat("Rectangle length (mm)", 80, 1, 1, 20, "%.0f")
    
    if gui.Button("Generate Calibration G-code"):
      self.generate_calibration(self.num_rectangles, self.num_layers, self.layer_height)

  def show(self, V, F):
    ps.set_give_focus_on_show(True)
    ps.init()
    ps.set_up_dir("neg_y_up")

    ps.load_color_map("twilight", "data/twilight_colormap.png");
    ps.set_ground_plane_mode("shadow_only")

    self.V = V
    self.F = F

    self.param_screen()

    if self.leave:
      return

    self.optim_screen()

    if self.leave:
      return
    
    self.traj_screen()


  # Directions optimization
  def optim_screen(self):
    self.leave = True
    ps.remove_all_structures()
    ps.register_surface_mesh("Input mesh", self.V, self.F)

    self.targetV = self.V.copy()
    self.theta2 = np.zeros(self.V.shape[0])

    ps.set_user_callback(self.callback_optim)
    ps.reset_camera_to_home_view()
    ps.show()

  resolutions = ["Low", "Medium", "High"]
  resolution = resolutions[0]
  def callback_optim(self):
      
    gui.PushItemWidth(100)
    changed, self.width = gui.InputFloat("Width", self.width, 0, 0, "%.0f")
    if changed and self.width > 0:
      scale = self.width / (np.max(self.P) - np.min(self.P))
      self.P *= scale
      self.V *= scale
      self.targetV *= scale
      ps.get_surface_mesh("Input mesh").update_vertex_positions(self.targetV)

    if gui.Button("Simulation"):
      shrink_morph_py.simulation(self.V, self.P, self.F, self.theta2, self.E1, self.lambda1, self.lambda2, self.deltaLambda, self.thickness, self.width, self.n_iter, self.lim)
      ps.get_surface_mesh("Input mesh").set_transparency(0.5)
      ps.register_surface_mesh("Simulation", self.V, self.F)

    if gui.Button("Directions optimization"):
      self.theta2 = shrink_morph_py.directions_optimization(self.V, self.targetV, self.P, self.F, self.E1, self.lambda1, self.lambda2, self.deltaLambda, self.thickness, self.width,
                                          self.n_iter, self.lim, self.wM, self.wL)
      ps.get_surface_mesh("Input mesh").set_transparency(0.5)
      if ps.has_surface_mesh("Simulation"):
        ps.get_surface_mesh("Simulation").update_vertex_positions(self.V)
      else:
        ps.register_surface_mesh("Simulation", self.V, self.F)
      ps.get_surface_mesh("Simulation").add_scalar_quantity("theta2", self.theta2)
      ps.get_surface_mesh("Simulation").add_scalar_quantity("theta1", self.angles, defined_on='faces', vminmax=(-np.pi/2, np.pi/2), cmap='twilight')

    changed = gui.BeginCombo("Trajectory resolution", self.resolution)
    if changed:
      for val in self.resolutions:
        _, selected = gui.Selectable(val, self.resolution==val)
        if selected:
          self.resolution = val
      gui.EndCombo()

    if gui.Button("Generate trajectories"):
      self.leave = False
      ps.unshow()

  def traj_screen(self):
    # Trajectories & G-code generation
    if self.resolution == "Low":
      target_edge_length = 1
    elif self.resolution == "Medium":
      target_edge_length = 0.5
    elif self.resolution == "High":
      target_edge_length = 0.2
    self.V, self.P, self.F, self.theta2 = shrink_morph_py.subdivide(self.V, self.P, self.F, self.theta2, target_edge_length)
    self.trajectories = shrink_morph_py.generate_trajectories(self.V, self.P, self.F, self.theta2, self.printer.layer_height, self.printer.nozzle_width, self.n_layers)
    ps.remove_all_structures();
  
    self.display_trajectories()
    self.display_buildplate()

    ps.set_user_callback(self.callback_traj)
    ps.reset_camera_to_home_view()
    ps.show()

  def display_trajectories(self):
    nodes = self.trajectories[0]
    edges = np.empty([nodes.shape[0] - 1, 2])
    edges[:, 0] = np.arange(nodes.shape[0] - 1)
    edges[:, 1] = np.arange(1, nodes.shape[0])
    colors = np.zeros(nodes.shape[0])
    path_id = 1

    h = nodes[0,2]
    k = 1
    for traj in self.trajectories:
      if traj[0, 2] > h:
        ps_traj = ps.register_curve_network("Layer " + str(k), nodes, edges, enabled=False)
        ps_traj.set_radius(self.printer.nozzle_width / 2, relative=False)
        ps_traj.add_scalar_quantity("Ordering", colors, enabled=True, cmap="blues")
        nodes = traj
        edges = np.empty([nodes.shape[0] - 1, 2])
        edges[:, 0] = np.arange(nodes.shape[0] - 1)
        edges[:, 1] = np.arange(1, nodes.shape[0])
        colors = np.zeros(nodes.shape[0])

        h = nodes[0,2]
        path_id = 1
        k += 1

      new_edges = np.empty([traj.shape[0] - 1, 2])
      new_edges[:, 0] = np.arange(nodes.shape[0], nodes.shape[0] + traj.shape[0] - 1)
      new_edges[:, 1] = np.arange(nodes.shape[0] + 1, nodes.shape[0] + traj.shape[0])
      
      nodes = np.vstack((nodes, traj))
      edges = np.vstack((edges, new_edges))
      colors = np.hstack((colors, path_id * np.ones(traj.shape[0])))
      path_id += 1

    ps_traj = ps.register_curve_network("Layer " + str(k), nodes, edges, enabled=True, radius=self.printer.nozzle_width / 2)
    ps_traj.add_scalar_quantity("Ordering", colors, enabled=True, cmap="blues")
    ps_traj.set_radius(self.printer.nozzle_width / 2, relative=False)


  layer_id = 1
  def callback_traj(self):
    # global self.layer_id, self.V, self.P, self.F, self.theta2, self.printer_profile, self.trajectories, self.printer
    gui.PushItemWidth(200)
    changed = gui.BeginCombo("Select self.printer", self.printer_profile)
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val, self.printer_profile==val)
        if selected:
          self.printer_profile = val
          self.printer = togcode.Printer(self.printer_profile)
          self.display_buildplate()
      gui.EndCombo()
    gui.PopItemWidth()

    if gui.Button("Increase mesh resolution and reload trajectories"):
      self.V, self.P, self.F, self.theta2 = shrink_morph_py.subdivide(self.V, self.P, self.F, self.theta2)
      ps.remove_all_structures();
      self.trajectories = shrink_morph_py.generate_trajectories(self.V, self.P, self.F, self.theta2, self.printer.layer_height, self.printer.nozzle_width, self.n_layers)
      self.display_trajectories(self.trajectories)
      self.display_buildplate()

    gui.PushItemWidth(200)
    changed, self.layer_id = gui.SliderInt("Layer", self.layer_id, 1, self.n_layers)
    gui.PopItemWidth()
    if changed:
      for i in range(1, self.n_layers + 1):
        if i == self.layer_id:
          ps.get_curve_network("Layer " + str(i)).set_enabled(True)
        else:
          ps.get_curve_network("Layer " + str(i)).set_enabled(False)
    if gui.Button("Export to g-code"):
      filename = filedialog.asksaveasfilename(defaultextension='.gcode')
      self.printer.to_gcode(self.trajectories, filename)
  
  def read_trajectories(self, filename):
    print("reading file " + filename)

    with open(filename, "r") as file:
        paths = []
        line = file.readline()
        while line:
            num_vertices = int(line)
            path = []
            for _ in range(num_vertices):
                line = file.readline()
                vertex = line.split()
                path.append([float(x) for x in vertex])
            paths.append(np.array(path))
            line = file.readline()
    return paths

  def generate_calibration(self, num_rectangles, num_layers, layer_height):
    length = 80
    width = 20
    nb_layers = int(num_layers)
    layer_height = float(layer_height)
    jump_y = width+10
    jump_x = 0
    shift_y = 2 * jump_y
    shift_x = 0

    paths = []
    with open("sample.path", "w") as file:
        #layer_height = 0.08
        for j in range(int(num_rectangles)):
            nb_layers = 10-j
            z = 0
            for i in range(nb_layers):
                z += layer_height + i / (nb_layers - 1) * 2 * (0.8 / nb_layers - 0.08)
                y = -width/2 -j*jump_y + shift_y
                x = -length/2 + shift_x
                while(y < width/2 - j*jump_y + shift_y):
                    path = []
                    if x < 0 + shift_x:
                        file.write("2\n")
                        path.append([x, y, z])
                        x = length/2 + shift_x
                        path.append([x, y, z])
                    else:
                        file.write("2\n")
                        path.append([x, y, z])
                        x = -length/2 + shift_x
                        path.append([x, y, z])
                    y += 0.4
                    paths.append(np.array(path))
                print(f"{z:.4f}") # for debug purposes
    self.printer.to_gcode(paths, "data/calibration.gcode")

main = ShrinkMorph()
main.show(V, F)