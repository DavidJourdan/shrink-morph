import argparse
import os.path 
import numpy as np

class PrinterParameters:
    def __init__(self):
        self.printer_profile = None
        self.nozzle_width = None
        self.print_speed = None
        self.retract_speed = None
        self.first_layer_speed = None
        self.travel_speed = None
        self.flow_multiplier = None
        self.bed_temp = None
        self.extruder_temp = None
        self.layer_height = None
        self.filament_diameter = None
        self.filament_priming = None
        self.z_lift = None
        self.total_extrusion_length = None
        
        # constants
        self.z_lift_speed = 10
        self.travel_max_length_without_retract = 2
        self.linear_advance_factor = 30

    def get_print_feedrate(self) -> int:
        # mm/s to mm/min
        return int(self.print_speed * 60)

    def get_first_layer_feedrate(self) -> int:
        # mm/s to mm/min
        return int(self.first_layer_speed * 60)

    def get_travel_feedrate(self) -> int:
        # mm/s to mm/min
        return int(self.travel_speed * 60)

    def get_retract_feedrate(self) -> int:
        # mm/s to mm/min
        return int(self.retract_speed * 60)

    def get_z_lift_feedrate(self) -> int:
        # mm/s to mm/min
        return int(self.retract_speed * 60)

def compute_extrusion_length(width: float,
                             local_length: float,
                             filament_diameter: float,
                             layer_height: float,
                             flow_multiplier: float) -> float:
    """Compute the length of the filament to extrude.

    Notes
    -----
    https://3dprinting.stackexchange.com/questions/6289/how-is-the-e-argument-calculated-for-a-given-g1-command
    """
    crsec = np.pi * filament_diameter**2 / 4.0
    v_mm3 = local_length * layer_height * width
    return flow_multiplier * v_mm3 / crsec

def update_extrusion_length(local_length: float,
                            param: PrinterParameters):
    param.total_extrusion_length += compute_extrusion_length(
        param.nozzle_width, local_length, param.filament_diameter, param.layer_height, param.flow_multiplier)


def features_get_bed_size_xy(features_str: str) -> np.ndarray:

    bed_circular = features_str.find("bed_circular  = true") != -1
    bed_circular = bed_circular or features_str.find(
        "bed_circular = true") != -1
    if not bed_circular:
        # The bed is not circular

        features_bed_size_x_i = features_str.find("bed_size_x_mm")
        bed_size_x = features_str[features_bed_size_x_i:]
        bed_size_x_i_end = bed_size_x.find('\n')
        bed_size_x = bed_size_x[:bed_size_x_i_end]
        bed_size_x = float(bed_size_x.split()[2])

        features_bed_size_y_i = features_str.find("bed_size_y_mm")
        bed_size_y = features_str[features_bed_size_y_i:]
        bed_size_y_i_end = bed_size_y.find('\n')
        bed_size_y = bed_size_y[:bed_size_y_i_end]
        bed_size_y = float(bed_size_y.split()[2])

        return np.array([bed_size_x, bed_size_y])
    else:

        bed_radius = features_str.find("bed_radius")
        bed_radius = features_str[bed_radius:]
        bed_radius_i_end = bed_radius.find('\n')
        bed_radius = bed_radius[:bed_radius_i_end]
        bed_radius = float(bed_radius.split()[2])

        return np.array([bed_radius, bed_radius]) * 2.


def features_origin_is_bed_center(features_str: str) -> np.ndarray:
    origin_is_bed_center = False
    index = features_str.find("bed_origin_x = bed_size_x_mm / 2.0")
    origin_is_bed_center = index != -1 or origin_is_bed_center
    index = features_str.find("bed_origin_x = bed_size_x_mm/2")
    origin_is_bed_center = index != -1 or origin_is_bed_center

    return origin_is_bed_center


def travel(point_start: np.ndarray, point_end: np.ndarray, printer_param: PrinterParameters) -> str:
    dist = np.linalg.norm(point_end - point_start)

    travel_str = ""
    if dist < printer_param.travel_max_length_without_retract:
        travel_str += ";travel\n"
        travel_str += f"G0 F{printer_param.get_travel_feedrate()} X{point_end[0]:.3f} Y{point_end[1]:.3f} Z{point_end[2]:.2f}\n"
    else:
        point_z_lifted = max(point_end[2], point_start[2]) + printer_param.z_lift

        direction = (point_end - point_start) / dist
        intermediate1 = point_start + printer_param.z_lift * direction
        intermediate2 = point_end - printer_param.z_lift * direction

        # Retract
        travel_str += ";retract\n"
        printer_param.total_extrusion_length -= printer_param.filament_priming
        travel_str += f"G1 F{printer_param.get_retract_feedrate()} E{printer_param.total_extrusion_length:.5f}\n"
        # Travel
        travel_str += ";travel\n"
        travel_str += f"G0 F{printer_param.get_z_lift_feedrate()} X{intermediate1[0]:.3f} Y{intermediate1[1]:.3f} Z{point_z_lifted:.2f}\n"
        travel_str += f"G0 F{printer_param.get_travel_feedrate()} X{intermediate2[0]:.3f} Y{intermediate2[1]:.3f} Z{point_z_lifted:.2f}\n"
        travel_str += f"G0 F{printer_param.get_z_lift_feedrate()} X{point_end[0]:.3f} Y{point_end[1]:.3f} Z{point_end[2]:.2f}\n"
        # Prime
        travel_str += ";prime\n"
        printer_param.total_extrusion_length += printer_param.filament_priming
        travel_str += f"G1 F{printer_param.get_retract_feedrate()} E{printer_param.total_extrusion_length:.5f}\n"
    return travel_str


def travel_to(point_end: np.ndarray, 
              printer_param: PrinterParameters, 
              prime: bool = True) -> str:

    # Retract
    travel_str = ";retract\n"
    printer_param.total_extrusion_length -= printer_param.filament_priming
    travel_str += f"G1 F{printer_param.get_retract_feedrate()} E{printer_param.total_extrusion_length:.5f}\n"
    # Travel
    travel_str += ";travel\n"
    travel_str += f"G0 F{printer_param.get_print_feedrate()} X{point_end[0]:.3f} Y{point_end[1]:.3f} Z{point_end[2]:.2f}\n"
    if prime:
        # Prime
        travel_str += ";prime\n"
        printer_param.total_extrusion_length += printer_param.filament_priming
        travel_str += f"G1 F{printer_param.get_retract_feedrate()} E{printer_param.total_extrusion_length:.5f}\n"
    return travel_str

def read_trajectories(filename):
    print('reading file ' + filename)

    with open(filename, 'r') as file:
        paths = []
        line = file.readline()
        while line:
            num_vertices = int(line)
            path = []
            for i in range(num_vertices):
                line = file.readline()
                vertex = line.split()
                path.append([float(x) for x in vertex])
            paths.append(np.array(path))
            line = file.readline()
    return paths

def run():

    parser = argparse.ArgumentParser(
        description='Convert the cycle generated by fill_2d_shape.py to machine instructions. Input: The JSON file used by fill_2d_shape.py. Output: The G-code to print the generated cycle.'
    )
    parser.add_argument(
        "input_filename", help="File containing the trajectory data")
    parser.add_argument(
        "printer_profile", help="The printer profile. Only 'CR10S_Pro' was tested. In theory, the name of the folders at `src/ext/iceslprinters/fff` are valid inputs.")
    parser.add_argument(
        "-o",
        "--output",
        help="Name of the output gcode file",
        default="output.gcode")
    parser.add_argument(
        "-nw",
        "--nozzle_width",
        help="The nozzle diameter in millimeter. Default: 0.4mm.",
        type=float,
        default=0.4)
    parser.add_argument(
        "-s",
        "--print_speed",
        help="The speed of the moving head, in mm/s. Default: 30mm/s.",
        type=float,
        default=30)
    parser.add_argument(
        "-rs",
        "--retract_speed",
        help="Speed of retraction. Default: 45mm/s.",
        type=float,
        default=45)
    parser.add_argument(
        "-fs",
        "--first_layer_speed",
        help="The speed of the moving head for the first layer, in mm/s. Default: 30mm/s.",
        type=float,
        default=30)
    parser.add_argument(
        "-ts",
        "--travel_speed",
        help="Speed of travels, in mm/s. Default: 120mm/s.",
        type=float,
        default=120)
    parser.add_argument(
        "-fm",
        "--flow_multiplier",
        help="The flow multiplier. Default: 1",
        type=float,
        default=1.0)
    parser.add_argument(
        "-bt",
        "--bed_temp",
        help="The bed temperature in degree Celsius. Default: 55.",
        type=int,
        default=55)
    parser.add_argument(
        "-et",
        "--extruder_temp",
        help="The extruder temperature in degree Celsius. Default: 190.",
        type=int,
        default=190)
    parser.add_argument(
        "-lh",
        "--layer_height",
        help="The layer height. Default: the layer height associated with the polyline. Default: 0.08mm",
        type=float,
        default=0.08)
    parser.add_argument(
        "-fd",
        "--filament_diameter",
        help="The filament diameter of the filament used by the printer. Default: 1.75 mm.",
        type=float,
        default=1.75)
    parser.add_argument(
        "-fp",
        "--filament_priming",
        help="Retraction setting. Between 0.4mm and 0.8mm of retract/prime for direct-drive setup, between 3mm and 6mm for bowden (stock) setup. Default: 0.4 mm.",
        type=float,
        default=0.4)
    parser.add_argument(
        "-zl",
        "--z_lift",
        help="Distance to move the printhead up (or the build plate down), after each retraction, right before a travel move takes place. Default: 0.4",
        type=float,
        default=0.4)
    parser.add_argument(
        "--ncols",
        help="Object count along the x-axis. Default: 1.",
        type=int,
        default=1)
    parser.add_argument(
        "--nrows",
        help="Object count along the y-axis. Default: 1.",
        type=int,
        default=1)

    args = parser.parse_args()

    input_filename = args.input_filename
    printer_profile = args.printer_profile
    gcode_filename = args.output

    printer_param = PrinterParameters()

    printer_param.nozzle_width = args.nozzle_width
    printer_param.print_speed = args.print_speed
    printer_param.retract_speed = args.retract_speed
    printer_param.first_layer_speed = args.first_layer_speed
    printer_param.travel_speed = args.travel_speed
    printer_param.flow_multiplier = args.flow_multiplier
    printer_param.bed_temp = args.bed_temp
    printer_param.extruder_temp = args.extruder_temp
    printer_param.layer_height = args.layer_height
    printer_param.filament_diameter = args.filament_diameter
    printer_param.filament_priming = args.filament_priming
    printer_param.z_lift = args.z_lift

    trajectories = read_trajectories(input_filename)

    object_max_xy = np.array([0, 0, 0])
    object_min_xy = np.array([0, 0, 0])
    for trajectory in trajectories:      
        object_max_xy = np.maximum(np.amax(trajectory, axis=0), object_max_xy)
        object_min_xy = np.minimum(np.amin(trajectory, axis=0), object_min_xy)

    object_width_xy = object_max_xy[:2] - object_min_xy[:2]
    object_width_xy_plus_margin = object_width_xy + 5

    # add rectangle around object
    rect_trajectories = []
    for i in range(5):
        rect_trajectory = []
        X = object_width_xy_plus_margin[0] / 2 - printer_param.nozzle_width * i
        Y = object_width_xy_plus_margin[1] / 2 - printer_param.nozzle_width * i
        rect_trajectory.append([-X, -Y, printer_param.layer_height])
        rect_trajectory.append([X, -Y, printer_param.layer_height])
        rect_trajectory.append([X, Y, printer_param.layer_height])
        rect_trajectory.append([-X, Y, printer_param.layer_height])
        rect_trajectory.append([-X, -Y, printer_param.layer_height])
        rect_trajectories.append(np.array(rect_trajectory))
    trajectories = rect_trajectories + trajectories

    # Get the printer profile header and footer
    source_folder = os.path.dirname(__file__)
    printer_profile_path_head = f"{source_folder}/icesl-printers/fff/{printer_profile}/"
    printer_profile_header_path = printer_profile_path_head + "header.gcode"
    printer_profile_footer_path = printer_profile_path_head + "footer.gcode"
    printer_profile_features_path = printer_profile_path_head + "features.lua"
    printer_profile_printer_path = printer_profile_path_head + "printer.lua"

    # File to string
    with open(printer_profile_header_path, 'r') as f:
        header_str = f.read()
    with open(printer_profile_footer_path, 'r') as f:
        footer_str = f.read()
    with open(printer_profile_features_path, 'r') as f:
        features_str = f.read()
    with open(printer_profile_printer_path, 'r') as f:
        printer_str = f.read()

    # Extract the bed size from features string
    bed_size_xy = features_get_bed_size_xy(features_str)
    origin_is_bed_center = features_origin_is_bed_center(printer_str)

    if not origin_is_bed_center:
        for trajectory in trajectories:
            trajectory += np.array([bed_size_xy[0]/ 2, bed_size_xy[1] / 2, 0])

    # Put user defined parameters in the header
    header_str = header_str.replace(
        "<HBPTEMP>", str(int(printer_param.bed_temp)))
    header_str = header_str.replace(
        "<TOOLTEMP>", str(int(printer_param.extruder_temp)))
    header_str = header_str.replace("<BEDLVL>", "G0 F6200 X0 Y0")
    header_str = header_str.replace("<NOZZLE_DIAMETER>", str(printer_param.nozzle_width))
    header_str = header_str.replace("<ACCELERATIONS>", '''
M201 X1000 Y1000 Z1000 E5000 ; sets maximum accelerations, mm/sec^2
M203 X200 Y200 Z12 E120 ; sets maximum feedrates, mm/sec
M204 P1250 R1250 T1250 ; sets acceleration (P, T) and retract acceleration (R), mm/sec^2
M205 X8 Y8 Z0.4 E1.5 ; sets the jerk limits, mm/sec
M205 S0 T0 ; sets the minimum extruding and travel feed rate, mm/sec
    ''')
    header_str = header_str.replace("<FILAMENT>", str(printer_param.linear_advance_factor))

    max_dim = object_max_xy - object_min_xy
    header_str = header_str.replace("<MAX_X>", str(max_dim[0]))
    header_str = header_str.replace("<MAX_Y>", str(max_dim[1]))
    header_str = header_str.replace("<MAX_Z>", str(max_dim[2]))

    # Total extrusion length
    printer_param.total_extrusion_length = 0.0

    with open(gcode_filename, 'w', encoding="utf-8") as f:
        # Write header
        f.write(header_str)

        prev_point = None
        for trajectory in trajectories:
            if prev_point is None:
                f.write(travel_to(trajectory[0,:], printer_param))
            else:
                f.write(travel(prev_point, trajectory[0,:], printer_param))
            prev_point = trajectory[0,:]

            for point in trajectory:
                dist = np.linalg.norm(point - prev_point)
                update_extrusion_length(dist, printer_param)

                # point = point + world_to_object + object_to_bed[row_index, col_index]
                f.write(
                    f"G1 F{printer_param.get_print_feedrate()} X{point[0]:.3f} Y{point[1]:.3f} Z{point[2]:.2f} E{printer_param.total_extrusion_length:.5f}\n")
                prev_point = point

        # Write footer
        f.write(footer_str)


if __name__ == "__main__":
    run()
