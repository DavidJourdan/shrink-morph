-- Build Area dimensions
bed_size_x_mm = 310
bed_size_y_mm = 310
bed_size_z_mm = 400

-- Printer Extruder
extruder_count = 1 -- number of extruders. Change this value if you want to use the virtual extruders feature for "simple multi-material" (using firmware's filament swap)
nozzle_diameter_mm = 0.4
filament_diameter_mm = 1.75

-- Retraction Settings
-- between 0.5mm and 0.8mm of retract/prime for direct-drive setup, between 3mm and 6mm for bowden (stock) setup
filament_priming_mm = 0.4 
priming_mm_per_sec = 45
retract_mm_per_sec = 45

-- Layer height
z_layer_height_mm = 0.2
z_layer_height_mm_min = nozzle_diameter_mm * 0.10
z_layer_height_mm_max = nozzle_diameter_mm * 0.90

-- Printing temperatures
extruder_temp_degree_c = 210
extruder_temp_degree_c_min = 150
extruder_temp_degree_c_max = 270

bed_temp_degree_c = 50
bed_temp_degree_c_min = 0
bed_temp_degree_c_max = 120

-- Printing speeds
print_speed_mm_per_sec = 60
print_speed_mm_per_sec_min = 5
print_speed_mm_per_sec_max = 200

perimeter_print_speed_mm_per_sec = 45
perimeter_print_speed_mm_per_sec_min = 5
perimeter_print_speed_mm_per_sec_max = 150

cover_print_speed_mm_per_sec = 45
cover_print_speed_mm_per_sec_min = 5
cover_print_speed_mm_per_sec_max = 200

first_layer_print_speed_mm_per_sec = 10
first_layer_print_speed_mm_per_sec_min = 5
first_layer_print_speed_mm_per_sec_max = 60

travel_speed_mm_per_sec = 180
travel_speed_mm_per_sec_min = 5
travel_speed_mm_per_sec_max = 200

-- Misc default settings
add_brim = true
brim_distance_to_print_mm = 2.0
brim_num_contours = 3

travel_straight = true
enable_z_lift = true
z_lift_mm = 0.4

-- default filament infos (when using "custom" profile)
name_en = "PLA"
filament_density = 1.25 --g/cm3 PLA