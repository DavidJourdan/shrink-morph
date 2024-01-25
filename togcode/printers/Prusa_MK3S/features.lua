-- Original Prusa MK3S
-- 01/05/0219

-- Build Area dimensions
bed_size_x_mm = 250
bed_size_y_mm = 210
bed_size_z_mm = 210

-- Printer Extruder
extruder_count = 1 -- number of extruders. Change this value if you want to use the virtual extruders feature for "simple multi-material" (using firmware's filament swap)

nozzle_diameter_mm = 0.4 -- 0.25, 0.4, 0.6
filament_diameter_mm = 1.75
filament_linear_adv_factor = 0 -- default

-- Retraction Settings
filament_priming_mm = 0.8 -- min 0.5 - max 2
priming_mm_per_sec = 35
retract_mm_per_sec = 35

enable_z_lift = true
z_lift_mm = 0.6

-- Layer height
z_layer_height_mm = 0.2
z_layer_height_mm_min = nozzle_diameter_mm * 0.10
z_layer_height_mm_max = nozzle_diameter_mm * 0.90

-- Printing temperatures
extruder_temp_degree_c = 210
extruder_temp_degree_c_min = 150
extruder_temp_degree_c_max = 270

bed_temp_degree_c = 55
bed_temp_degree_c_min = 0
bed_temp_degree_c_max = 120

-- Printing speeds
print_speed_mm_per_sec = 60
print_speed_mm_per_sec_min = 5
print_speed_mm_per_sec_max = 120

perimeter_print_speed_mm_per_sec = 45
perimeter_print_speed_mm_per_sec_min = 5
perimeter_print_speed_mm_per_sec_max = 80

cover_print_speed_mm_per_sec = 45
cover_print_speed_mm_per_sec_min = 5
cover_print_speed_mm_per_sec_max = 80

first_layer_print_speed_mm_per_sec = 20
first_layer_print_speed_mm_per_sec_min = 1
first_layer_print_speed_mm_per_sec_max = 50

travel_speed_mm_per_sec = 180

-- Acceleration settings
x_max_speed = 200 -- mm/s
y_max_speed = 200 -- mm/s
z_max_speed = 12 -- mm/s
e_max_speed = 120 -- mm/s

x_max_acc = 1000 -- mm/s²
y_max_acc = 1000 -- mm/s²
z_max_acc = 1000 -- mm/s²
e_max_acc = 5000 -- mm/s²
ex_max_acc = 1250 -- mm/s²
e_prime_max_acc = 1250 -- mm/s²

perimeter_acc = 800 -- mm/s²
infill_acc = 1250 -- mm/s²
default_acc = 1000 -- mm/s²

x_max_jerk = 8.00 -- mm/s
y_max_jerk = 8.00 -- mm/s
z_max_jerk = 0.40 -- mm/s
e_max_jerk = 1.50 -- mm/s