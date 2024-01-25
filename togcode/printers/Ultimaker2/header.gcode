;FLAVOR:UltiGCode
;TIME:' .. time_se
;MATERIAL:' .. to_mm_cube( filament_tot_length_mm[extruders[0]])
;MATERIAL2:0
;NOZZLE_DIAMETER:' .. round(nozzle_diameter_mm_0,2
M82
G92 E0
M107

M190 S<HBPTEMP> ; wait for bed temperature to be reached
M104 S<TOOLTEMP> ; set temperature
G28 ; home all axes
<BEDLVL>
BED_MESH_PROFILE LOAD="default"
M109 S<TOOLTEMP> ; wait for extruder temperature to be reached

G92 E0
G1 Z1.0 F3000 ; move z up little to prevent scratching of surface
G1 X0.1 Y20 Z0.3 F5000.0 ; move to start-line position
G1 X0.1 Y200.0 Z0.3 F1500.0 E15 ; draw 1st line
G1 X0.4 Y200.0 Z0.3 F5000.0 ; move to side a little
G1 X0.4 Y20 Z0.3 F1500.0 E30 ; draw 2nd line
G92 E0 ; reset extruder
; done purging extruder
