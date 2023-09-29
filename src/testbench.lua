-- custom settings before calling header()
z_layer_height_mm = 0.08
print_speed = 30
first_layer_speed = 30
flow_multiplier = 0.93
bed_temp_degree_c = 55
extruder_temp_degree_c_0 = 190
nozzle_diameter_mm_0 = 0.4
travel_max_length_without_retract = 2
min_z_lift_length = 5
travel_speed = 120
z_lift = 0.4
z_lift_speed = 10
filament_priming = 0.4
retract_speed = 25

number_of_layers = 10
previous = {x=0., y=0., z=0.}

function read_trajectories(filename)
    print('reading file ' .. filename)
    local file = io.open(filename, "r")
    io.input(file)
    local paths = {}

    local line = io.read()
    while line do
        local num_vertices = string.match(line,'%d+')
        local path = {}
        for j=1,num_vertices do
            local n = 1
            local line = io.read()
            for v in string.gmatch(line,'%S+') do
              if n == 1     then vx = tonumber(v)
              elseif n == 2 then vy = tonumber(v)
              elseif n == 3 then vz = tonumber(v)
              end
              n = n + 1
            end
            vx = vx + bed_size_x_mm/2
            vy = vy + bed_size_y_mm/2
            table.insert(path, {x=vx, y=vy, z=vz})
        end

        table.insert(paths,path)
        line = io.read()
    end
    return paths
end

total_e = 0

function e(w,l,height) -- w: width, l: length
    height = height or z_layer_height_mm
    local  crsec = math.pi * filament_diameter_mm_0 * filament_diameter_mm_0 / 4.0
    local  v_mm3 = l * height * w
    total_e      = total_e + flow_multiplier * v_mm3 / crsec
    return total_e
end

function travel(from, to)
    local dist = math.sqrt(math.pow(from.x - to.x, 2) + math.pow(from.y - to.y, 2))
    local dir = {x=(to.x - from.x) / dist, y=(to.y - from.y) / dist}

    if dist > min_z_lift_length then

        if dist > travel_max_length_without_retract then
            set_feedrate(retract_speed*60)
            move_e(total_e - filament_priming)
        end

        z = math.max(from.z, to.z) + z_lift

        set_feedrate(z_lift_speed*60)
        move_xyz(from.x + z_lift * dir.x, from.y + z_lift * dir.y, z)
        
        set_feedrate(travel_speed*60)
        move_xyz(to.x - z_lift * dir.x, to.y - z_lift * dir.y, z)

        set_feedrate(z_lift_speed*60)
        move_xyz(to.x, to.y, to.z)

        if dist > travel_max_length_without_retract then
            set_feedrate(retract_speed*60)
            move_e(total_e)
        end
    else
        set_feedrate(travel_speed*60)
        move_xyz(to.x, to.y, to.z)
    end

    if to.z == z_layer_height_mm then
        set_feedrate(first_layer_speed*60)
    else
        set_feedrate(print_speed*60)
    end
end



function add_trajectories(paths)
    for k,path in pairs(paths) do
        travel(previous, path[1])

        for i,v in pairs(path) do
            if i > 1 then
                local l  = math.sqrt(math.pow(v.x - path[i-1].x, 2) + math.pow(v.y - path[i-1].y, 2))
                move_xyze(v.x, v.y, v.z, e(nozzle_diameter_mm_0, l))
            end
        end
        previous = path[#path]
    end
end

function print_square(width, height, size)
    X = width / 2
    Y = height / 2

    travel(previous, {x=-X + bed_size_x_mm/2, y=-Y + bed_size_y_mm/2, z=z_layer_height_mm})

    for i=1,size do
        move_xyz(-X + bed_size_x_mm/2, -Y + bed_size_y_mm/2, z_layer_height_mm)
        move_xyze(X + bed_size_x_mm/2, -Y + bed_size_y_mm/2, z_layer_height_mm, e(nozzle_diameter_mm_0, 2 * X))
        move_xyze(X + bed_size_x_mm/2, Y + bed_size_y_mm/2, z_layer_height_mm, e(nozzle_diameter_mm_0, 2 * Y))
        move_xyze(-X + bed_size_x_mm/2, Y + bed_size_y_mm/2, z_layer_height_mm, e(nozzle_diameter_mm_0, 2 * X))
        move_xyze(-X + bed_size_x_mm/2, -Y + bed_size_y_mm/2, z_layer_height_mm, e(nozzle_diameter_mm_0, 2 * Y))

        X = X - nozzle_diameter_mm_0
        Y = Y - nozzle_diameter_mm_0
    end

    previous = {x=-X + bed_size_x_mm/2, y=-Y + bed_size_y_mm/2, z=z_layer_height_mm}
end

function print_sample(width, height, centerX, centerY, layer_height, n)
    n = n or number_of_layers
    layer_height = layer_height or z_layer_height_mm
    x = -width/2 + centerX
    y = -height/2 + centerY
    travel(previous, {x=x, y=y+3*nozzle_diameter_mm_0/2, z=z_layer_height_mm})
    for layer=1,n do
        if layer > 1 then
            set_feedrate(print_speed*60.0)
        end
    
        z = layer*layer_height

        if layer%2 == 0 then
            y = -height/2 + centerY
        else
            y = -height/2 + nozzle_diameter_mm_0/2 + centerY
        end
    
        while y < height/2 + centerY do
            travel({x=x, y=y, z=z}, {x=x, y=y+nozzle_diameter_mm_0, z=z})
            y = y + nozzle_diameter_mm_0
    
            if x > centerX then
                move_xyze(x - width, y, z, e(nozzle_diameter_mm_0, width, layer_height))
                x = x - width
            else
                move_xyze(x + width, y, z, e(nozzle_diameter_mm_0, width, layer_height))
                x = x + width
            end
        end
    
        if layer < n then
            if layer%2 == 0 then
                travel({x=x, y=y, z=z}, {x=x, y=-height/2 + centerY + 3*nozzle_diameter_mm_0/2, z=z+layer_height})
            else
                travel({x=x, y=y, z=z}, {x=x, y=-height/2 + centerY + nozzle_diameter_mm_0, z=z+layer_height})
            end
        end
    end
    previous = {x=x, y=height/2 + centerY, z=z}
end

header()
select_extruder(0)

print('nozzle_diameter_mm_0   = ' .. nozzle_diameter_mm_0)
print('filament_diameter_mm_0 = ' .. filament_diameter_mm_0)
print('z_layer_height_mm      = ' .. z_layer_height_mm)

set_feedrate(first_layer_speed*60.0)

print_square(170, 110, 5)
add_trajectories(read_trajectories("/home/djourdan/code/shrink-morph/data/saddle.path"))

footer()

print('total length      = ' .. total_e / (math.pi * filament_diameter_mm_0 * filament_diameter_mm_0 / 4.0))
