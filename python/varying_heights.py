length = 80
width = 20
nb_layers = 10
layer_height = 0.08
jump_y = width+10
jump_x = 0
shift_y = 2 * jump_y
shift_x = 0

with open("sample.path", "w") as file:
    nb_layers = 9
    layer_height = 0.08
    for j in range(5):
        z = 0
        for i in range(nb_layers):
            z += layer_height + i / (nb_layers - 1) * 2 * (0.8 / nb_layers - 0.08)
            y = -width/2 -j*jump_y + shift_y
            x = -length/2 + shift_x
            while(y < width/2 - j*jump_y + shift_y):
                if x < 0 + shift_x:
                    file.write("2\n")
                    file.write(f"{x:.2f} {y:.2f} {z:.4f}\n")
                    x = length/2 + shift_x
                    file.write(f"{x:.2f} {y:.2f} {z:.4f}\n")
                else:
                    file.write("2\n")
                    file.write(f"{x:.2f} {y:.2f} {z:.4f}\n")
                    x = -length/2 + shift_x
                    file.write(f"{x:.2f} {y:.2f} {z:.4f}\n")
                # y += 0.3 + i / (nb_layers - 1) * 0.5
                y += 0.4
            print(f"{z:.4f}")
