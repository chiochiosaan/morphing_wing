import numpy as np

def analysis_of_a_deflection_by_coords(airfoil_points, P, alpha_deg, E, I):
    airfoil_points = np.array(airfoil_points)
    force_application_point_x, force_application_point_y = airfoil_points[-1]    
    #Разложение силы на проекции
    alpha = np.radians(alpha_deg)
    P_x = P * np.cos(alpha)
    P_y = P * np.sin(alpha)
    
    EI = E * I
    delta_x = 0.0
    delta_y = 0.0
    M_p_array = []  # Массив для эпюры моментов
    
    for i in range(len(airfoil_points) - 1):
        x1, y1 = airfoil_points[i]
        x2, y2 = airfoil_points[i+1]
        
        x_mid = (x1 + x2) / 2.0
        y_mid = (y1 + y2) / 2.0
        
        current_length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        arm_x = force_application_point_x - x_mid
        arm_y = force_application_point_y - y_mid
        
        M_p = arm_x * P_y - arm_y * P_x
        M_p_array.append(M_p)
        
        M_1x = arm_x * 0 - arm_y * 1
        
        M_1y = arm_x * 1 - arm_y * 0
        
        # Интеграл Мора
        delta_x += (M_p * M_1x / EI) * current_length
        delta_y += (M_p * M_1y / EI) * current_length
        
    delta_total = np.sqrt(delta_x**2 + delta_y**2)
    return delta_x, delta_y, delta_total, M_p_array


# Константы материала
P = 1.0    #Н
E = 2000   #МПа
I = 42.67  # мм^4

Radius = 40.0
n_points = 100
phi = np.linspace(-np.pi/2, np.pi/2, n_points)

x_semicircle = Radius * np.cos(phi)
y_semicircle = Radius * np.sin(phi)

airfoil_points = np.column_stack((x_semicircle, y_semicircle))


# Продольное (угол 0)
dx_0, dy_0, dtotal_0, bending_moment_0 = analysis_of_a_deflection_by_coords(airfoil_points, P, 0, E, I)
print("продольное нагружение (alpha=0°):")
print(f"   Смещение по X: {dx_0:+.3f} мм")
print(f"   Смещение по Y: {dy_0:+.3f} мм")
print(f"   Полное перемещение: {dtotal_0:.3f} мм\n")

# Поперечное (угол 90)
dx_90, dy_90, dtotal_90, bending_moment_90 = analysis_of_a_deflection_by_coords(airfoil_points, P, 90, E, I)
print("поперечное нагружение (alpha=90°):")
print(f"   Смещение по X: {dx_90:+.3f} мм")
print(f"   Смещение по Y: {dy_90:+.3f} мм")
print(f"   Полное перемещение: {dtotal_90:.3f} мм\n")
