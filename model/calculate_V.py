import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull

def calculate_V_S(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    data = [line.strip().split('\t') for line in lines]
    list_all = []
    for row in data:
        extracted_values = row[0].split()
        row_new = [float(value) if '.' in value else int(value) for value in extracted_values]

        list_all.append(row_new)

    x = [row[23] for row in list_all]
    y = [row[24] for row in list_all]
    z = [row[25] for row in list_all]
    points = np.array([x, y, z]).T
    
    hull = Delaunay(points)
    simplices = hull.simplices

    for s in simplices:
        s = np.append(s, s[0])  # Close the triangle
        #ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-")

    tetrahedra = points[simplices]
    volumes = np.abs(np.linalg.det(tetrahedra[:, 1:] - tetrahedra[:, 0][:, np.newaxis]))
    total_volume = np.sum(volumes) / 6
    V_changfang = (max(x)-min(x))*(max(y)-min(y))*(max(z)-min(z))
    V_tuoqiu = V_changfang*3.14159/6
    """
    print("Total volume:", total_volume)
    print("Bounding Box Volume", V_changfang)
    print("Approximate Ellipsoid Volume", V_tuoqiu)
    """


    # 将三维点组合成点云
    points_S = np.column_stack((x, y, z))

    # 计算凸壳
    hull = ConvexHull(points_S)

    # 绘制凸壳表面
    for simplex in hull.simplices:
        simplex = np.append(simplex, simplex[0])  # Close the triangle
        #ax.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], "r-")

    # 计算凸壳的表面积
    convex_hull_surface_area = hull.area
    a = (max(x)-min(x))/2
    b = (max(y)-min(y))/2
    c = (max(z)-min(z))/2
    S_tuoqiu = 4*3.1415926*(a*b + b*c + c*a)/3

    return total_volume, V_changfang, V_tuoqiu, convex_hull_surface_area, S_tuoqiu


