
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
"""
with open('C:/Users/zhr/Desktop/tomo.txt', 'r') as file:
    lines = file.readlines()

xs = []
ys = []
zs = []
for line in lines:
    elements = line.split()
    if elements:
        x_coordinate = float(elements[0])
        y_coordinate = float(elements[1])
        z_coordinate = float(elements[2])
        xs.append(x_coordinate)
        ys.append(y_coordinate)
        zs.append(z_coordinate)

x_bar = sum(xs)/len(xs)
y_bar = sum(ys)/len(ys)
z_bar = sum(zs)/len(zs)

x_bar = 259.2108423913043
y_bar = 308.7937771739129
z_bar = 199.6304347826087

x = [value - x_bar for value in xs]
y = [value - y_bar for value in ys]
z = [value - z_bar for value in zs]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s = [200]*len(xs), c='b', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()


"""
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
#from scipy.spatial.qhull import QhullError

current_file_path = os.getcwd()
folder_path = current_file_path+"/surfacedata"

with open(folder_path+'/surface117.tbl', 'r') as file:
    lines = file.readlines()

# 处理每一行数据，将其转换为列表
data = [line.strip().split('\t') for line in lines]
list_all = []
# 打印数据
for row in data:
    extracted_values = row[0].split()
    row_new = [float(value) if '.' in value else int(value) for value in extracted_values]

    list_all.append(row_new)

x = [row[23] for row in list_all]
y = [row[24] for row in list_all]
z = [row[25] for row in list_all]
x_bar = (max(x)+min(x))/2
y_bar = (max(y)+min(y))/2
z_bar = (max(z)+min(z))/2
xs = [value - x_bar for value in x]
ys = [value - y_bar for value in y]
zs = [value - z_bar for value in z]
x_new = [0.981*x + 0.191*y for x, y in zip(xs, ys)]
y_new = [-0.191*x + 0.981*y for x, y in zip(xs, ys)]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制三维点图
ax.scatter(x_new, y_new, zs, s = [15]*len(xs), c='b', marker='o')

# 设置坐标轴标签
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
print(len(x_new))
# 显示图形
plt.show()