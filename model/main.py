
from utils import*
import random
G = (np.sqrt(5)-1)/2

#input_l =95
#l=input_l*10/12.35

l = 68
#n = math.floor(l/10)-2
ax = plt.subplot(projection='3d')
#n0-n:展示第n0层到第n层
n0 = 7
n = 7

#xs, ys, zs = getVertex(l).T
pts = getVertex(l)
dMat = getDisMat(pts)
# 由于存在舍入误差，所以得到的边的数值可能不唯一
ix, jx = np.where((dMat-np.min(dMat))<0.01)
edges = [pts[[i,j]] for i,j in zip(ix, jx)]

for pt in edges:
    ax.plot(*pt.T, color='red')
#plt.show()
faces = [es for es in combinations(edges, 3) if isFace(*es)]
pt_list_all = []
for f in faces:
    pt = np.unique(np.vstack(f), axis=0)
    pt_list = pt.tolist()
    pt_list_all.append(pt_list)
    try:
        ax.plot_trisurf(*pt.T, cmap='viridis', alpha=0.2)
    except:
        pass
"""
#5 7 10
pt = np.unique(np.vstack(faces[6]), axis=0)
#ax.plot_trisurf(*pt.T, alpha=0.8)
"""






"""

with open('C:/Users/zhr/Desktop/surface128.tbl', 'r') as file:
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
#ax.scatter(x_new, y_new, zs, s = [100]*len(xs), c='b', marker='o')

"""
"""
points_all = []
for p in pt_list_all:
    points = generate_points(p[0], p[1], p[2], n-1)
    points2 = generate_points_edge(p[0], p[1], p[2], n-1)
    print(len(points))
    points_all = points_all + points + points2 +[p[0], p[1], p[2]]

result = remove_duplicate_sublists(points_all)
print(len(result))
data = result

# 提取 x、y 和 z 坐标以绘制
x_values = [point[0] for point in data]
y_values = [point[1] for point in data]
z_values = [point[2] for point in data]
"""



prob = prob_surface(n)
points_all = [[0],[0],[0]]

for i in range(n0-1,n):
    print(i)
    l_true = l*(i+1)/n-1
    points = cal_surface(l_true, i+1)
    num_points_to_delete = int(len(points[0]) * prob[i])
    indices_to_delete = random.sample(range(len(points[0])), num_points_to_delete)
    for sublist in points:
        sublist[:] = [value for index, value in enumerate(sublist) if index not in indices_to_delete]
    points_all[0] = points_all[0]+points[0]
    points_all[1] = points_all[1]+points[1]
    points_all[2] = points_all[2]+points[2]


x_values = points_all[0]
y_values = points_all[1]
z_values = points_all[2]

ax.scatter(x_values, y_values, z_values,  s = [100]*len(x_values), color='blue', marker='o')

# 创建 XY 平面的网格坐标
x_flat = y_flat = range(-120, 120)
X, Y = np.meshgrid(x_flat, y_flat)
Z = np.zeros_like(X)  # XY 平面的 Z 坐标为 0
ax.plot_surface(X, Y, Z, color='yellow', alpha=1)






plt.axis('equal')
# 设置坐标轴标签
"""

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.legend([scatter], ['Points'])
"""

plt.show()

x_fit = np.linspace(1, n, 1000)
y_fit = new_prob(x_fit,n)

plt.plot(x_fit, y_fit, label="Y:The probability of \n filling each layer of Rubisco along \n the edges of an icosahedron.")
plt.xlabel("The nth layer")
plt.ylabel("Y")
plt.legend()
plt.show()
