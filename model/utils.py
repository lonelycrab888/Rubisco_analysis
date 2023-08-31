import numpy as np
from itertools import product
from itertools import combinations
import matplotlib.pyplot as plt
import math
G = (np.sqrt(5)-1)/2
def getVertex(l):
    G = (np.sqrt(5) - 1) / 2
    pt2 =  [(a,b) for a,b in product([l,-l], [G*l, -G*l])]

    """
    pts =  [(a,b,0) for a,b in pt2]
    pts += [(0,a,b) for a,b in pt2]
    pts += [(b,0,a) for a,b in pt2]
    """
    pts =  [(a,0.934173*b,-0.35682*b) for a,b in pt2]
    pts += [(0,0.934173*a + 0.35682*b, -0.35682*a + 0.934173*b) for a,b in pt2]
    pts += [(b,0.35682*a,0.934173*a) for a,b in pt2]
    """
    pts =  [(0.934173*b,-a,-0.35682*b) for a,b in pt2]
    pts += [(0.934173*a + 0.35682*b,0, -0.35682*a + 0.934173*b) for a,b in pt2]
    pts += [(0.35682*a,-b,0.934173*a) for a,b in pt2]
"""
    return np.array(pts)


def getDisMat(pts):
    N = len(pts)
    dMat = np.ones([N,N])*np.inf
    for i in range(N):
        for j in range(i):
            dMat[i,j] = np.linalg.norm([pts[i]-pts[j]])
    return dMat

def isFace(e1, e2, e3):
    pts = np.vstack([e1, e2, e3])
    pts = np.unique(pts, axis=0)
    return len(pts)==3

def remove_duplicate_sublists(lst):
    seen = set()
    result = []
    for sublist in lst:
        sublist_tuple = tuple(sublist)
        if sublist_tuple not in seen:
            seen.add(sublist_tuple)
            result.append(sublist)
    return result

def generate_points(p1, p2, p3, rows):
    points1 = []
    points2 = []
    r = rows-1
    for i in range(0,r+1):
        x1 = p1[0]+(p2[0]-p1[0])*i/rows
        y1 = p1[1]+(p2[1]-p1[1])*i/rows
        z1 = p1[2]+(p2[2]-p1[2])*i/rows
        flag1 = [x1, y1, z1]


        x2 = p1[0]+(p3[0]-p1[0])*i/rows
        y2 = p1[1]+(p3[1]-p1[1])*i/rows
        z2 = p1[2]+(p3[2]-p1[2])*i/rows
        flag2 = [x2,y2,z2]

        points1.append(flag1)
        points2.append(flag2)
    points = points1 + points2
    for j in range(0,rows-1):
        R_new = rows-j
        R_new = R_new -1
        for i in range(1,R_new):
            x3 = points1[R_new][0] + (points2[R_new][0] - points1[R_new][0]) * i / (R_new)
            y3 = points1[R_new][1] + (points2[R_new][1] - points1[R_new][1]) * i / (R_new)
            z3 = points1[R_new][2] + (points2[R_new][2] - points1[R_new][2]) * i / (R_new)
            flag3 = [x3,y3,z3]
            points.append(flag3)
    return points

def cal_surface(l,n):
    pts = getVertex(l)
    dMat = getDisMat(pts)
    # 由于存在舍入误差，所以得到的边的数值可能不唯一
    ix, jx = np.where((dMat - np.min(dMat)) < 0.01)
    edges = [pts[[i, j]] for i, j in zip(ix, jx)]

    faces = [es for es in combinations(edges, 3) if isFace(*es)]
    pt_list_all = []
    for f in faces:
        pt = np.unique(np.vstack(f), axis=0)
        pt_list = pt.tolist()
        pt_list_all.append(pt_list)
        """
        try:
            ax.plot_trisurf(*pt.T, cmap='viridis', alpha=0.2)
        except:
            pass
        """
    points_all = []
    for p in pt_list_all:
        points = generate_points(p[0], p[1], p[2], n - 1)
        points2 = generate_points_edge(p[0], p[1], p[2], n - 1)
        points_all = points_all + points + points2 + [p[0], p[1], p[2]]

    result = remove_duplicate_sublists(points_all)

    data = result
    # 提取 x、y 和 z 坐标以绘制
    x_values = [point[0] for point in data]
    y_values = [point[1] for point in data]
    z_values = [point[2] for point in data]

    return [x_values, y_values, z_values]


def generate_points_edge(p1, p2, p3, rows):
    points = []
    r = rows-1
    for i in range(0,r+1):
        x1 = p1[0]+(p2[0]-p1[0])*i/rows
        y1 = p1[1]+(p2[1]-p1[1])*i/rows
        z1 = p1[2]+(p2[2]-p1[2])*i/rows
        flag1 = [x1, y1, z1]


        x2 = p1[0]+(p3[0]-p1[0])*i/rows
        y2 = p1[1]+(p3[1]-p1[1])*i/rows
        z2 = p1[2]+(p3[2]-p1[2])*i/rows
        flag2 = [x2,y2,z2]

        x3 = p2[0] + (p3[0] - p2[0]) * i / rows
        y3 = p2[1] + (p3[1] - p2[1]) * i / rows
        z3 = p2[2] + (p3[2] - p2[2]) * i / rows
        flag3 = [x3,y3,z3]

        points.append(flag1)
        points.append(flag2)
        points.append(flag3)

    return points

def add_sphere(x,y,z):
    # 创建球体的参数
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x_true = 5 * np.outer(np.cos(u), np.sin(v))+x
    y_true = 5 * np.outer(np.sin(u), np.sin(v))+y
    z_true = 5 * np.outer(np.ones(np.size(u)), np.cos(v))+z
    return x_true, y_true, z_true

def plot_sphere(center, radius, ax):

    # 创建球体的参数
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # 绘制球体表面
    ax.plot_surface(x, y, z, color='b', alpha=0.3)



# 高斯函数模型
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean)**2) / (2 * stddev**2))
def new_fit(x, x_max):
    x = x*(8/x_max)
    y = (2*gaussian(x, 22.37303119232446 ,15.434947425627788 ,5.2657893542305265)+14*gaussian(x, 8.327540078996142 ,9.139206392757732 ,-3.160802872169968)+ 23*gaussian(x,8.284031865985432 ,9.018061563717959 ,-3.207703935302629)+4*gaussian(x,8.331598620914686 ,8.849643761216303 ,-3.050929770598077)+2*gaussian(x,8.300190861076823 ,8.745882439171481 ,-2.9736606479652945))/45
    y = y*x_max/8
    return y
def new_fit_true(n, x):
    x_max = x*1.07
    x_max = math.ceil(x_max)-1
    y = new_fit(n,x_max)
    return y

def new_prob(x, x_max):
    y = new_fit_true(x, x_max)
    #y = round(y)
    y = y/x
    return y
def new_prob_list(x_max):
    list = [0]
    for i in range(1, x_max+1):
        list.append(new_prob(i, x_max))
    return list

def cal_num(n, x_max):
    num = 10*n*n - 20*n +12
    prob = new_prob_list(x_max)[n]
    return round(num*prob)


def cal_num_all(x_max):
    all = 0
    for i in range(1,x_max+1):
        all = all+cal_num(i,x_max)
    return all

def prob_surface(x_max):
    list = []
    for i in range(1,x_max+1):
        num = 10 * i * i - 20 * i + 12
        list.append(1-(cal_num(i,x_max)/num))
    return list