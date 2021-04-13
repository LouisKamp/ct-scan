import matplotlib.pyplot as plt
import numpy as np
import math

from numpy.core.function_base import linspace


# calculates rho for a givin angle phi and intersection point (x,y)
def calc_rho(x: float, y: float, phi: float):
    return round(math.cos(phi) * x + math.sin(phi) * y, 4)


# calculates the boundaries for a givin cell
# N^2 is the number of cells
# a*2/N is the length of the cells
def calc_cell(N, a, i, j):
    x = [-a+2*a*j/N, -a+2*a*(j+1)/N]
    y = [a-2*a*(i+1)/N, a-2*a*i/N]
    return (x, y)


# checks if a line r_l(rho,phi) intersects a cell
# phi lies in the interval [0,pi]
def intersects_cell(i: int, j: int, a: int, N: int, rho: float, phi: float):
    # calcs the cell boundaries
    cell = calc_cell(N, a, i, j)

    q0x = cell[0][0]
    q0y = cell[1][0]

    q1x = cell[0][1]
    q1y = cell[1][1]

    # initialices lower and upper bounds
    rho_lower = None
    rho_upper = None

    # calcs rho values based on the angle
    if (phi <= math.pi/2):
        rho_lower = calc_rho(q0x, q0y, phi)
        rho_upper = calc_rho(q1x, q1y, phi)
        # returns true if r_l(rho,phi) passes throw the cell
        return (rho_lower < rho and rho <= rho_upper)
    else:
        rho_lower = calc_rho(q1x, q0y, phi)
        rho_upper = calc_rho(q0x, q1y, phi)
        # returns true if r_l(rho,phi) passes throw the cell
        return (rho_lower <= rho and rho < rho_upper)


def get_y_from_x(x, rho, phi):
    return round(-(x * math.cos(phi) - rho) / math.sin(phi), 4)


def get_x_from_y(y, rho, phi):
    return round(-(y * math.sin(phi) - rho) / math.cos(phi), 4)


# calculates the length a ray will travel in a given cell
def get_length(i: int, j: int, a: int, N: int, rho: float, phi: float):
    # if the ray doesn't go throw the cell return 0
    if intersects_cell(i, j, a, N, rho, phi) == False:
        return 0

    # if the ray goes throw and has an angle of 0, 90 or 180 we know the length must be 2*a/N
    if (phi == 0 or phi == math.pi or phi == math.pi / 2):
        return 2*a/N

    # define the list of possible points
    points = []

    # calcs the cells x,y-intervals
    cell = calc_cell(N, a, i, j)

    # defines the intervals
    q0x = cell[0][0]
    q0y = cell[1][0]

    q1x = cell[0][1]
    q1y = cell[1][1]

    # calc line intersections and check if they lie on the edge of the cell

    # q0x
    point1_y = get_y_from_x(q0x, rho, phi)

    if (q0y <= point1_y and point1_y <= q1y):
        points.append((q0x, point1_y))

    # q1x
    point2_y = get_y_from_x(q1x, rho, phi)

    if (q0y <= point2_y and point2_y <= q1y):
        points.append((q1x, point2_y))

    # q0y
    point3_x = get_x_from_y(q0y, rho, phi)

    if (q0x <= point3_x and point3_x <= q1x):
        points.append((point3_x, q0y))

    # q1y
    point4_x = get_x_from_y(q1y, rho, phi)

    if (q0x < point4_x and point4_x <= q1x):
        points.append((point4_x, q1y))

    # if there are more than two intersection points, the line goes throw a cornor and we need to sort them
    if (2 < len(points)):
        # find duplicate points
        points = list(set([i for i in points]))

    # in the unfortunate case if there aren't any points, return a length of 0
    if (len(points) == 0 or len(points) == 1):
        return 0

    # calcs the vector between point one and point two
    delta_x = points[0][0] - points[1][0]
    delta_y = points[0][1] - points[1][1]

    # calcs the length of the vector
    l = math.sqrt(delta_x**2+delta_y**2)

    return l


def calc_pos(i: int, j: int, N: int):
    return N*i+j


def calculate_matrix(a: int, N: int, rhos, phis):
    # define matrix
    A = []

    # for every ray
    for rho in rhos:
        for phi in phis:
            # declare new row
            row = np.zeros(N**2, dtype=np.float64)
            for i in range(N):
                for j in range(N):
                    pos = calc_pos(i, j, N)
                    # get the length the ray travels in every cell and assign it to the appropriate position
                    row[pos] = get_length(i, j, a, N, rho, phi)
            # push row on matrix
            A.append(row)
    # return matrix
    return np.array(A)


temp = math.sqrt((0.5**2)*2)

rhos = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
phis = [0]
A1 = calculate_matrix(a=5.5, N=11, rhos=rhos, phis=phis)

rhos = [1, 0, -1]
phis = [-math.pi/2]
A2 = calculate_matrix(a=5.5, N=11, rhos=rhos, phis=phis)

rhos = [-temp*4, -temp*3, -temp*2, -temp*1, 0, temp*1, temp*2, temp*3, temp*4]
phis = [math.pi/4]
A3 = calculate_matrix(a=5.5, N=11, rhos=rhos, phis=phis)

rhos = [-temp*4, -temp*3, -temp*2, -temp*1, 0, temp*1, temp*2, temp*3, temp*4]
phis = [-math.pi/4]
A4 = calculate_matrix(a=5.5, N=11, rhos=rhos, phis=phis)

A = np.concatenate((A1, A2, A3, A4))

# remove exess points
for ray in A:
    for i in range(44):
        ray[i] = 0
for ray in A:
    for i in range(77, 121):
        ray[i] = 0

b1 = [12, 9, 6, 3, 6, 12, 6, 3, 12, 6, 12]
b2 = [32, 26, 29]
b3 = [8.5, 12.7, 4.2, 8.5, 12.7, 8.5, 12.7, 12.7, 12.7]
b4 = [8.5, 12.7, 8.5, 8.5, 12.7, 8.5, 8.5, 8.5, 12.7]

b = np.concatenate((b1, b2, b3, b4))

# find best x
x = np.dot(np.linalg.inv(np.dot(A.T, A) +
                         np.identity(121)*0.01), np.dot(A.T, b))

plt.matshow(np.reshape(x, (11, 11)), vmin=2, vmax=5, cmap='hot')
plt.title("DTU")
plt.show()
