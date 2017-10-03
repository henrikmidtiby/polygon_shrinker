# MIT License
#
# Copyright (c) 2017 Henrik Skov Midtiby
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import matplotlib.pyplot as plt
from pylab import plot, title, xlabel, ylabel, show, axis, grid
import collections
import numpy as np
import math

TwoEdgesAndAVertex = collections.namedtuple('TwoEdgesAndAVertex', ['edge1', 'edge2', 'vertex'])

# Test polygons
ekrod1 = [[10.25,55.39167,914.4], [10.14583,55.49167,914.4], [10.40000,55.57500,914.4], [10.50417,55.47500,914.4]]
ekrod2 = [[10.14583,55.49167,914.4], [10.05417,55.58333,914.4], [10.20000,55.58333,914.4], [10.30000,55.60833,914.4], [10.36667,55.60833,914.4], [10.40000,55.57500,914.4]]
ekrod3 = [[9.93333,55.58333,914.4], [9.93333,55.73333,914.4], [10.40833,55.73333,914.4], [10.46667,55.63333,914.4], [10.43333,55.60833,914.4], [10.30000,55.60833,914.4], [10.20000,55.58333,914.4]]


def calculate_area_of_polygon(polygon):
    temp = 0
    last_coordinate = polygon[-1]
    for coord in polygon:
        sum_of_x_coords = coord[0] + last_coordinate[0]
        difference_of_y_coords = coord[1] - last_coordinate[1]
        temp += sum_of_x_coords*difference_of_y_coords/2
        last_coordinate = coord
    return temp


def make_polygon_go_counter_clockwise(polygon):
    area = calculate_area_of_polygon(polygon)
    if area < 0:
        return np.array(list(reversed(polygon)), dtype=np.float32)
    return polygon


def get_edges_around_each_vertex(polygon):
    second_last_coordinate = polygon[-2]
    last_coordinate = polygon[-1]
    for coord in polygon:
        edge_one = [second_last_coordinate, last_coordinate]
        edge_two = [last_coordinate, coord]
        yield TwoEdgesAndAVertex(edge_one, edge_two, last_coordinate)
        second_last_coordinate = last_coordinate
        last_coordinate = coord


def get_edges_of_polygon(polygon):
    last_coordinate = polygon[-1]
    for coord in polygon:
        edge_one = [last_coordinate, coord]
        yield edge_one
        last_coordinate = coord


def side_shift_edge(edge, shift):
    start_point = edge[0]
    end_point = edge[1]
    direction = end_point - start_point
    normalised_direction = direction / np.linalg.norm(direction)
    hat_direction = np.array([-normalised_direction[1], normalised_direction[0]])
    shifted_start_point = start_point + shift*hat_direction
    shifted_end_point = end_point + shift*hat_direction
    return [shifted_start_point, shifted_end_point]


def get_edge_intersection_point(edge_one, edge_two):
    coefficient_matrix = np.array([edge_one[1] - edge_one[0], edge_two[1]-edge_two[0]]).transpose()
    constants = edge_two[0] - edge_one[0]
    solution = np.linalg.solve(coefficient_matrix, constants)
    intersection_point_one = edge_one[0] + solution[0] * (edge_one[1] - edge_one[0])
    #intersection_point_two = edge_two[0] + solution[1] * (edge_two[1] - edge_two[0])
    return intersection_point_one


def get_orientation_of_edges(polygon):
    orientations = []
    for edge in get_edges_of_polygon(polygon):
        orientation = direction_of_edge(edge)
        orientations.append(orientation)
    return np.array(orientations, dtype=np.float32)


def plot_edge(edge, color='blue'):
    e1 = []
    e1.append(edge[0][0])
    e1.append(edge[1][0])
    n1 = []
    n1.append(edge[0][1])
    n1.append(edge[1][1])
    plot(e1, n1, color)


def plot_point(point, color='blue'):
    e1 = []
    e1.append(point[0])
    n1 = []
    n1.append(point[1])
    plot(e1, n1, color)


def shift_and_plot_edge(edge, shift, color='blue'):
    edge = side_shift_edge(edge, shift)
    plot_edge(edge, color)


def direction_of_edge(edge):
    direction = math.atan2(edge[1][1] - edge[0][1], edge[1][0] - edge[0][0])
    return direction


def direction_change(edge1, edge2):
    direction_one = direction_of_edge(edge1)
    direction_two = direction_of_edge(edge2)
    change_in_direction = direction_one - direction_two
    if change_in_direction < math.pi:
        change_in_direction += 2*math.pi
    if change_in_direction > math.pi:
        change_in_direction -= 2*math.pi
    return change_in_direction


def shrink_polygon(polygon, offset):
    shrunken_polygon = []
    for edge_pair in get_edges_around_each_vertex(polygon):
        shifted_edge_one = side_shift_edge(edge_pair.edge1, offset)
        shifted_edge_two = side_shift_edge(edge_pair.edge2, offset)
        point = get_edge_intersection_point(shifted_edge_one, shifted_edge_two)
        shrunken_polygon.append(point)
    return shrunken_polygon


def plot_polygon(polygon, color='blue'):
    e1 = []
    n1 = []
    for point in polygon:
        e1.append(point[0])
        n1.append(point[1])
    e1.append(polygon[0][0])
    n1.append(polygon[0][1])

    plot(e1, n1, color)

    title('Testing')
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    axis('equal')
    grid(True)


def main():
    plt.figure(num=1)
    polygon = make_polygon_go_counter_clockwise(np.asarray(ekrod3, dtype=np.float32))
    original_edge_orientations = get_orientation_of_edges(shrink_polygon(polygon[:, 0:2], 0))

    plot_polygon(polygon, 'yellow')

    for shift_magnitude in [0.01, 0.02, 0.03, 0.04]:
        print("shift: %f" % shift_magnitude)
        shrunken_polygon = shrink_polygon(polygon[:, 0:2], shift_magnitude)
        plot_polygon(shrunken_polygon)
        error = np.linalg.norm(get_orientation_of_edges(shrunken_polygon) - original_edge_orientations)
        if error > 0.1:
            print("Error")

    axis('equal')
    show()

main()
