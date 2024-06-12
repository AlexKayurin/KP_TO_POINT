'''
Calculates KP of points based on reference line KP

:copyright: 2024 Aleksandr Kaiurin <akayurin@gmail.com>
:license: MIT, see LICENSE.txt for more details.
'''

import numpy as np
import math


class T_calc:
    '''
    Calculates grid bearing in radians from Delta E, Delta N
    Input format is numpy array [[de, dn], ..., [de,dn]]
    Output format is numpy array [t, ..., t]
    '''
    @staticmethod
    def t(de, dn):
        grbrg = np.zeros((len(de), 3))
        grbrg[:, 0] = de
        grbrg[:, 1] = dn

        # N
        grbrg[:, 2][(grbrg[:, 0] == 0) & (grbrg[:, 1] >= 0)] = 0
        # S
        grbrg[:, 2][(grbrg[:, 0] == 0) & (grbrg[:, 1] < 0)] = math.radians(180)
        # NE
        grbrg[:, 2][(grbrg[:, 0] >= 0) & (grbrg[:, 1] > 0)] = \
            np.arctan(grbrg[:, 0][(grbrg[:, 0] >= 0) & (grbrg[:, 1] > 0)] / grbrg[:, 1][
                (grbrg[:, 0] >= 0) & (grbrg[:, 1] > 0)])
        # SE & SW
        grbrg[:, 2][(grbrg[:, 0] >= 0) & (grbrg[:, 1] < 0)] = \
            np.arctan(grbrg[:, 0][(grbrg[:, 0] >= 0) & (grbrg[:, 1] < 0)] / grbrg[:, 1][
                (grbrg[:, 0] >= 0) & (grbrg[:, 1] < 0)]) + math.radians(180)

        grbrg[:, 2][(grbrg[:, 0] < 0) & (grbrg[:, 1] < 0)] = \
            np.arctan(grbrg[:, 0][(grbrg[:, 0] < 0) & (grbrg[:, 1] < 0)] / grbrg[:, 1][
                (grbrg[:, 0] < 0) & (grbrg[:, 1] < 0)]) + math.radians(180)
        # NW
        grbrg[:, 2][(grbrg[:, 0] < 0) & (grbrg[:, 1] > 0)] = \
            np.arctan(grbrg[:, 0][(grbrg[:, 0] < 0) & (grbrg[:, 1] > 0)] / grbrg[:, 1][
                (grbrg[:, 0] < 0) & (grbrg[:, 1] > 0)]) + math.radians(360)

        return grbrg[:, 2]


class RefLineCalc:
    '''
    Calculates reference line segments parameters: Delta KP, length, KP scale, grid bearing
    Input format is numpy array [[e, n, kp], ..., [e, n, kp]]
    Output format is list of 3 x numpy arrays [length, ..., length], [kp scale, ..., kp scale], [t, ..., t]
    '''
    def __init__(self, refarr):
        self.refarr = refarr

    def __sort_refarr(self):
        # sort ref array by kp
        self.refarr = self.refarr[self.refarr[:, 2].argsort()]

    def calcline(self):
        # reverse / sort
        self.__sort_refarr()

        # easting diff in ref array
        self.ediff = np.diff(self.refarr[:, 0])
        # northing diff in ref array
        self.ndiff = np.diff(self.refarr[:, 1])
        # delta kp of ref segments
        self.segm_deltakp = np.diff(self.refarr[:, 2])

        # segment length in ref array
        self.segm_l = np.power((np.power(self.ediff, 2) + np.power(self.ndiff, 2)), 0.5)
        # kp scale in ref array
        self.segm_kpscale = self.segm_deltakp / self.segm_l
        # t of segment in ref array
        self.segm_t = T_calc.t(self.ediff, self.ndiff)

        return self.segm_l, self.segm_kpscale, self.segm_t


class Kp_to_Point:
    '''
    Calculates KP of given points based on reference line
    Input format is:
        reference line numpy array [[e, n, kp], ..., [e, n, kp]]
        3 x reference line segments parameters numpy arrays (from RefLineCalc) [length, ..., length], [kp scale, ..., kp scale], [t, ..., t]
        single point numpy array [e, n]
        *optional - Max point offset - if point is closer to segment than this value, it is considered to belong to segment / disregards if == 0 (omitted)
    Output format is single point numpy array [e, n, kp]
    '''
    def __init__(self, refarr, ref, point, maxoff):
        self.refarr = refarr
        self.segm_l = ref[0]
        self.segm_kpscale = ref[1]
        self.segm_t = ref[2]
        self.e = point[0]
        self.n = point[1]
        self.maxoff = maxoff

    def calckp(self):
        # sides form point to start/end of ref segment
        self.side1 = ((self.e - self.refarr[:-1, 0]) ** 2 + (self.n - self.refarr[:-1, 1]) ** 2) ** 0.5
        self.side2 = ((self.e - self.refarr[1:, 0]) ** 2 + (self.n - self.refarr[1:, 1]) ** 2) ** 0.5

        # distance from point to ref segment
        self.dist_to_segm = abs((self.refarr[1:, 1] - self.refarr[:-1, 1]) * self.e -
                                (self.refarr[1:, 0] - self.refarr[:-1, 0]) * self.n +
                                self.refarr[1:, 0] * self.refarr[:-1, 1] - self.refarr[1:, 1] * self.refarr[:-1, 0]) / \
                            ((self.refarr[1:, 0] - self.refarr[:-1, 0]) ** 2 + (
                                        self.refarr[1:, 1] - self.refarr[:-1, 1]) ** 2) ** 0.5

        # check if both triangle angles are acute (point belongs to ref segment)
        self.p1 = np.sign(self.segm_l ** 2 + self.side1 ** 2 - self.side2 ** 2)
        self.p2 = np.sign(self.segm_l ** 2 + self.side2 ** 2 - self.side1 ** 2)

        # changing 0 to 1 (0 to positive)
        self.p1[self.p1[:] == 0] = 1
        self.p2[self.p2[:] == 0] = 1
        # changing -1 to 0 (negative to 0)
        self.p1[self.p1[:] == -1] = 0
        self.p2[self.p2[:] == -1] = 0

        # points indexes where p1 => 0 & p2 => 0
        self.ixs = np.where(np.logical_and(self.p1, self.p2))[0]
        # set maxoff as maxoff if provided or max(dist_to_segm) - disregard if omitted
        self.maxoff = np.max(self.dist_to_segm) if self.maxoff == 0 else self.maxoff

        if self.ixs.size != 0 and np.any(self.dist_to_segm[self.ixs] <= self.maxoff):
            # within ref extents - (p1 => 0 & p2 => 0) & dist_to_segment <= max_off
            self.near_segm = self.ixs[np.argmin(self.dist_to_segm[self.ixs])]
        else:
            # outside ref extents - nearest ref segment start (first or last)
            self.near_segm = np.argmin(self.side1)

            # start of ref segment to point de & dn
        self.de = np.array([self.e - self.refarr[self.near_segm, 0]])
        self.dn = np.array([self.n - self.refarr[self.near_segm, 1]])
        # start of ref segment to point t & dt
        self.vector_t = T_calc.t(self.de, self.dn)[0]
        self.dt = self.segm_t[self.near_segm] - self.vector_t
        # start of ref segment to point distance
        self.vector_l = ((self.de ** 2 + self.dn ** 2) ** 0.5)[0]
        # start of ref segment to point projection length (to ref segment)
        self.vector_proj = self.vector_l * math.cos(self.dt) * self.segm_kpscale[self.near_segm]
        # point kp
        self.point_kp = self.refarr[self.near_segm, 2] + self.vector_proj

        return (np.array([self.e, self.n, self.point_kp]))


def go(refarr, pointarr, maxoff=0):
    '''
    Input format is:
        reference line numpy array [[e, n, kp], ..., [e, n, kp]]
        points numpy array [[e, n], ..., [e, n]]
        *optional - Max point offset - if point is closer to segment than this value, it is considered to belong to segment / disregards offset if omitted
    Output format is numpy array [[e, n, kp], ..., [e, n, kp]]
    '''
    ref = RefLineCalc(refarr).calcline()

    newarr = []
    for p in pointarr:
        newarr.append(Kp_to_Point(refarr, ref, p, maxoff).calckp())
    return (np.array(newarr))

