# Copyright (c) 2021, 2022, 2023, 2026.

# Author(s):

#   R. Dove <admin@wx-star.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from heregoes.core import heregoes_njit
from heregoes.navigation._funcs import el2za, norm_az, za2el


@heregoes_njit
def solve_sun(facet_slope_rad, facet_az_rad, sat_za_rad, sat_az_rad):
    # unit vector of facet normal
    n_x = np.sin(facet_slope_rad) * np.cos(facet_az_rad)
    n_y = np.sin(facet_slope_rad) * np.sin(facet_az_rad)
    n_z = np.cos(facet_slope_rad)

    # unit vector of satellite
    r_x = np.sin(sat_za_rad) * np.cos(sat_az_rad)
    r_y = np.sin(sat_za_rad) * np.sin(sat_az_rad)
    r_z = np.cos(sat_za_rad)

    omega = np.arccos((n_x * r_x) + (n_y * r_y) + (n_z * r_z))

    # unit vector of sun
    s_x = 2.0 * np.cos(omega) * n_x - r_x
    s_y = 2.0 * np.cos(omega) * n_y - r_y
    s_z = 2.0 * np.cos(omega) * n_z - r_z

    sun_za_rad = np.arccos(s_z)
    sun_az_rad = norm_az(np.arctan2(s_y, s_x))

    return sun_za_rad, sun_az_rad, omega


@heregoes_njit
def altaz2hadec(za_rad, az_rad, lat_deg):
    # adapted from vallado 2013 algorithm 28
    decl = np.arcsin(
        np.cos(za_rad) * np.sin(np.deg2rad(lat_deg))
        + np.sin(za_rad) * np.cos(np.deg2rad(lat_deg)) * np.cos(az_rad)
    )

    sin_ha = -np.sin(az_rad) * np.sin(za_rad) / np.cos(decl)
    cos_ha = (
        np.cos(np.deg2rad(lat_deg)) * np.cos(za_rad)
        - np.sin(np.deg2rad(lat_deg)) * np.cos(az_rad) * np.sin(za_rad)
    ) / np.cos(decl)

    ha = np.arctan2(sin_ha, cos_ha)

    return ha, decl


@heregoes_njit
def angular_separation(lon1, lat1, lon2, lat2):
    # copied from Astropy (BSD) but its just Vincenty's formula on a sphere
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.arctan2(np.hypot(num1, num2), denominator)


@heregoes_njit
def solve_hsat(sun_alt_rad, sun_az_rad, facet_axis_azimuth_rad):
    # https://www.nrel.gov/docs/fy13osti/58891.pdf Section 5 eqs. 8-10
    # this calculates the rotation angle of a horizontal axis tracking facet,
    # but neglects how real power plants might track slightly differently for energy efficiency / shading

    # NREL takes beta from the zenith, we want elevation/altitude
    r = np.arctan(
        np.tan(el2za(sun_alt_rad)) * np.sin(sun_az_rad - facet_axis_azimuth_rad)
    )
    beta = np.abs(r)
    gamma = facet_axis_azimuth_rad + np.arcsin(np.sin(r) / np.sin(beta))

    beta = za2el(beta)
    gamma = norm_az(gamma)

    return beta, gamma
