# Copyright (c) 2023, 2026.

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

import calendar
import datetime
import logging

import astropy.coordinates as coordinates
import astropy.units as u
import numpy as np
from astroplan import Observer
from astropy.time import Time
from heregoes.navigation._funcs import (
    el2za,
    fractional_jd,
    norm_az,
    za2el,
)
from heregoes.navigation._orbital import get_alt_az

from _funcs import altaz2hadec, angular_separation, solve_sun

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


time_fmt_string = "%Y-%m-%dT%H:%M:%SZ"
decl_min_deg = -23.45
decl_max_deg = 23.45
sun_diameter_deg = 0.53


class FastSparkleAlmanac:
    """
    For fixed-tilt facets or dual-axis tracking solar panels
    """

    almanac_by_time: list[dict]
    almanac_by_separation: list[dict]
    closest_approaches: list[dict]
    sparkle_durations: list[dict]

    def __init__(
        self,
        lat_deg,
        lon_deg,
        sat_za_deg,
        sat_az_deg,
        facet_slope_deg=None,
        facet_azimuth_deg=None,
        facet_height_m=0.0,
        separation_sun_diameters=4.0,
        year_range=range(2021, 2022),
    ):
        self.lat_deg = np.atleast_1d(lat_deg)
        self.lon_deg = np.atleast_1d(lon_deg)
        self.sat_za_deg = np.atleast_1d(sat_za_deg)
        self.sat_za_rad = np.deg2rad(self.sat_za_deg)
        self.sat_az_deg = np.atleast_1d(sat_az_deg)
        self.sat_az_rad = np.deg2rad(self.sat_az_deg)
        self.facet_slope_deg = facet_slope_deg
        self.facet_azimuth_deg = facet_azimuth_deg
        self.facet_height_m = np.atleast_1d(facet_height_m)
        self.year_range = year_range
        self.separation_sun_diameters = separation_sun_diameters

        # settings
        self.sep_thresh_deg = sun_diameter_deg * self.separation_sun_diameters
        self.doy_search_buffer = 90
        self.bisection_time_slice_minutes = 60
        self.bisection_epsilon_seconds = 1
        self.bisection_max_iter = int(
            np.ceil(
                np.log2(
                    self.bisection_time_slice_minutes
                    * 60
                    / self.bisection_epsilon_seconds
                )
            )
        )
        self.fine_iter_time_slice_minutes = 30
        self.fine_time_step_seconds = 30

        # setup
        self.obs = Observer(
            longitude=self.lon_deg * u.deg,
            latitude=self.lat_deg * u.deg,
            elevation=self.facet_height_m * u.m,
        )
        self.earth_position = coordinates.EarthLocation.from_geodetic(
            lon=self.lon_deg * u.deg,
            lat=self.lat_deg * u.deg,
            height=self.facet_height_m * u.m,
            ellipsoid="GRS80",
        )

        if self.facet_slope_deg is not None and self.facet_azimuth_deg is not None:
            # static facet or fixed-tilt solar panel: solve for specular Sun position
            self.tracking_axes = 0
            self.facet_slope_deg = np.atleast_1d(self.facet_slope_deg)
            self.facet_slope_rad = np.deg2rad(self.facet_slope_deg)

            self.facet_azimuth_deg = np.atleast_1d(self.facet_azimuth_deg)
            self.facet_azimuth_rad = np.deg2rad(self.facet_azimuth_deg)

            derived_sun_za_rad, derived_sun_az_rad, derived_sun_omega_rad = solve_sun(
                self.facet_slope_rad,
                self.facet_azimuth_rad,
                self.sat_za_rad,
                self.sat_az_rad,
            )

            self.sol_i_za_rad = derived_sun_za_rad
            self.sol_i_alt_rad = za2el(self.sol_i_za_rad)
            self.sol_i_az_rad = derived_sun_az_rad
            self.sol_i_omega_rad = derived_sun_omega_rad

        else:
            # dual-axis tracking facet: solve for satellite eclipses of the Sun
            self.tracking_axes = 2
            self.sol_i_za_rad = self.sat_za_rad
            self.sol_i_alt_rad = za2el(self.sol_i_za_rad)
            self.sol_i_az_rad = self.sat_az_rad

        self.sol_i_hour_angle, self.sol_i_declination = altaz2hadec(
            self.sol_i_za_rad, self.sol_i_az_rad, self.lat_deg
        )

        if not (
            decl_min_deg - self.sep_thresh_deg
            <= np.rad2deg(self.sol_i_declination)
            <= decl_max_deg + self.sep_thresh_deg
        ):
            raise Exception(
                f"Derived sol_i_declination out of range: {np.rad2deg(self.sol_i_declination)}"
            )

        print("Calculating optimized sparkle almanac for:")
        print("Years:", self.year_range)
        print("Facet latitude:", self.lat_deg)
        print("Facet longitude:", self.lon_deg)
        print("Facet slope:", self.facet_slope_deg)
        print("Facet azimuth:", self.facet_azimuth_deg)
        print("Satellite zenith angle:", self.sat_za_deg)
        print("Satellite azimuth:", self.sat_az_deg)
        print("")

        self.calc_almanac()

    def calc_almanac(self):

        self.almanac_pass = 0
        almanac = []
        self.closest_approaches = []
        self.sparkle_durations = []
        self.trig_bisection_differences = {}
        self.bisection_iteration_differences = {}

        for year in self.year_range:
            self.doys = decl2doy(year, np.rad2deg(self.sol_i_declination))

            # if we got a declination-day solution, construct the target days here
            if len(self.doys) >= 1:
                doy_ranges = []
                for doy in self.doys:
                    doy_start = doy - self.doy_search_buffer
                    doy_stop = doy + self.doy_search_buffer
                    doy_range = np.linspace(
                        doy_start, doy_stop, np.abs(doy_stop - doy_start) + 1
                    ).astype(np.int16)
                    doy_ranges.append(doy_range)

                self.doy_range = np.unique(np.hstack(doy_ranges))

            # otherwise, if the declination is still technically in range, target the entire year
            elif (
                (decl_min_deg - self.sep_thresh_deg)
                <= np.rad2deg(self.sol_i_declination)
                <= (decl_max_deg + self.sep_thresh_deg)
            ):
                year_length = days_in_year(year)
                self.doy_range = np.linspace(
                    1, year_length, np.abs(year_length - 1) + 1
                ).astype(np.int16)

            else:
                raise Exception("Could not construct day list")

            self.date_range = sorted(
                set(
                    [
                        (datetime.datetime.strptime(f"{year} {i}", "%Y %j")).replace(
                            tzinfo=datetime.timezone.utc
                        )
                        for i in (1 + (self.doy_range % days_in_year(year)))
                    ]
                )
            )

            for date in self.date_range:
                self.almanac_pass += 1
                day_almanac = []
                sunrise_time, solar_noon_time, sunset_time = self.get_sun_times(date)

                # this check probably isn't necessary, we're not saving any time
                # abort this day if the altitude of sol_r can never reach the altitude needed for detection
                solar_noon_alt_rad, solar_noon_az_rad = self.get_sun_altaz(
                    solar_noon_time, precise=False
                )
                if self.sol_i_alt_rad > solar_noon_alt_rad + np.deg2rad(
                    self.sep_thresh_deg
                ):
                    logger.debug(
                        "sol_i altitude not in range of day %s",
                        date.strftime(time_fmt_string),
                    )
                    continue

                # solution time "t1" for target hour angle
                noon_offset_hours = (np.rad2deg(self.sol_i_hour_angle) / 15.0).item()
                trig_solution = solar_noon_time + datetime.timedelta(
                    hours=noon_offset_hours
                )
                sol_r_t_alt_rad, sol_r_t_az_rad = self.get_sun_altaz(trig_solution)

                # we're only interested in this for debugging
                # sol_r_hour_angle, sol_r_declination = altaz2hadec(el2za(sol_r_t_alt_rad), sol_r_t_az_rad, self.lat_deg)

                logger.debug(
                    "Got trig solution for date %s: %s",
                    date.strftime(time_fmt_string),
                    trig_solution.strftime(time_fmt_string),
                )  # , np.rad2deg(sol_r_hour_angle), np.rad2deg(sol_r_declination))

                sep = np.rad2deg(
                    angular_separation(
                        sol_r_t_az_rad,
                        sol_r_t_alt_rad,
                        self.sol_i_az_rad,
                        self.sol_i_alt_rad,
                    )
                )

                if sep > self.sep_thresh_deg:
                    logger.debug(
                        "Trig solution for date %s not within %s degrees of sol_i",
                        date.strftime(time_fmt_string),
                        self.sep_thresh_deg,
                    )
                    continue

                # bisection yields a time solution "t2"
                try:
                    self._reset_bisection()
                    bisection_t_start = trig_solution - datetime.timedelta(
                        minutes=self.bisection_time_slice_minutes / 2
                    )
                    bisection_t_stop = trig_solution + datetime.timedelta(
                        minutes=self.bisection_time_slice_minutes / 2
                    )
                    logger.debug(
                        "Starting bisection on trig solution %s from %s to %s",
                        trig_solution.strftime(time_fmt_string),
                        bisection_t_start.strftime(time_fmt_string),
                        bisection_t_stop.strftime(time_fmt_string),
                    )
                    bisection_time_solution = self._bisect(
                        bisection_t_start, bisection_t_stop
                    )
                except Exception as e:
                    logger.debug(e)
                    continue

                # if bisection did not converge, abort
                if self.bisection_converged:
                    trig_bisection_difference = (
                        trig_solution - bisection_time_solution
                    ).total_seconds()
                    self.trig_bisection_differences[date] = trig_bisection_difference
                    logger.debug(
                        "Bisection converged on %s, %s seconds from trig solution",
                        bisection_time_solution.strftime(time_fmt_string),
                        trig_bisection_difference,
                    )

                else:
                    logger.debug(
                        "No bisection solution for day %s",
                        date.strftime(time_fmt_string),
                    )
                    continue

                # iterate over t2 +/- 15 minutes
                logger.debug(
                    "Starting fine iteration on bisection_time_solution=%s",
                    bisection_time_solution.strftime(time_fmt_string),
                )
                day_almanac = self._fine_iter(bisection_time_solution)

                if len(day_almanac) > 0:
                    day_almanac = sorted(day_almanac, key=lambda i: i["time"])
                    almanac.extend(day_almanac)

                    # this assumes only 1 continuous sparkle will occur per day, which should almost always be true
                    sparkle_duration = {
                        "start_time": day_almanac[0]["time"],
                        "stop_time": day_almanac[-1]["time"],
                        "duration_s": (
                            day_almanac[-1]["time"] - day_almanac[0]["time"]
                        ).seconds,
                    }
                    self.sparkle_durations.append(sparkle_duration)

                    closest_approach = sorted(
                        day_almanac, key=lambda i: i["separation"]
                    )[0]
                    closest_approach["duration_s"] = sparkle_duration["duration_s"]
                    self.closest_approaches.append(closest_approach)

                    self.bisection_iteration_differences[date] = (
                        bisection_time_solution - closest_approach["time"]
                    ).total_seconds()
                    if closest_approach["time"] == bisection_time_solution:
                        logger.debug(
                            "Closest approach in the fine iteration (%s) was bisection_time_solution (%s)",
                            closest_approach["time"].strftime(time_fmt_string),
                            bisection_time_solution.strftime(time_fmt_string),
                        )

                    else:
                        logger.debug(
                            "Closest approach in the fine iteration (%s) was NOT bisection_time_solution (%s)",
                            closest_approach["time"].strftime(time_fmt_string),
                            bisection_time_solution.strftime(time_fmt_string),
                        )

        self.almanac_by_time = sorted(almanac, key=lambda i: i["time"])
        self.almanac_by_separation = sorted(almanac, key=lambda i: i["separation"])
        self.closest_approaches = sorted(
            self.closest_approaches, key=lambda i: i["separation"]
        )
        self.sparkle_durations = sorted(
            self.sparkle_durations, key=lambda i: i["duration_s"], reverse=True
        )

        # not needed after sorting
        del almanac

    def _reset_bisection(self):
        self.bisection_iter = 0
        self.bisection_converged = False

    def _bisect(self, time_start, time_stop):
        # trying a bisection on hour angle between -180 and +180
        total_minutes = (time_stop - time_start).total_seconds() / 60
        time_mid = time_start + datetime.timedelta(minutes=total_minutes / 2)

        if (
            time_mid - time_start
        ).total_seconds() <= self.bisection_epsilon_seconds and (
            time_stop - time_mid
        ).total_seconds() <= self.bisection_epsilon_seconds:
            self.bisection_converged = True

        sun_alt_mid, sun_az_mid = self.get_sun_altaz(time_mid, precise=False)
        sol_r_hour_angle_mid, sol_r_declination_mid = altaz2hadec(
            el2za(sun_alt_mid), sun_az_mid, self.lat_deg
        )

        if (
            self.bisection_iter < self.bisection_max_iter
            and not self.bisection_converged
        ):

            self.bisection_iter += 1

            logger.debug(
                "Bisection: Iteration #%s | time_start: %s | time_mid: %s | time_stop: %s",
                self.bisection_iter,
                time_start,
                time_mid,
                time_stop,
            )

            logger.debug(
                "Bisection: Iteration #%s | sol_i_hour_angle: %s | sol_r_mid_hour_angle: %s",
                self.bisection_iter,
                np.rad2deg(self.sol_i_hour_angle),
                np.rad2deg(sol_r_hour_angle_mid),
            )

            # the hour angle always advances clockwise. if the bisection went too far ahead of sol_i, then the difference is positive. otherwise it is negative
            hour_angle_diff = sol_r_hour_angle_mid - self.sol_i_hour_angle

            if np.sign(hour_angle_diff) == np.sign(1):
                time_mid = self._bisect(time_start, time_mid)

            elif np.sign(hour_angle_diff) == np.sign(-1):
                time_mid = self._bisect(time_mid, time_stop)

            else:
                raise Exception("Sun hour angle not found by bisection")

        return time_mid

    def _fine_iter(self, time_mid):
        time_start = time_mid - datetime.timedelta(
            minutes=self.fine_iter_time_slice_minutes / 2
        )
        time_stop = time_mid + datetime.timedelta(
            minutes=self.fine_iter_time_slice_minutes / 2
        )
        total_seconds = (time_stop - time_start).total_seconds()

        def check_positions(check_time, precise=False):

            real_sun_alt_rad, real_sun_az_rad = self.get_sun_altaz(
                check_time, precise=precise
            )

            _ideal_sun_alt_rad = self.sol_i_alt_rad
            _ideal_sun_az_rad = self.sol_i_az_rad

            separation_deg = np.rad2deg(
                angular_separation(
                    real_sun_az_rad,
                    real_sun_alt_rad,
                    _ideal_sun_az_rad,
                    _ideal_sun_alt_rad,
                )
            )

            logger.debug(
                "Time: %s | precise: %s | separation_deg: %s",
                check_time.strftime(time_fmt_string),
                precise,
                separation_deg,
            )

            return separation_deg

        logger.debug(
            "Starting fine iteration for %s seconds between %s and %s",
            str(total_seconds),
            time_start.strftime(time_fmt_string),
            time_stop.strftime(time_fmt_string),
        )

        # "current time" is "t3"
        day_almanac = []
        current_time = time_start
        for second in range(
            0,
            int(total_seconds + 1 + self.fine_time_step_seconds),
            self.fine_time_step_seconds,
        ):
            try:
                separation_deg = check_positions(current_time, precise=True)

            except Exception as e:
                logger.debug(e)
                pass

            else:
                if separation_deg <= self.sep_thresh_deg:
                    # found sparkle
                    day_almanac.append(
                        {
                            "time": current_time,
                            "separation": separation_deg,
                            "almanac_pass": self.almanac_pass,
                        }
                    )

            finally:
                # advance the time by the fine timestep
                current_time += datetime.timedelta(seconds=self.fine_time_step_seconds)

        return day_almanac

    def get_sun_altaz(self, time, precise=False):
        if precise:
            real_sun = coordinates.get_sun(Time(time)).transform_to(
                coordinates.AltAz(obstime=Time(time), location=self.earth_position)
            )

            real_sun_alt_rad = real_sun.alt.rad
            real_sun_az_rad = real_sun.az.rad

        else:
            real_sun_alt_rad, real_sun_az_rad = get_alt_az(
                jdays2000=fractional_jd(time), lon=self.lon_deg, lat=self.lat_deg
            )

            real_sun_az_rad = norm_az(real_sun_az_rad)

        return real_sun_alt_rad, real_sun_az_rad

    def get_sun_times(self, date):
        utc_offset_hours = (self.lon_deg / 15).item()
        # this gets us to roughly the middle of the day at the location. we use this to determine unambiguously the sunrise time in UTC
        solar_noon_guess = Time(
            datetime.datetime(
                year=date.year,
                month=date.month,
                day=date.day,
                hour=12,
                tzinfo=datetime.timezone.utc,
            )
            + datetime.timedelta(hours=-utc_offset_hours)
        )

        sunrise_time = self.obs.sun_rise_time(time=solar_noon_guess, which="previous")
        solar_noon = self.obs.noon(sunrise_time, which="next")
        sunset_time = self.obs.sun_set_time(time=solar_noon, which="next")

        return (
            sunrise_time.datetime.replace(tzinfo=datetime.timezone.utc),
            solar_noon.datetime.replace(tzinfo=datetime.timezone.utc),
            sunset_time.datetime.replace(tzinfo=datetime.timezone.utc),
        )


def days_in_year(year):
    return 365 + calendar.isleap(year)


def doy2decl(year, doy):
    year_length = days_in_year(year)
    return decl_min_deg * np.cos((2 * np.pi / year_length) * ((doy % year_length) + 10))


def decl2doy(year, decl):
    year_length = days_in_year(year)
    doys = []

    nom = year_length * np.arccos(decl / decl_min_deg)
    denom = 2 * np.pi
    res = (nom / denom) - 10
    if ~np.isnan(res):
        doys.append(res)

    nom = year_length * ((2 * np.pi) - np.arccos(decl / decl_min_deg))
    denom = 2 * np.pi
    res = (nom / denom) - 10
    if ~np.isnan(res):
        doys.append(res)

    doys = sorted(set([int(np.asarray(np.floor(i)).item()) for i in doys]))

    return doys
