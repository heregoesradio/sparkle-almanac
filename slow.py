# Copyright (c) 2021, 2022, 2026.

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

import datetime
import logging
import traceback

import astropy.units as u
import numpy as np
from astroplan import Observer
from astropy import coordinates
from astropy.time import Time
from heregoes.navigation._funcs import (
    fractional_jd,
    norm_az,
    za2el,
)
from heregoes.navigation._orbital import get_alt_az

from _funcs import angular_separation, solve_hsat, solve_sun

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class SparkleAlmanac:
    """
    Slower than `FastSparkleAlmanac` but supports predicting HSAT reflections
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
        year_range=range(2022, 2023),
        day_range=range(1, 366),
        tracking_axes=0,
    ):
        self.lat_deg = np.atleast_1d(lat_deg)
        self.lon_deg = np.atleast_1d(lon_deg)

        self.facet_height_m = np.atleast_1d(facet_height_m)

        self.sat_za_deg = np.atleast_1d(sat_za_deg)
        self.sat_za_rad = np.deg2rad(self.sat_za_deg)
        self.sat_az_deg = np.atleast_1d(sat_az_deg)
        self.sat_az_rad = np.deg2rad(self.sat_az_deg)

        self.facet_slope_deg = facet_slope_deg
        if self.facet_slope_deg is not None:
            self.facet_slope_deg = np.atleast_1d(self.facet_slope_deg)
            self.facet_slope_rad = np.deg2rad(self.facet_slope_deg)

        self.facet_azimuth_deg = facet_azimuth_deg
        if self.facet_azimuth_deg is not None:
            self.facet_azimuth_deg = np.atleast_1d(self.facet_azimuth_deg)
            self.facet_azimuth_rad = np.deg2rad(self.facet_azimuth_deg)

        self.year_range = year_range
        self.day_range = day_range
        self.tracking_axes = tracking_axes

        self.obs = Observer(
            longitude=self.lon_deg * u.deg,
            latitude=self.lat_deg * u.deg,
            elevation=facet_height_m * u.m,
        )
        self.earth_position = coordinates.EarthLocation.from_geodetic(
            lon=self.lon_deg * u.deg,
            lat=self.lat_deg * u.deg,
            height=self.facet_height_m * u.m,
            ellipsoid="GRS80",
        )

        self._sun_alt_range = None
        self._sun_az_range = None

        if self.tracking_axes == 0:
            if self.facet_slope_rad is not None and self.facet_azimuth_rad is not None:
                # for a stationary facet, derive the sun position and use it to solve the almanac
                derived_sun_za_rad, derived_sun_az_rad, derived_sun_omega_rad = (
                    solve_sun(
                        self.facet_slope_rad,
                        self.facet_azimuth_rad,
                        self.sat_za_rad,
                        self.sat_az_rad,
                    )
                )

                self.ideal_sun_alt_rad = za2el(derived_sun_za_rad)
                self.ideal_sun_az_rad = derived_sun_az_rad

                # print('sol_i_alt: %s', np.rad2deg(self.ideal_sun_alt_rad))
                # print('sol_i_az: %s', np.rad2deg(self.ideal_sun_az_rad))

                self.check_sun_solution(self.ideal_sun_alt_rad, self.ideal_sun_az_rad)

            else:
                raise Exception(
                    f"facet_slope_deg and facet_azimuth_deg must be defined for tracking_axes={self.tracking_axes}"
                )

        elif self.tracking_axes == 1:
            # for a horizontal azimuthally tracking facet, we only need the facet axis azimuth to calculate specular reflection at each time step
            if self.facet_azimuth_deg is None:
                raise Exception(
                    f"facet_azimuth_deg must be defined for tracking_axes={self.tracking_axes}"
                )

        elif self.tracking_axes == 2:
            # for a dual-axis tracking facet, derive nothing and solve the almanac for satellite eclipses of the sun
            self.ideal_sun_alt_rad = za2el(self.sat_za_rad)
            self.ideal_sun_az_rad = self.sat_az_rad

            # print('sol_i_alt: %s', np.rad2deg(self.ideal_sun_alt_rad))
            # print('sol_i_az: %s', np.rad2deg(self.ideal_sun_az_rad))

        else:
            raise Exception("tracking_axes must be between 0 and 2 inclusive")

        print("Calculating sparkle almanac for:")
        print("Years:", self.year_range)
        print("Days:", self.day_range)
        print("Facet latitude:", self.lat_deg)
        print("Facet longitude:", self.lon_deg)
        print("Facet slope:", self.facet_slope_deg)
        print("Facet azimuth:", self.facet_azimuth_deg)
        print("Facet tracking axes:", self.tracking_axes)
        print("Satellite zenith angle:", self.sat_za_deg)
        print("Satellite azimuth:", self.sat_az_deg)
        print("")

        self.calculate_almanac()

        print("Done!")
        if len(self.closest_approaches) > 0:
            print("Closest approaches:")
            for i in self.closest_approaches[0:10]:
                print(
                    " | ".join(
                        (
                            "Time: " + i["time"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                            "Separation: " + str(np.round(i["separation"], 6)),
                        )
                    )
                )
        else:
            print("No sparkle could be found.")
        print("")

    def calculate_almanac(self, coarse_sun_diameters=20, fine_sun_diameters=15):

        self.almanac_pass = 0
        sun_diameter_deg = 0.53
        self.coarse_separation_threshold_deg = sun_diameter_deg * coarse_sun_diameters
        self.fine_separation_threshold_deg = sun_diameter_deg * fine_sun_diameters

        if self.tracking_axes == 1:
            sun_speed_deg_per_min = 0.50

        else:
            sun_speed_deg_per_min = 0.25

        self.coarse_timestep_minutes = int(
            np.floor(self.coarse_separation_threshold_deg / (sun_speed_deg_per_min * 2))
        )
        self.fine_timestep_seconds = 30

        def check_positions(check_time, precise=False):

            real_sun_alt_rad, real_sun_az_rad = self.get_sun_altaz(
                check_time, precise=precise
            )

            if self.tracking_axes == 1:

                beta, gamma = solve_hsat(
                    real_sun_alt_rad, real_sun_az_rad, self.facet_azimuth_rad
                )

                logger.debug(
                    "tracking_axes: %s | precise: %s | %s | beta_deg: %s | gamma_deg: %s",
                    self.tracking_axes,
                    precise,
                    check_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    np.rad2deg(beta),
                    np.rad2deg(gamma),
                )

                derived_sun_za_rad, derived_sun_az_rad, derived_sun_omega_rad = (
                    solve_sun(beta, gamma, self.sat_za_rad, self.sat_az_rad)
                )

                _ideal_sun_alt_rad = za2el(derived_sun_za_rad)
                _ideal_sun_az_rad = derived_sun_az_rad

                # this will raise an ImpossibleSunSolution if the solution is not possible
                self.check_sun_solution(_ideal_sun_alt_rad, _ideal_sun_az_rad)

            else:
                _ideal_sun_alt_rad = self.ideal_sun_alt_rad
                _ideal_sun_az_rad = self.ideal_sun_az_rad

            separation_deg = np.rad2deg(
                angular_separation(
                    real_sun_az_rad,
                    real_sun_alt_rad,
                    _ideal_sun_az_rad,
                    _ideal_sun_alt_rad,
                )
            )

            logger.debug(
                "tracking_axes: %s | precise: %s | %s | separation_deg: %s",
                self.tracking_axes,
                precise,
                check_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
                separation_deg,
            )
            logger.debug(
                "tracking_axes: %s | precise: %s | %s | sol_r_altaz: (%s,%s)",
                self.tracking_axes,
                precise,
                check_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
                np.rad2deg(real_sun_alt_rad),
                np.rad2deg(real_sun_az_rad),
            )

            return separation_deg

        def iter_day(date):
            sunrise_time, solar_noon, sunset_time = self.get_sun_times(
                date.year, date.month, date.day
            )
            daylight_length_minutes = int(
                np.rint(float((sunset_time - sunrise_time).seconds / 60))
            )

            # by default, iterate over the entire day from sunrise to sunset
            start_minutes = 0
            stop_minutes = daylight_length_minutes
            start_time = sunrise_time

            # #this is a shortcut to skipping the half of the day which does not contain sol_i
            # if self.tracking_axes == 0 or self.tracking_axes == 2:
            #     sunrise_sun_alt_rad, sunrise_sun_az_rad = self.get_sun_altaz(sunrise_time)
            #     solar_noon_sun_alt_rad, solar_noon_sun_az_rad = self.get_sun_altaz(solar_noon)
            #     sunset_sun_alt_rad, sunset_sun_az_rad = self.get_sun_altaz(sunset_time)
            #     buffer_minutes = self.coarse_timestep_minutes * 3

            #     #this is not well determined near the equator. fortunately there are not many sparkle sources there for G16 or G17
            #     if self.lat_deg > 0.0:
            #         before_noon_az_range = (np.rad2deg(sunrise_sun_az_rad), np.rad2deg(solar_noon_sun_az_rad))
            #         after_noon_az_range = (np.rad2deg(solar_noon_sun_az_rad), np.rad2deg(sunset_sun_az_rad))
            #     else:
            #         before_noon_az_range = (np.rad2deg(solar_noon_sun_az_rad), np.rad2deg(sunrise_sun_az_rad))
            #         after_noon_az_range = (np.rad2deg(sunset_sun_az_rad), np.rad2deg(solar_noon_sun_az_rad))

            #     #if the azimuth of sol_i is found before noon, only run on the first half of the day, plus some buffer to catch sparkles that straddle solar noon
            #     if self._angle_between(np.rad2deg(self.ideal_sun_az_rad), *before_noon_az_range):
            #         start_minutes = 0
            #         stop_minutes = (daylight_length_minutes / 2) + buffer_minutes
            #         start_time = sunrise_time

            #     #otherwise if the azimuth of sol_i is found after noon, only run on the second half of the day
            #     elif self._angle_between(np.rad2deg(self.ideal_sun_az_rad), *after_noon_az_range):
            #         start_minutes = (daylight_length_minutes / 2) - buffer_minutes
            #         stop_minutes = daylight_length_minutes
            #         start_time = solar_noon - datetime.timedelta(minutes=buffer_minutes)

            start_minutes = int(start_minutes)
            stop_minutes = int(stop_minutes) + 1 + self.coarse_timestep_minutes

            logger.debug(
                "Starting day iteration on %s for %s minutes from %s to %s",
                date.strftime("%Y-%m-%d"),
                stop_minutes - start_minutes,
                start_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
                (
                    start_time
                    + datetime.timedelta(minutes=stop_minutes - start_minutes)
                ).strftime("%Y-%m-%dT%H:%M:%SZ"),
            )

            # iterate over minutes of daylight (plus coarse_timestep_minutes to ensure full coverage) in steps of coarse_timestep_minutes
            day_almanac = []
            current_time = start_time
            for minute in range(
                start_minutes, stop_minutes, self.coarse_timestep_minutes
            ):

                try:
                    separation_deg = check_positions(current_time, precise=False)

                except ImpossibleSunSolution as e:
                    logger.debug(e)
                    pass

                else:

                    if separation_deg <= self.coarse_separation_threshold_deg:
                        # if we're within coarse distance, reverse 1 coarse timestep
                        current_time -= datetime.timedelta(
                            minutes=self.coarse_timestep_minutes
                        )

                        # then iterate over a coarse timestep in steps of fine_timestep_seconds
                        # this will result in a repetition of the previous coarse check but with precise=True
                        for second in range(
                            0,
                            (60 * self.coarse_timestep_minutes) + 1,
                            self.fine_timestep_seconds,
                        ):

                            try:
                                separation_deg = check_positions(
                                    current_time, precise=True
                                )

                            except ImpossibleSunSolution as e:
                                logger.debug(e)
                                pass

                            else:
                                if separation_deg <= self.fine_separation_threshold_deg:
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
                                current_time += datetime.timedelta(
                                    seconds=self.fine_timestep_seconds
                                )

                finally:
                    # advance the time by the coarse timestep
                    current_time += datetime.timedelta(
                        minutes=self.coarse_timestep_minutes
                    )

            # add this day's almanac to the full almanac, and log the closest approach and totality length
            if len(day_almanac) > 0:
                day_almanac = sorted(day_almanac, key=lambda i: i["time"])
                self.almanac.extend(day_almanac)

                closest_approach = sorted(day_almanac, key=lambda i: i["separation"])[0]
                self.closest_approaches.append(closest_approach)

                # this assumes only 1 continuous sparkle will occur per day, which should almost always be true
                self.sparkle_durations.append(
                    {
                        "start_time": day_almanac[0]["time"],
                        "stop_time": day_almanac[-1]["time"],
                        "duration_s": (
                            day_almanac[-1]["time"] - day_almanac[0]["time"]
                        ).seconds,
                    }
                )

        self.almanac = []
        self.closest_approaches = []
        self.sparkle_durations = []
        try:
            for year in self.year_range:
                for doy in self.day_range:
                    self.almanac_pass += 1
                    date = datetime.datetime.strptime(
                        " ".join((str(year), str(doy))), "%Y %j"
                    )
                    iter_day(date)

        except (Exception, KeyboardInterrupt) as e:
            print(traceback.format_exc())
            raise e

        self.almanac_by_time = sorted(self.almanac, key=lambda i: i["time"])
        self.almanac_by_separation = sorted(self.almanac, key=lambda i: i["separation"])
        self.closest_approaches = sorted(
            self.closest_approaches, key=lambda i: i["separation"]
        )
        self.sparkle_durations = sorted(
            self.sparkle_durations, key=lambda i: i["duration_s"], reverse=True
        )

        # not needed after sorting
        del self.almanac

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

    def check_sun_solution(self, sun_alt_rad, sun_az_rad):
        if self.lat_deg > 23.5 or self.lat_deg < -23.5:
            if (
                not self._angle_between(np.rad2deg(sun_alt_rad), *self.sun_alt_range)
            ) or (not self._angle_between(np.rad2deg(sun_az_rad), *self.sun_az_range)):
                raise ImpossibleSunSolution
        else:
            # the Sun could be at any azimuth or altitude in a given year in the tropics. This is just a shortcut to fail early on an impossible extratropical sparkle
            pass

    def get_sun_times(self, year, month, day):
        utc_offset_hours = (self.lon_deg / 15).item()
        # this gets us to roughly the middle of the day at the location. we use this to determine unambiguously the sunrise time in UTC
        solar_noon_guess = Time(
            datetime.datetime(
                year=int(year),
                month=int(month),
                day=int(day),
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

    @property
    def sun_alt_range(self):
        if self._sun_alt_range is None:
            self._get_sun_altaz_range()

        return self._sun_alt_range

    @property
    def sun_az_range(self):
        if self._sun_az_range is None:
            self._get_sun_altaz_range()

        return self._sun_az_range

    def _get_sun_altaz_range(self):
        # for each possible solstice and equinox, get the location of the sun at sunrise, solar noon, and sunset
        sun_altitudes = []
        sun_azimuths = []
        for year in self.year_range:
            for month in [3, 6, 9, 12]:
                for day in [20, 21, 22, 23]:
                    sunrise_time, solar_noon, sunset_time = self.get_sun_times(
                        year, month, day
                    )

                    sunrise_alt, sunrise_az = self.get_sun_altaz(sunrise_time)
                    solar_noon_alt, solar_noon_az = self.get_sun_altaz(solar_noon)
                    sunset_alt, sunset_az = self.get_sun_altaz(sunset_time)

                    sun_altitudes.append(0)
                    sun_altitudes.append(np.rad2deg(solar_noon_alt))
                    sun_azimuths.append(np.rad2deg(sunrise_az))
                    sun_azimuths.append(np.rad2deg(sunset_az))

        # this should be the range of extreme solar angles at a given location, but it doesn't make sense in the tropics
        if self.lat_deg > 0.0:
            self._sun_az_range = (min(sun_azimuths), max(sun_azimuths))
        else:
            self._sun_az_range = (max(sun_azimuths), min(sun_azimuths))
        self._sun_alt_range = (min(sun_altitudes), max(sun_altitudes))

    @staticmethod
    def _angle_between(n, a, b):
        # returns True if angle n is found between a and b sweeping clockwise
        n = (360 + (n % 360)) % 360
        a = (360 + a) % 360
        b = (360 + b) % 360

        if a < b:
            return (a <= n) & (n <= b)
        else:
            return (a <= n) | (n <= b)


class ImpossibleSunSolution(Exception):
    def __init__(self):
        message = "Could not compute valid solar position for almanac"
        super(ImpossibleSunSolution, self).__init__(message)
