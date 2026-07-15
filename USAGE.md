```python
from fast import FastSparkleAlmanac

# the Luxor hotel (pyramid) reflecting into GOES-16 from Las Vegas, NV
luxor_g16_almanac = FastSparkleAlmanac(
    lat_deg=36.09559418641608,
    lon_deg=-115.17581583550043,
    sat_za_deg=59.19038009643555,
    sat_az_deg=125.06707000732422,
    facet_slope_deg=47.297,
    facet_azimuth_deg=179.5,
    facet_height_m=636.8504275155227,
)

# Fixed-tilt PV power plant near Boulder, CO reflecting into GOES-17
boulder_g17_almanac = FastSparkleAlmanac(
    lat_deg=40.04809850954546,
    lon_deg=-105.18007220380804,
    sat_za_deg=56.78099917346962,
    sat_az_deg=224.20993560620028,
    facet_slope_deg=31.5,
    facet_azimuth_deg=180.0,
    facet_height_m=1548.09798006095,
)

# greenhouse in Deerfield Township, NJ
nj_g16_almanac = FastSparkleAlmanac(
    lat_deg=39.49810255544305,
    lon_deg=-75.16813121459857,
    sat_za_deg=45.68356172014903,
    sat_az_deg=180.05090522072413,
    facet_slope_deg=48.06531412167497,
    facet_azimuth_deg=180.5605425145889,
    facet_height_m=-3.496409320276328,
)

# agriculture in Chile
chile_g16_almanac = FastSparkleAlmanac(
    lat_deg=-33.91305388344763,
    lon_deg=-71.57313527802935,
    sat_za_deg=39.589080469892245,
    sat_az_deg=353.51315868902327,
    facet_slope_deg=25.0,
    facet_azimuth_deg=1.61269363016629,
    facet_height_m=303.4207563420813,
)

# Fixed-tilt PV solar power plant in Puerto Rico
pr_g16_almanac = FastSparkleAlmanac(
    lat_deg=17.979165659596198,
    lon_deg=-66.2206110036751,
    sat_za_deg=23.46898651123047,
    sat_az_deg=207.13157653808594,
    facet_slope_deg=10,
    facet_azimuth_deg=180.2,
    facet_height_m=-30.509709395211168,
)

# Dual-axis tracking power plant (Alamosa) in Colorado for G16 and G17
dualaxis_g16_almanac = FastSparkleAlmanac(
    lat_deg=37.598538966330395,
    lon_deg=-105.95200431560147,
    sat_za_deg=54.10700988769531,
    sat_az_deg=135.69297790527344,
    facet_slope_deg=None,
    facet_azimuth_deg=None,
    facet_height_m=2295.036997255185,
)
dualaxis_g17_almanac = FastSparkleAlmanac(
    lat_deg=37.598538966330395,
    lon_deg=-105.95200431560147,
    sat_za_deg=54.40988695315823,
    sat_az_deg=224.86881742038847,
    facet_slope_deg=None,
    facet_azimuth_deg=None,
    facet_height_m=2295.036997255185,
)
```

```python
from slow import SparkleAlmanac

#HSAT power plant in Mexico - Parque Solar Don Jose
hsat_g16_almanac = SparkleAlmanac(
    lat_deg=21.34127591454728,
    lon_deg=-100.59481599991655,
    sat_za_deg=38.04798889160156,
    sat_az_deg=127.44476318359375,
    facet_slope_deg=None,
    facet_azimuth_deg=357.64,
    facet_height_m=2011.643882918253,
    tracking_axes=1,
)
```