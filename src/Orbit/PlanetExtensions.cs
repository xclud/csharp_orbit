namespace System;

public static class PlanetExtensions
{
    /// <summary>
    /// Creates a instance of the class from geodetic coordinates.
    /// </summary>
    /// <param name="geocentric">The geocentric coordinates.</param>
    /// <param name="date">The Julian date.</param>
    /// <remarks>
    /// Assumes the Earth is an oblate spheroid.
    /// Reference: The 1992 Astronomical Almanac, page K11
    /// Reference: www.celestrak.com (Dr. T.S. Kelso)
    /// </remarks>
    internal static OrbitalState<double> GetOrbitalState(this IPlanet planet, Geocentric<double> geocentric, Julian date)
    {
        double lat = geocentric.Latitude;
        double lon = geocentric.Longitude;
        double alt = geocentric.Altitude;

        // Calculate Local Mean Sidereal Time (theta)
        double theta = date.ToLmst(lon);
        var F = planet.Flattening;
        var Xkmper = planet.Radius;

        double c = 1.0 / Math.Sqrt(1.0 + (F * (F - 2.0) *
                         Globals.Sqr(Math.Sin(lat))));
        double s = Globals.Sqr(1.0 - F) * c;
        double achcp = ((Xkmper * c) + alt) * Math.Cos(lat);

        var X = achcp * Math.Cos(theta);            // km
        var Y = achcp * Math.Sin(theta);            // km
        var Z = ((Xkmper * s) + alt) * Math.Sin(lat);   // km

        // range, km
        var W = Math.Sqrt(Globals.Sqr(X) + Globals.Sqr(Y) + Globals.Sqr(Z));


        double mfactor = Globals.TwoPi * (Globals.OmegaE / Globals.SecPerDay);

        var vX = -mfactor * Y;               // km / sec
        var vY = mfactor * X;               // km / sec
        var vZ = 0.0;                                 // km / sec
        var vW = Math.Sqrt(Globals.Sqr(X) +  // range rate km/sec^2
                               Globals.Sqr(Y));

        var position = new Vector4<double>(X, Y, Z, W);
        var velocity = new Vector4<double>(vX, vY, vZ, vW);

        return new OrbitalState<double>(position, velocity);
    }

    public static OrbitalState<double> GetOrbitalState(this IPlanet planet, Geocentric<double> geocentric, DateTime utc)
    {
        return GetOrbitalState(planet, geocentric, new Julian(utc));
    }

    /// <summary>
    /// Returns the topo-centric (azimuth, elevation, etc.) coordinates for
    /// a target object described by the given ECI coordinates.
    /// </summary>
    /// <param name="satellite">The ECI coordinates of the target object.</param>
    /// <returns>The look angle to the target object.</returns>
    public static Topocentric<double> GetLookAngle(this IPlanet planet, Geocentric<double> observer, OrbitalState<double> satellite, DateTime utc)
    {
        return GetLookAngle(planet, observer, satellite, new Julian(utc));
    }

    public static Topocentric<double> GetLookAngle(this IPlanet planet, Geocentric<double> observer, Orbit satellite, DateTime utc)
    {
        return GetLookAngle(planet, observer, satellite.GetOrbitalState(utc), new Julian(utc));
    }
    /// <summary>
    /// Returns the topo-centric (azimuth, elevation, etc.) coordinates for
    /// a target object described by the given ECI coordinates.
    /// </summary>
    /// <param name="satellite">The ECI coordinates of the target object.</param>
    /// <returns>The look angle to the target object.</returns>
    internal static Topocentric<double> GetLookAngle(this IPlanet planet, Geocentric<double> observer, OrbitalState<double> satellite, Julian date)
    {
        // Calculate the ECI coordinates for this Site object at the time
        // of interest.
        var eciSite = planet.GetOrbitalState(observer, date);
        var vecRgRate = new Vector4<double>(satellite.Velocity.X - eciSite.Velocity.X,
                                      satellite.Velocity.Y - eciSite.Velocity.Y,
                                      satellite.Velocity.Z - eciSite.Velocity.Z, 1);

        double x = satellite.Position.X - eciSite.Position.X;
        double y = satellite.Position.Y - eciSite.Position.Y;
        double z = satellite.Position.Z - eciSite.Position.Z;
        double w = Math.Sqrt(Globals.Sqr(x) + Globals.Sqr(y) + Globals.Sqr(z));

        var vecRange = new Vector4<double>(x, y, z, w);

        // The site's Local Mean Sidereal Time at the time of interest.
        double theta = date.ToLmst(observer.Longitude);

        double sin_lat = Math.Sin(observer.Latitude);
        double cos_lat = Math.Cos(observer.Latitude);

        double sin_theta = Math.Sin(theta);
        double cos_theta = Math.Cos(theta);

        double top_s = (sin_lat * cos_theta * vecRange.X) +
                       (sin_lat * sin_theta * vecRange.Y) -
                       (cos_lat * vecRange.Z);
        double top_e = (-sin_theta * vecRange.X) +
                        (cos_theta * vecRange.Y);
        double top_z = (cos_lat * cos_theta * vecRange.X) +
                       (cos_lat * sin_theta * vecRange.Y) +
                       (sin_lat * vecRange.Z);
        double az = Math.Atan(-top_e / top_s);

        if (top_s > 0.0)
        {
            az += Globals.Pi;
        }

        if (az < 0.0)
        {
            az += 2.0 * Globals.Pi;
        }

        double el = Math.Asin(top_z / vecRange.W);
        double rate = ((vecRange.X * vecRgRate.X) +
                       (vecRange.Y * vecRgRate.Y) +
                       (vecRange.Z * vecRgRate.Z)) / vecRange.W;

        var topo = new Topocentric<double>(az, el, vecRange.W, rate);

#if WANT_ATMOSPHERIC_CORRECTION
  // Elevation correction for atmospheric refraction.
  // Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104
  // Note:  Correction is meaningless when apparent elevation is below horizon
  topo.ElevationRad += Globals.ToRadians((1.02 / 
                                Math.Tan(Globals.ToRadians(Globals.ToDegrees(el) + 10.3 / 
                                (Globals.ToDegrees(el) + 5.11)))) / 60.0);
  if (topo.ElevationRad < 0.0)
  {
     topo.ElevationRad = el;    // Reset to true elevation
  }

  if (topo.ElevationRad > (Math.PI / 2.0))
  {
     topo.ElevationRad = (Math.PI / 2.0);
  }
#endif
        return topo;
    }

    internal static double xke(this IPlanet p)
    {
        var earthRadius = p.Radius;
        var mu = p.Mu;

        return 60.0 / Math.Sqrt(earthRadius * earthRadius * earthRadius / mu);
    }
}