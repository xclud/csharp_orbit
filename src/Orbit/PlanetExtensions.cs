namespace System;

public static class PlanetExtensions
{
    private const double twoPi = Math.PI * 2;
    private const double pi = Math.PI;
    private const double deg2rad = pi / 180.0;
    private const double rad2deg = 180 / pi;
    private const double x2o3 = 2.0 / 3.0;
    private const double xpdotp = 1440.0 / (2.0 * pi); // 229.1831180523293;

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
    internal static OrbitalState<double> GetOrbitalState(this IPlanet planet, Geodetic<double> geocentric, Julian date)
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
        //var W = Math.Sqrt(Globals.Sqr(X) + Globals.Sqr(Y) + Globals.Sqr(Z));


        double mfactor = Globals.TwoPi * (Globals.OmegaE / Globals.SecPerDay);

        var vX = -mfactor * Y;               // km / sec
        var vY = mfactor * X;               // km / sec
        var vZ = 0.0; // km / sec

        // range rate km/sec^2                            
        //var vW = Math.Sqrt(Globals.Sqr(X) + Globals.Sqr(Y));

        var position = new EarthCenteredInertial<double>(X, Y, Z);
        var velocity = new EarthCenteredInertial<double>(vX, vY, vZ);

        return new OrbitalState<double>(position, velocity);
    }

    public static OrbitalState<double> GetOrbitalState(this IPlanet planet, Geodetic<double> geocentric, DateTime utc)
    {
        return GetOrbitalState(planet, geocentric, new Julian(utc.ToUniversalTime()));
    }

    internal static Geodetic<double> GetGeocentric(this IPlanet planet, EarthCenteredInertial<double> position, Julian date)
    {
        var gmst = date.ToGmst();
        var longitude = (gmst - Math.Atan2(position.Y, position.X)) % (Math.PI * 2.0);
        if (longitude > Math.PI)
        {
            longitude -= Math.PI * 2.0;
        }
        if (longitude < -Math.PI)
        {
            longitude += Math.PI * 2.0;
        }


        double num = 1.0;
        var _latitude = Math.PI;
        double num2 = Math.Sqrt(position.X * position.X + position.Y * position.Y);
        _latitude = Math.Atan2(position.Z, num2);
        double latitude;
        do
        {
            latitude = _latitude;
            double num3 = Math.Sin(latitude);
            num = 1.0 / Math.Sqrt(1.0 - 0.0066943177782667227 * num3 * num3);
            _latitude = Math.Atan((position.Z + 6378.136658 * num * 0.0066943177782667227 * num3) / num2);
        }
        while (Math.Abs(_latitude - latitude) > 1E-07);
        var altitude = (num2 < 0.001 ? Math.Abs(position.Z) - 6356.7521724571861 : num2 / Math.Cos(_latitude) - planet.Radius * num);

        return new Geodetic<double>(latitude, longitude, altitude);
    }

    public static Geodetic<double> GetGeocentric(this IPlanet planet, EarthCenteredInertial<double> position, DateTime utc)
    {
        return GetGeocentric(planet, position, new Julian(utc));
    }

    /// <summary>
    /// Returns the topo-centric (azimuth, elevation, etc.) coordinates for
    /// a target object described by the given ECI coordinates.
    /// </summary>
    /// <param name="satellite">The ECI coordinates of the target object.</param>
    /// <returns>The look angle to the target object.</returns>
    public static LookAngle<double> GetLookAngle(this IPlanet planet, Geodetic<double> observer, OrbitalState<double> satellite, DateTime utc)
    {
        return GetLookAngle(planet, observer, satellite, new Julian(utc));
    }

    public static LookAngle<double> GetLookAngle(this IPlanet planet, Geodetic<double> observer, Orbit satellite, DateTime utc)
    {
        return GetLookAngle(planet, observer, satellite.GetPosition(utc), new Julian(utc));
    }
    /// <summary>
    /// Returns the topo-centric (azimuth, elevation, etc.) coordinates for
    /// a target object described by the given ECI coordinates.
    /// </summary>
    /// <param name="satellite">The ECI coordinates of the target object.</param>
    /// <returns>The look angle to the target object.</returns>
    internal static LookAngle<double> GetLookAngle(this IPlanet planet, Geodetic<double> observer, OrbitalState<double> satellite, Julian date)
    {
        // Calculate the ECI coordinates for this Site object at the time
        // of interest.
        var eciSite = planet.GetOrbitalState(observer, date);

        var vx = satellite.Velocity.X - eciSite.Velocity.X;
        var vy = satellite.Velocity.Y - eciSite.Velocity.Y;
        var vz = satellite.Velocity.Z - eciSite.Velocity.Z;

        double rx = satellite.Position.X - eciSite.Position.X;
        double ry = satellite.Position.Y - eciSite.Position.Y;
        double rz = satellite.Position.Z - eciSite.Position.Z;
        double rw = Math.Sqrt(Globals.Sqr(rx) + Globals.Sqr(ry) + Globals.Sqr(rz));

        // The site's Local Mean Sidereal Time at the time of interest.
        double theta = date.ToLmst(observer.Longitude);

        double sin_lat = Math.Sin(observer.Latitude);
        double cos_lat = Math.Cos(observer.Latitude);

        double sin_theta = Math.Sin(theta);
        double cos_theta = Math.Cos(theta);

        double top_s = (sin_lat * cos_theta * rx) +
                       (sin_lat * sin_theta * ry) -
                       (cos_lat * rz);
        double top_e = (-sin_theta * rx) +
                        (cos_theta * ry);
        double top_z = (cos_lat * cos_theta * rx) +
                       (cos_lat * sin_theta * ry) +
                       (sin_lat * rz);
        double az = Math.Atan(-top_e / top_s);

        if (top_s > 0.0)
        {
            az += Globals.Pi;
        }

        if (az < 0.0)
        {
            az += 2.0 * Globals.Pi;
        }

        double el = Math.Asin(top_z / rw);
        double rate = ((rx * vx) + (ry * vy) + (rz * vz)) / rw;

        var topo = new LookAngle<double>(az, el, rw, rate);

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

    //export function radiansToDegrees(radians) {
    //  return radians * rad2deg;
    //}

    //export function degreesToRadians(degrees) {
    //  return degrees * deg2rad;
    //}

    //export function degreesLat(radians) {
    //  if (radians < (-pi / 2) || radians > (pi / 2)) {
    //    throw new RangeError('Latitude radians must be in range [-pi/2; pi/2].');
    //  }
    //  return radiansToDegrees(radians);
    //}

    //export function degreesLong(radians) {
    //  if (radians < -pi || radians > pi) {
    //    throw new RangeError('Longitude radians must be in range [-pi; pi].');
    //  }
    //  return radiansToDegrees(radians);
    //}

    //export function radiansLat(degrees) {
    //  if (degrees < -90 || degrees > 90) {
    //    throw new RangeError('Latitude degrees must be in range [-90; 90].');
    //  }
    //  return degreesToRadians(degrees);
    //}

    //export function radiansLong(degrees) {
    //  if (degrees < -180 || degrees > 180) {
    //    throw new RangeError('Longitude degrees must be in range [-180; 180].');
    //  }
    //  return degreesToRadians(degrees);
    //}

    public static Topocentric<double> Topocentric(this IPlanet planet, Geodetic<double> observer, EarthCenteredEarthFixed<double> satellite)
    {
        // http://www.celestrak.com/columns/v02n02/
        // TS Kelso's method, except I'm using ECF frame
        // and he uses ECI.


        var longitude = observer.Longitude;
        var latitude = observer.Latitude;


        var observerEcf = observer.ToEcf(planet);

        var rx = satellite.X - observerEcf.X;
        var ry = satellite.Y - observerEcf.Y;
        var rz = satellite.Z - observerEcf.Z;

        var south = ((Math.Sin(latitude) * Math.Cos(longitude) * rx) + (Math.Sin(latitude) * Math.Sin(longitude) * ry)) - (Math.Cos(latitude) * rz);

        var topE = (-Math.Sin(longitude) * rx) + (Math.Cos(longitude) * ry);

        var topZ = (Math.Cos(latitude) * Math.Cos(longitude) * rx) + (Math.Cos(latitude) * Math.Sin(longitude) * ry) + (Math.Sin(latitude) * rz);

        return new Topocentric<double>(south, topE, topZ);
    }


    public static LookAngle<double> GetLookAngle(this IPlanet planet, Geodetic<double> observer, EarthCenteredEarthFixed<double> satellite)
    {
        var topocentricCoords = Topocentric(planet, observer, satellite);
        return topocentricCoords.ToLookAngle();
    }

    public static LookAngle<double> GetLookAngle(this IPlanet planet, Geodetic<double> observer, EarthCenteredInertial<double> satellite, double gmst)
    {
        var topocentricCoords = Topocentric(planet, observer, satellite.ToEcf(gmst));
        return topocentricCoords.ToLookAngle();
    }


    public static LookAngle<double> GetLookAngle(this IPlanet planet, Geodetic<double> observer, EarthCenteredInertial<double> satellite, DateTime utc)
    {
        var gmst = new Julian(utc).ToGmst();

        var topocentricCoords = Topocentric(planet, observer, satellite.ToEcf(gmst));
        return topocentricCoords.ToLookAngle();
    }
}