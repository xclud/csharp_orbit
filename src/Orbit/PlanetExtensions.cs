namespace System;

public static class PlanetExtensions
{
    private const double twoPi = Math.PI * 2;
    private const double pi = Math.PI;
    private const double deg2rad = pi / 180.0;
    private const double rad2deg = 180 / pi;
    private const double x2o3 = 2.0 / 3.0;
    private const double xpdotp = 1440.0 / (2.0 * pi); // 229.1831180523293;

    internal static double xke(this IPlanet p)
    {
        var earthRadius = p.Radius;
        var mu = p.Mu;

        return 60.0 / Math.Sqrt(earthRadius * earthRadius * earthRadius / mu);
    }

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