namespace System;

internal static class Globals
{
    #region Constants

    public const double Pi = 3.141592653589793;
    public const double TwoPi = 2.0 * Pi;
    public const double RadsPerDegree = Pi / 180.0;
    public const double DegreesPerRad = 180.0 / Pi;

    public const double Gm = 398601.2;   // Earth gravitational constant, km^3/sec^2
    public const double GeoSyncAlt = 42241.892;  // km
    public const double EarthDiam = 12800.0;    // km
    public const double DaySidereal = (23 * 3600) + (56 * 60) + 4.09;  // sec
    public const double DaySolar = 24 * 3600;   // sec

    public const double Ae = 1.0;
    public const double Au = 149597870.0;  // Astronomical unit (km) (IAU 76)

    public const double HoursPerDay = 24.0;          // Hours per day   (solar)
    public const double MinPerDay = 1440.0;        // Minutes per day (solar)
    public const double SecPerDay = 86400.0;       // Seconds per day (solar)
    public const double OmegaE = 1.00273790934; // Earth rotation per sidereal day

    #endregion

    #region Utility

    // ///////////////////////////////////////////////////////////////////////////
    public static double Sqr(double x)
    {
        return x * x;
    }

    // ///////////////////////////////////////////////////////////////////////////
    public static double Fmod2p(double arg)
    {
        double modu = arg % TwoPi;

        if (modu < 0.0)
        {
            modu += TwoPi;
        }

        return modu;
    }

    // ///////////////////////////////////////////////////////////////////////////
    // Globals.AcTan()
    // ArcTangent of sin(x) / cos(x). The advantage of this function over arctan()
    // is that it returns the correct quadrant of the angle.
    public static double AcTan(double sinx, double cosx)
    {
        double ret = cosx == 0.0 ? sinx > 0.0 ? Pi / 2.0 : 3.0 * Pi / 2.0 : cosx > 0.0 ? Math.Atan(sinx / cosx) : Pi + Math.Atan(sinx / cosx);
        return ret;
    }

    // ///////////////////////////////////////////////////////////////////////////
    public static double ToDegrees(double radians)
    {
        return radians * DegreesPerRad;
    }

    // ///////////////////////////////////////////////////////////////////////////
    public static double ToRadians(double degrees)
    {
        return degrees * RadsPerDegree;
    }

    #endregion
}
