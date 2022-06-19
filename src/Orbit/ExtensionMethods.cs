using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace System;

public static class ExtensionMethods
{
    public static EarthCenteredInertial<double> ToEci(this EarthCenteredEarthFixed<double> ecf, double gmst)
    {
        // ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
        //
        // [X]     [C -S  0][X]
        // [Y]  =  [S  C  0][Y]
        // [Z]eci  [0  0  1][Z]ecf
        //
        var x = (ecf.X * Math.Cos(gmst)) - (ecf.Y * Math.Sin(gmst));
        var y = (ecf.X * (Math.Sin(gmst))) + (ecf.Y * Math.Cos(gmst));
        var z = ecf.Z;

        return new EarthCenteredInertial<double>(x: x, y: y, z: z);
    }


    public static EarthCenteredEarthFixed<double> ToEcf(this EarthCenteredInertial<double> eci, double gmst)
    {
        // ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
        //
        // [X]     [C -S  0][X]
        // [Y]  =  [S  C  0][Y]
        // [Z]eci  [0  0  1][Z]ecf
        //
        //
        // Inverse:
        // [X]     [C  S  0][X]
        // [Y]  =  [-S C  0][Y]
        // [Z]ecf  [0  0  1][Z]eci

        var x = (eci.X * Math.Cos(gmst)) + (eci.Y * Math.Sin(gmst));
        var y = (eci.X * (-Math.Sin(gmst))) + (eci.Y * Math.Cos(gmst));
        var z = eci.Z;

        return new EarthCenteredEarthFixed<double>(x, y, z);
    }


    public static EarthCenteredEarthFixed<double> ToEcf(this Geodetic<double> geodetic, IPlanet planet)
    {
        var longitude = geodetic.Longitude;
        var latitude = geodetic.Latitude;
        var height = geodetic.Altitude;

        var a = planet.Radius;
        var f = planet.Flattening;
        var e2 = ((2 * f) - (f * f));
        var normal = a / Math.Sqrt(1 - (e2 * (Math.Sin(latitude) * Math.Sin(latitude))));

        var x = (normal + height) * Math.Cos(latitude) * Math.Cos(longitude);
        var y = (normal + height) * Math.Cos(latitude) * Math.Sin(longitude);
        var z = ((normal * (1 - e2)) + height) * Math.Sin(latitude);

        return new EarthCenteredEarthFixed<double>(x, y, z);
    }

    ///**
    // * @param {Object} tc
    // * @param {Number} tc.topS Positive horizontal vector S due south.
    // * @param {Number} tc.topE Positive horizontal vector E due east.
    // * @param {Number} tc.topZ Vector Z normal to the surface of the earth (up).
    // * @returns {Object}
    // */
    public static LookAngle<double> ToLookAngle(this Topocentric<double> tc)
    {
        var topS = tc.South;
        var topE = tc.East;
        var topZ = tc.Normal;

        // Range in km.
        var range = Math.Sqrt((topS * topS) + (topE * topE) + (topZ * topZ));

        var elevation = Math.Asin(topZ / range);
        var azimuth = Math.Atan2(-topE, topS) + Math.PI;

        return new LookAngle<double>(azimuth: azimuth, elevation: elevation, range: range, rate: 0);
    }

    public static double Length(this EarthCenteredInertial<double> v)
    {
        var s = (v.X * v.X) + (v.Y * v.Y) + (v.Z * v.Z);

        return Math.Sqrt(s);
    }
}
