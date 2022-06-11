using System.Numerics;

namespace System;

/// <summary>
/// Encapsulates geocentric coordinates.
/// </summary>
public sealed class Geocentric<T> where T : INumber<T>, IFloatingPoint<T>
{
    /// <summary>
    /// Creates a new instance of the class with the given components.
    /// </summary>
    /// <param name="latitude">Latitude, in radians. Negative values indicate latitude south.</param>
    /// <param name="longitude">Longitude, in radians. Negative value indicate longitude
    /// west.</param>
    /// <param name="altitude">Altitude above the ellipsoid model, in kilometers.</param>
    public Geocentric(T latitude, T longitude, T altitude)
    {
        Latitude = latitude;
        Longitude = longitude;
        Altitude = altitude;
    }

    /// <summary>
    /// Creates a Geo object from a source Geo object.
    /// </summary>
    /// <param name="geo">The source Geo object.</param>
    public Geocentric(Geocentric<T> geo)
    {
        Latitude = geo.Latitude;
        Longitude = geo.Longitude;
        Altitude = geo.Altitude;
    }


    /// <summary>
    /// Latitude, in radians. A negative value indicates latitude south.
    /// </summary>
    public readonly T Latitude;

    /// <summary>
    /// Longitude, in radians. A negative value indicates longitude west.
    /// </summary>
    public readonly T Longitude;

    /// <summary>
    /// Altitude, in kilometers, above the ellipsoid model.
    /// </summary>
    public readonly T Altitude;


    /// <summary>
    /// Converts to a string representation of the form "38.0N 045.0W 500m".
    /// </summary>
    /// <returns>The formatted string.</returns>
    public override string ToString()
    {
        bool latNorth = Latitude >= T.Zero;
        bool lonEast = Longitude >= T.Zero;

        var lat = 180.0 / Math.PI * double.Parse(Latitude.ToString() ?? "0");
        var lng = 180.0 / Math.PI * double.Parse(Longitude.ToString() ?? "0");
        var alt = 1000 * double.Parse(Altitude.ToString() ?? "0");

        // latitude, longitude in degrees and elevation in meters.
        var u = latNorth ? 'N' : 'S';
        var v = lonEast ? 'E' : 'W';

        string str = $"{Math.Abs(lat):00.0}{u} {Math.Abs(lng):000.0}{v} {alt:F0}m";

        return str;
    }
}
