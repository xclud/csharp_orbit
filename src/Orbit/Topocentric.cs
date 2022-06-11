namespace System;

/// <summary>
/// Encapsulates topo-centric coordinates.
/// </summary>
public readonly ref struct Topocentric<T>
{
    /// <summary>
    /// Creates a new instance of the class from the given components.
    /// </summary>
    /// <param name="azimuth">Azimuth, in radians.</param>
    /// <param name="elevation">Elevation, in radians.</param>
    /// <param name="range">Range, in kilometers.</param>
    /// <param name="rate">Range rate, in kilometers per second. A negative
    /// range rate means "towards the observer".</param>
    public Topocentric(T azimuth, T elevation, T range, T rate)
    {
        Azimuth = azimuth;
        Elevation = elevation;
        Range = range;
        Rate = rate;
    }

    /// <summary>
    /// The azimuth, in radians.
    /// </summary>
    public readonly T Azimuth;

    /// <summary>
    /// The elevation, in radians.
    /// </summary>
    public readonly T Elevation;

    /// <summary>
    /// The range, in kilometers.
    /// </summary>
    public readonly T Range;

    /// <summary>
    /// The range rate, in kilometers per second. 
    /// A negative value means "towards observer".
    /// </summary>
    public readonly T Rate;
}
