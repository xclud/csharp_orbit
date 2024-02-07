using System.Numerics;

namespace System.Astronomy;

/// <summary>
/// Encapsulates topo-centric coordinates.
/// </summary>
/// <remarks>
/// Creates a new instance of the class from the given components.
/// </remarks>
/// <param name="azimuth">Azimuth, in radians.</param>
/// <param name="elevation">Elevation, in radians.</param>
/// <param name="range">Range, in kilometers.</param>
/// <param name="rate">Range rate, in kilometers per second. A negative
/// range rate means "towards the observer".</param>
public sealed class LookAngle<T>(Angle<T> azimuth, Angle<T> elevation, T range, T rate) where T : INumber<T>, IFloatingPoint<T>
{
    /// <summary>
    /// The azimuth, in radians.
    /// </summary>
    public readonly Angle<T> Azimuth = azimuth;

    /// <summary>
    /// The elevation, in radians.
    /// </summary>
    public readonly Angle<T> Elevation = elevation;

    /// <summary>
    /// The range, in kilometers.
    /// </summary>
    public readonly T Range = range;

    /// <summary>
    /// The range rate, in kilometers per second. 
    /// A negative value means "towards observer".
    /// </summary>
    public readonly T Rate = rate;
}
