namespace System;

public interface IPlanet
{
    /// <summary>
    /// Radius of the planet in Kilometers.
    /// </summary>
    double Radius { get; }

    double Mu { get; }
    double J2 { get; }
    double J3 { get; }
    double J4 { get; }
    double Flattening { get; }
    double B { get; }
}
