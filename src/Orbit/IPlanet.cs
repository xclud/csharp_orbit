namespace System.Astronomy;

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

    public static Orbit operator &(IKeplerianElements<double> keplerianElements, IPlanet planet)
    {
        return new Orbit(keplerianElements, planet);
    }

    public static Orbit operator &(IPlanet planet, IKeplerianElements<double> keplerianElements)
    {
        return new Orbit(keplerianElements, planet);
    }
}
