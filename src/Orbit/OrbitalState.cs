using System.Numerics;

namespace System;

/// <summary>
/// In astrodynamics and celestial dynamics, the orbital state vectors of an orbit are Cartesian vectors of position (r) and velocity (v) that together with their time (epoch) (t) uniquely determine the trajectory of the orbiting body in space.
/// 
/// <seealso cref="https://en.wikipedia.org/wiki/Orbital_state_vectors"/>
/// </summary>
public sealed class OrbitalState<T> where T : INumber<T>, IFloatingPoint<T>
{
    public readonly EarthCenteredInertial<T> Position;
    public readonly EarthCenteredInertial<T> Velocity;

    public OrbitalState(EarthCenteredInertial<T> position, EarthCenteredInertial<T> velocity)
    {
        Position = position;
        Velocity = velocity;
    }
}