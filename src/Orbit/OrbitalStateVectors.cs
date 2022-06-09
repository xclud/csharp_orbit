using System.Numerics;

namespace System;

/// <summary>
/// In astrodynamics and celestial dynamics, the orbital state vectors of an orbit are Cartesian vectors of position (r) and velocity (v) that together with their time (epoch) (t) uniquely determine the trajectory of the orbiting body in space.
/// 
/// <seealso cref="https://en.wikipedia.org/wiki/Orbital_state_vectors"/>
/// </summary>
public ref struct OrbitalStateVectors
{
    public readonly Vector3 Position;
    public readonly Vector3 Velocity;

    private OrbitalStateVectors(Vector3 position, Vector3 velocity)
    {
        Position = position;
        Velocity = velocity;
    }

    internal static OrbitalStateVectors Create(Vector3 position, Vector3 velocity)
    {
        return new OrbitalStateVectors(position, velocity);
    }
}
