using System.Numerics;

namespace System.Astronomy;

public sealed class KeplerianElements<T>(T epoch, T eccentricity, T inclination, T rightAscensionOfAscendingNode, T argumentOfPeriapsis, T meanMotion, T meanAnomaly, T drag) : IKeplerianElements<T> where T : INumber<T>, IFloatingPoint<T>
{
    public readonly T Epoch = epoch;

    public readonly T Eccentricity = eccentricity;

    public readonly T Inclination = inclination;

    public readonly T RightAscensionOfAscendingNode = rightAscensionOfAscendingNode;

    public readonly T ArgumentOfPeriapsis = argumentOfPeriapsis;

    public readonly T MeanMotion = meanMotion;

    public readonly T MeanAnomaly = meanAnomaly;

    public readonly T Drag = drag;

    T IKeplerianElements<T>.Epoch => Epoch;

    T IKeplerianElements<T>.Eccentricity => Eccentricity;

    T IKeplerianElements<T>.MeanMotion => MeanMotion;

    T IKeplerianElements<T>.Inclination => Inclination;

    T IKeplerianElements<T>.RightAscensionOfAscendingNode => RightAscensionOfAscendingNode;

    T IKeplerianElements<T>.MeanAnomaly => MeanAnomaly;

    T IKeplerianElements<T>.ArgumentOfPeriapsis => ArgumentOfPeriapsis;

    T IKeplerianElements<T>.Drag => Drag;
}