using System.Numerics;

namespace System;

public sealed class KeplerianElements<T> : IKeplerianElements<T> where T : INumber<T>, IFloatingPoint<T>
{
    public KeplerianElements(T epoch, T eccentricity, T inclination, T rightAscensionOfAscendingNode, T argumentOfPeriapsis, T meanMotion, T meanAnomaly, T drag)
    {
        Epoch = epoch;
        Eccentricity = eccentricity;
        Inclination = inclination;
        RightAscensionOfAscendingNode = rightAscensionOfAscendingNode;
        ArgumentOfPeriapsis = argumentOfPeriapsis;
        MeanMotion = meanMotion;
        MeanAnomaly = meanAnomaly;
        Drag = drag;
    }

    public readonly T Epoch;

    public readonly T Eccentricity;

    public readonly T Inclination;

    public readonly T RightAscensionOfAscendingNode;

    public readonly T ArgumentOfPeriapsis;

    public readonly T MeanMotion;

    public readonly T MeanAnomaly;

    public readonly T Drag;

    T IKeplerianElements<T>.Epoch => Epoch;

    T IKeplerianElements<T>.Eccentricity => Eccentricity;

    T IKeplerianElements<T>.MeanMotion => MeanMotion;

    T IKeplerianElements<T>.Inclination => Inclination;

    T IKeplerianElements<T>.RightAscensionOfAscendingNode => RightAscensionOfAscendingNode;

    T IKeplerianElements<T>.MeanAnomaly => MeanAnomaly;

    T IKeplerianElements<T>.ArgumentOfPeriapsis => ArgumentOfPeriapsis;

    T IKeplerianElements<T>.Drag => Drag;
}