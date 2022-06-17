namespace System;

internal interface ICartesianElements
{
    OrbitalState<double> GetPosition(double tsince);
}
