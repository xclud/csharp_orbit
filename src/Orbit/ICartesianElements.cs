namespace System.Astronomy;

internal interface ICartesianElements
{
    OrbitalState<double> GetPosition(double tsince);
}
