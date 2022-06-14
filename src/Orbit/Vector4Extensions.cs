namespace System;

public static class Vector4Extensions
{
    public static double Length(this Vector4<double> v)
    {
        var s = (v.X * v.X) + (v.Y * v.Y) + (v.Z * v.Z) + (v.W * v.W);

        return Math.Sqrt(s);
    }
}
