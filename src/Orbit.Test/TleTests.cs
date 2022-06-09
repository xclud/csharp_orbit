namespace Orbit.Test;

[TestClass]
public class TleTests
{
    [TestMethod]
    public void TestMethod1()
    {
        string str1 = "SGP4 Test";
        string str2 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8";
        string str3 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   105";

        var tle = TwoLineElement<decimal>.Parse(str1, str2, str3);

        Assert.AreEqual(tle.Epoch, 80275.98708465m);
        Assert.AreEqual(tle.BStarDrag, 0.000066816m);
        Assert.AreEqual(tle.Inclination, 72.8435m);
        Assert.AreEqual(tle.Eccentricity, 0.0086731m);
        Assert.AreEqual(tle.RightAscensionOfAscendingNode, 115.9689m);
        Assert.AreEqual(tle.ArgumentOfPeriapsis, 52.6988m);
        Assert.AreEqual(tle.MeanAnomaly, 110.5714m);
        Assert.AreEqual(tle.MeanMotion, 16.05824518m);
    }
}