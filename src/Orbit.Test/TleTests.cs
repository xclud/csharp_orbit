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
        var ke = tle.KeplerianElements;

        Assert.AreEqual(ke.Epoch, 80275.98708465m);
        Assert.AreEqual(ke.Inclination, 72.8435m);
        Assert.AreEqual(ke.Eccentricity, 0.0086731m);
        Assert.AreEqual(ke.RightAscensionOfAscendingNode, 115.9689m);
        Assert.AreEqual(ke.ArgumentOfPeriapsis, 52.6988m);
        Assert.AreEqual(ke.MeanAnomaly, 110.5714m);
        Assert.AreEqual(ke.MeanMotion, 16.05824518m);
        Assert.AreEqual(ke.Drag, 0.000066816m);
    }

    [TestMethod]
    public void TestMethod2()
    {
        string str1 = "SGP4 Test";
        string str2 = "1 22824U 93061B   13001.82735048  .00001220  00000-0  51624-3 0  3663";
        string str3 = "2 22824  98.6242 303.1417 0007629 108.0132 252.1969 14.27315768  3864";

        var tle = TwoLineElement<decimal>.Parse(str1, str2, str3);
        var ke = tle.KeplerianElements;

        Assert.AreEqual(ke.Epoch, 13001.82735048m);
        Assert.AreEqual(ke.Drag, 0.00051624m);
        Assert.AreEqual(ke.Inclination, 98.6242m);
    }

    [TestMethod]
    public void TestMethod3()
    {
        string str1 = "TERRA";

        string str2 = "1 26959U 01049C   13002.18143479  .00002415  00000-0  12016-3 0  4325";
        string str3 = "2 26959  97.8817 215.5979 0013092 320.2664 188.6125 15.18597754618041";

        var tle = TwoLineElement<decimal>.Parse(str1, str2, str3);
        var ke = tle.KeplerianElements;

        Assert.AreEqual(ke.MeanAnomaly, 188.6125m);

    }

    [TestMethod]
    public void TBA()
    {
        string str1 = "0 TBA - TO BE ASSIGNED";
        string str2 = "1 T0038U          23343.89355734  .00002702  00000-0  12997-1 0  9990";
        string str3 = "2 T0038 102.5718 355.7923 0205926 353.6386   6.2087 12.68552858186995";

        var tle = TwoLineElement<decimal>.Parse(str1, str2, str3);
        var ke = tle.KeplerianElements;

        Assert.AreEqual(tle.Id, "T0038");

    }

    //[TestMethod]
    //public void LAN()
    //{
    //    string str1 = "Bird 2";

    //    string str2 = "1 26959U 01049C   13002.18143479  .00002415  00000-0  12016-3 0  4325";
    //    string str3 = "2 26959  97.8817 215.5979 0013092 320.2664 188.6125 15.18597754618041";

    //    var tle = TwoLineElement<decimal>.Parse(str1, str2, str3);
    //    var ke = tle;

    //    var lng = ke.LongitudeOfAscendingNode(283.05529587063938m);

    //    Assert.AreEqual(lng, -28.463978093814383m);

    //}
}