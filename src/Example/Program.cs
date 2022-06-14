string str1 = "SGP4 Test";
string str2 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8";
string str3 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   105";

var tle1 = TwoLineElement<double>.Parse(str1, str2, str3);

PrintPosVel(tle1);

Console.WriteLine();

// Test SDP4
str1 = "SDP4 Test";
str2 = "1 11801U          80230.29629788  .01431103  00000-0  14311-1       8";
str3 = "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     6";

var tle2 = TwoLineElement<double>.Parse(str1, str2, str3);

PrintPosVel(tle2);

Console.WriteLine("\nExample output:");

// Example: Define a location on the earth, then determine the look-angle
// to the SDP4 satellite defined above.

// Create an orbit object using the SDP4 TLE object.
var satellite = tle2 & Earth.WGS72;
//var satellite = new Orbit(tle2, Earth.WGS72);

// Get the location of the satellite from the Orbit object. The 
// earth-centered inertial information is placed into eciSDP4.
// Here we ask for the location of the satellite 90 minutes after
// the TLE epoch.

var when = satellite.Epoch.AddMinutes(90);

// Now create a site object. Site objects represent a location on the 
// surface of the earth. Here we arbitrarily select a point on the
// equator.
// 0.00 N, 100.00 W, 0 km altitude
var observer = new Geocentric<double>(0.0, -100.0 / 180.0 * Math.PI, 0);

// Now get the "look angle" from the site to the satellite. 
var topoLook = Earth.WGS72.GetLookAngle(observer, satellite, when);

// Print out the results. Note that the Azimuth and Elevation are
// stored in the Topocentric object as radians. Here we convert
// to degrees using 180 / Math.PI.
Console.Write($"AZ: {topoLook.Azimuth * 180 / Math.PI:f3}  EL: {topoLook.Elevation * 180 / Math.PI:f3}");


void PrintPosVel(IKeplerianElements<double> tle)
{
    const int Step = 360;

    var sat = new Orbit(tle, Earth.WGS72);
    var coords = new List<OrbitalState<double>>();

    // Calculate position, velocity
    // mpe = "minutes past epoch"
    for (int mpe = 0; mpe <= (Step * 4); mpe += Step)
    {
        // Get the position of the satellite at time "mpe".
        // The coordinates are placed into the variable "eci".
        var eci = sat.GetOrbitalState(mpe);

        // Add the coordinate object to the list
        coords.Add(eci);
    }

    // Print TLE data
    //Console.Write("{0}\n", tle.Name);
    //Console.Write("{0}\n", tle.Line1);
    //Console.Write("{0}\n", tle.Line2);

    // Header
    Console.Write("\n  TSINCE            X                Y                Z\n\n");

    // Iterate over each of the ECI position objects pushed onto the
    // coordinate list, above, printing the ECI position information
    // as we go.
    for (int i = 0; i < coords.Count; i++)
    {
        var e = coords[i];

        Console.Write("{0,8}.00 {1,16:f8} {2,16:f8} {3,16:f8}\n",
                      i * Step,
                      e.Position.X,
                      e.Position.Y,
                      e.Position.Z);
    }

    Console.Write("\n                  XDOT             YDOT             ZDOT\n\n");

    // Iterate over each of the ECI position objects in the coordinate
    // list again, but this time print the velocity information.
    for (int i = 0; i < coords.Count; i++)
    {
        var e = coords[i];

        Console.Write("{0,24:f8} {1,16:f8} {2,16:f8}\n",
                      e.Velocity.X,
                      e.Velocity.Y,
                      e.Velocity.Z);
    }
}