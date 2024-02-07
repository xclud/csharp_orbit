namespace System.Astronomy;

public sealed class SGP4 : ICartesianElements
{
    private const double twoPi = Math.PI * 2;
    private const double pi = Math.PI;
    private const double deg2rad = pi / 180.0;
    private const double rad2deg = 180 / pi;
    private const double x2o3 = 2.0 / 3.0;
    private const double xpdotp = 1440.0 / (2.0 * pi); // 229.1831180523293;

    private readonly IKeplerianElements<double> keplerianElements;
    private readonly IPlanet planet;

    public SGP4(IKeplerianElements<double> keplerianElements, IPlanet planet)
    {
        this.keplerianElements = keplerianElements;
        this.planet = planet;

        this.no = keplerianElements.MeanMotion / xpdotp;
        this.bstar = keplerianElements.Drag;

        this.inclo = keplerianElements.Inclination.Degrees * deg2rad;
        this.nodeo = keplerianElements.RightAscensionOfAscendingNode.Degrees * deg2rad;
        this.argpo = keplerianElements.ArgumentOfPeriapsis.Degrees * deg2rad;
        this.mo = keplerianElements.MeanAnomaly.Degrees * deg2rad;
        this.ecco = keplerianElements.Eccentricity;

        var year = (int)(keplerianElements.Epoch / 1000.0);
        var doy = keplerianElements.Epoch - (year * 1000);

        year += year < 57 ? 2000 : 1900;
        var j = new Julian(year, doy);

        double epoch = j.Value - 2433281.5;


        var earthRadius = planet.Radius;
        var j2 = planet.J2;
        var j3 = planet.J3;
        var j4 = planet.J4;

        var j3oj2 = j3 / j2;

        //  sgp4fix add opsmode
        this.operationmode = 'i';
        this.method = 'n';


        /* ------------------------ initialization --------------------- */
        // sgp4fix divisor for divide by zero check on inclination
        // the old check used 1.0 + Math.Cos(pi-1.0e-9), but then compared it to
        // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
        const double temp4 = 1.5e-12;



        // ------------------------ earth constants -----------------------
        // sgp4fix identify constants and allow alternate values

        var ss = (78.0 / earthRadius) + 1.0;
        // sgp4fix use multiply for speed instead of pow
        var qzms2ttemp = (120.0 - 78.0) / earthRadius;
        var qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;

        var t = 0.0;

        var initlResult = Initl(planet, ecco: this.ecco, epoch: epoch, this.inclo, this.no, this.operationmode);


        var ao = initlResult.ao;
        var con42 = initlResult.con42;
        var cosio = initlResult.cosio;
        var cosio2 = initlResult.cosio2;
        var eccsq = initlResult.eccsq;
        var omeosq = initlResult.omeosq;
        var posq = initlResult.posq;
        var rp = initlResult.rp;
        var rteosq = initlResult.rteosq;
        var sinio = initlResult.sinio;

        this.no = initlResult.no;
        this.con41 = initlResult.con41;
        this.gsto = initlResult.gsto;

        // sgp4fix remove this check as it is unnecessary
        // the mrt check in sgp4 handles decaying satellite cases even if the starting
        // condition is below the surface of te earth
        // if (rp < 1.0)
        // {
        //   printf("// *** satn%d epoch elts sub-orbital ***\n", satn);
        //   this.error = 5;
        // }

        if (omeosq >= 0.0 || this.no >= 0.0)
        {
            this.isimp = 0;
            if (rp < ((220.0 / earthRadius) + 1.0))
            {
                this.isimp = 1;
            }
            var sfour = ss;
            var qzms24 = qzms2t;
            var perige = (rp - 1.0) * earthRadius;

            // - for perigees below 156 km, s and qoms2t are altered -
            if (perige < 156.0)
            {
                sfour = perige - 78.0;
                if (perige < 98.0)
                {
                    sfour = 20.0;
                }

                // sgp4fix use multiply for speed instead of pow
                var qzms24temp = (120.0 - sfour) / earthRadius;
                qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
                sfour = (sfour / earthRadius) + 1.0;
            }
            var pinvsq = 1.0 / posq;

            var tsi = 1.0 / (ao - sfour);
            this.eta = ao * this.ecco * tsi;
            var etasq = this.eta * this.eta;
            var eeta = this.ecco * this.eta;
            var psisq = Math.Abs(1.0 - etasq);
            var coef = qzms24 * (tsi * tsi * tsi * tsi);
            var coef1 = coef / Math.Pow(psisq, 3.5);
            var cc2 = coef1 * this.no * ((ao * (1.0 + (1.5 * etasq) + (eeta * (4.0 + etasq))))
              + (0.375 * j2 * tsi / psisq * this.con41
                * (8.0 + (3.0 * etasq * (8.0 + etasq)))));
            this.cc1 = this.bstar * cc2;
            var cc3 = 0.0;
            if (this.ecco > 1.0e-4)
            {
                cc3 = -2.0 * coef * tsi * j3oj2 * this.no * sinio / this.ecco;
            }
            this.x1mth2 = 1.0 - cosio2;
            this.cc4 = 2.0 * this.no * coef1 * ao * omeosq * (
              (this.eta * (2.0 + (0.5 * etasq)))
            + (this.ecco * (0.5 + (2.0 * etasq)))
              - (j2 * tsi / (ao * psisq)
                * ((-3.0 * this.con41 * (1.0 - (2.0 * eeta) + (etasq * (1.5 - (0.5 * eeta)))))
                  + (0.75 * this.x1mth2
            * ((2.0 * etasq) - (eeta * (1.0 + etasq)))
                    * Math.Cos(2.0 * this.argpo))))
            );
            this.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + (2.75 * (etasq + eeta)) + (eeta * etasq));
            var cosio4 = cosio2 * cosio2;
            var temp1 = 1.5 * j2 * pinvsq * this.no;
            var temp2 = 0.5 * temp1 * j2 * pinvsq;
            var temp3 = -0.46875 * j4 * pinvsq * pinvsq * this.no;
            this.mdot = this.no + (0.5 * temp1 * rteosq * this.con41)
              + (0.0625 * temp2 * rteosq * (13.0 - (78.0 * cosio2) + (137.0 * cosio4)));
            this.argpdot = (-0.5 * temp1 * con42)
              + (0.0625 * temp2 * (7.0 - (114.0 * cosio2) + (395.0 * cosio4)))
              + (temp3 * (3.0 - (36.0 * cosio2) + (49.0 * cosio4)));
            var xhdot1 = -temp1 * cosio;
            this.nodedot = xhdot1 + (((0.5 * temp2 * (4.0 - (19.0 * cosio2)))
              + (2.0 * temp3 * (3.0 - (7.0 * cosio2)))) * cosio);
            var xpidot = this.argpdot + this.nodedot;
            this.omgcof = this.bstar * cc3 * Math.Cos(this.argpo);
            this.xmcof = 0.0;
            if (this.ecco > 1.0e-4)
            {
                this.xmcof = -x2o3 * coef * this.bstar / eeta;
            }
            this.nodecf = 3.5 * omeosq * xhdot1 * this.cc1;
            this.t2cof = 1.5 * this.cc1;

            // sgp4fix for divide by zero with xinco = 180 deg
            this.xlcof = Math.Abs(cosio + 1.0) > 1.5e-12
                ? -0.25 * j3oj2 * sinio * (3.0 + (5.0 * cosio)) / (1.0 + cosio)
                : -0.25 * j3oj2 * sinio * (3.0 + (5.0 * cosio)) / temp4;
            this.aycof = -0.5 * j3oj2 * sinio;

            // sgp4fix use multiply for speed instead of pow
            var delmotemp = 1.0 + (this.eta * Math.Cos(this.mo));
            this.delmo = delmotemp * delmotemp * delmotemp;
            this.sinmao = Math.Sin(this.mo);
            this.x7thm1 = (7.0 * cosio2) - 1.0;

            // --------------- deep space initialization -------------
            if (2 * pi / this.no >= 225.0)
            {
                this.method = 'd';
                this.isimp = 1;
                var tc = 0.0;
                var inclm = this.inclo;

                var dscomr = new DeepSpaceCommon(epoch, ep: this.ecco, this.argpo, tc, this.inclo, this.nodeo, this.no);

                this.e3 = dscomr.e3;
                this.ee2 = dscomr.ee2;
                this.peo = dscomr.peo;
                this.pgho = dscomr.pgho;
                this.pho = dscomr.pho;
                this.pinco = dscomr.pinco;
                this.plo = dscomr.plo;
                this.se2 = dscomr.se2;
                this.se3 = dscomr.se3;
                this.sgh2 = dscomr.sgh2;
                this.sgh3 = dscomr.sgh3;
                this.sgh4 = dscomr.sgh4;
                this.sh2 = dscomr.sh2;
                this.sh3 = dscomr.sh3;
                this.si2 = dscomr.si2;
                this.si3 = dscomr.si3;
                this.sl2 = dscomr.sl2;
                this.sl3 = dscomr.sl3;
                this.sl4 = dscomr.sl4;
                this.xgh2 = dscomr.xgh2;
                this.xgh3 = dscomr.xgh3;
                this.xgh4 = dscomr.xgh4;
                this.xh2 = dscomr.xh2;
                this.xh3 = dscomr.xh3;
                this.xi2 = dscomr.xi2;
                this.xi3 = dscomr.xi3;
                this.xl2 = dscomr.xl2;
                this.xl3 = dscomr.xl3;
                this.xl4 = dscomr.xl4;
                this.zmol = dscomr.zmol;
                this.zmos = dscomr.zmos;

                var sinim = dscomr.sinim;
                var cosim = dscomr.cosim;
                var em = dscomr.em;
                var emsq = dscomr.emsq;
                var s1 = dscomr.s1;
                var s2 = dscomr.s2;
                var s3 = dscomr.s3;
                var s4 = dscomr.s4;
                var s5 = dscomr.s5;
                var ss1 = dscomr.ss1;
                var ss2 = dscomr.ss2;
                var ss3 = dscomr.ss3;
                var ss4 = dscomr.ss4;
                var ss5 = dscomr.ss5;
                var sz1 = dscomr.sz1;
                var sz3 = dscomr.sz3;
                var sz11 = dscomr.sz11;
                var sz13 = dscomr.sz13;
                var sz21 = dscomr.sz21;
                var sz23 = dscomr.sz23;
                var sz31 = dscomr.sz31;
                var sz33 = dscomr.sz33;

                var nm = dscomr.nm;
                var z1 = dscomr.z1;
                var z3 = dscomr.z3;
                var z11 = dscomr.z11;
                var z13 = dscomr.z13;
                var z21 = dscomr.z21;
                var z23 = dscomr.z23;
                var z31 = dscomr.z31;
                var z33 = dscomr.z33;


                var dpperResult = this.Dpper(t: t, init: true, ep: this.ecco, inclp: this.inclo, nodep: this.nodeo, argpp: this.argpo, mp: this.mo, opsmode: this.operationmode);


                this.ecco = dpperResult.ep;
                this.inclo = dpperResult.inclp;
                this.nodeo = dpperResult.nodep;
                this.argpo = dpperResult.argpp;
                this.mo = dpperResult.mp;

                var argpm = 0.0;
                var nodem = 0.0;
                var mm = 0.0;

                var dsinitResult = DsInit(planet, cosim: cosim, emsq: emsq, argpo: this.argpo, s1: s1, s2: s2, s3: s3, s4: s4, s5: s5, sinim: sinim, ss1: ss1, ss2: ss2, ss3: ss3, ss4: ss4, ss5: ss5, sz1: sz1, sz3: sz3, sz11: sz11, sz13: sz13, sz21: sz21, sz23: sz23, sz31: sz31, sz33: sz33, t: t, tc: tc, gsto: this.gsto, mo: this.mo, mdot: this.mdot, no: this.no, nodeo: this.nodeo, nodedot: this.nodedot, xpidot: xpidot, z1: z1, z3: z3, z11: z11, z13: z13, z21: z21, z23: z23, z31: z31, z33: z33, ecco: this.ecco, eccsq: eccsq, em: em, argpm: argpm, inclm: inclm, mm: mm, nm: nm, nodem: nodem, irez: this.irez, atime: this.atime, d2201: this.d2201, d2211: this.d2211, d3210: this.d3210, d3222: this.d3222, d4410: this.d4410, d4422: this.d4422, d5220: this.d5220, d5232: this.d5232, d5421: this.d5421, d5433: this.d5433, dedt: this.dedt, didt: this.didt, dmdt: this.dmdt, dnodt: this.dnodt, domdt: this.domdt, del1: this.del1, del2: this.del2, del3: this.del3, xfact: this.xfact, xlamo: this.xlamo, xli: this.xli, xni: this.xni);

                this.irez = dsinitResult.irez;
                this.atime = dsinitResult.atime;
                this.d2201 = dsinitResult.d2201;
                this.d2211 = dsinitResult.d2211;

                this.d3210 = dsinitResult.d3210;
                this.d3222 = dsinitResult.d3222;
                this.d4410 = dsinitResult.d4410;
                this.d4422 = dsinitResult.d4422;
                this.d5220 = dsinitResult.d5220;

                this.d5232 = dsinitResult.d5232;
                this.d5421 = dsinitResult.d5421;
                this.d5433 = dsinitResult.d5433;
                this.dedt = dsinitResult.dedt;
                this.didt = dsinitResult.didt;

                this.dmdt = dsinitResult.dmdt;
                this.dnodt = dsinitResult.dnodt;
                this.domdt = dsinitResult.domdt;
                this.del1 = dsinitResult.del1;

                this.del2 = dsinitResult.del2;
                this.del3 = dsinitResult.del3;
                this.xfact = dsinitResult.xfact;
                this.xlamo = dsinitResult.xlamo;
                this.xli = dsinitResult.xli;

                this.xni = dsinitResult.xni;
            }

            // ----------- set variables if not deep space -----------
            if (this.isimp != 1)
            {
                var cc1sq = this.cc1 * this.cc1;
                this.d2 = 4.0 * ao * tsi * cc1sq;
                var temp = this.d2 * tsi * this.cc1 / 3.0;
                this.d3 = ((17.0 * ao) + sfour) * temp;
                this.d4 = 0.5 * temp * ao * tsi * ((221.0 * ao) + (31.0 * sfour)) * this.cc1;
                this.t3cof = this.d2 + (2.0 * cc1sq);
                this.t4cof = 0.25 * ((3.0 * this.d3)
                  + (this.cc1 * ((12.0 * this.d2) + (10.0 * cc1sq))));
                this.t5cof = 0.2 * (
                  (3.0 * this.d4)
                  + (12.0 * this.cc1 * this.d3)
                  + (6.0 * this.d2 * this.d2)
                  + (15.0 * cc1sq * ((2.0 * this.d2) + cc1sq))
                );
            }

            /* finally propogate to zero epoch to initialize all others. */
            // sgp4fix take out check to let satellites process until they are actually below earth surface
            // if(this.error == 0)
        }
    }

    internal readonly char method;
    internal readonly double aycof;
    internal readonly double con41;
    internal readonly double cc1;
    internal readonly double cc4;
    internal readonly int isimp;
    internal readonly double cc5;
    internal readonly double d2;
    internal readonly double d3;
    internal readonly double d4;
    internal readonly double delmo;
    internal readonly double eta;
    internal readonly double sinmao;
    internal readonly double argpdot;
    internal readonly double omgcof;

    internal readonly double x1mth2;
    internal readonly double xlcof;
    internal readonly double x7thm1;

    internal readonly double t2cof;
    internal readonly double t3cof;
    internal readonly double t4cof;
    internal readonly double t5cof;
    internal readonly double mdot;
    internal readonly double nodedot;
    internal readonly double xmcof;
    internal readonly double nodecf;
    internal readonly int irez;
    internal readonly char operationmode;
    internal readonly double ecco;
    internal readonly double no;
    internal readonly double gsto;
    internal readonly double d2201;
    internal readonly double d2211;
    internal readonly double d3210;
    internal readonly double d3222;
    internal readonly double d4410;
    internal readonly double d4422;
    internal readonly double d5220;
    internal readonly double d5232;
    internal readonly double d5421;
    internal readonly double d5433;
    internal readonly double dedt;
    internal readonly double del1;
    internal readonly double del2;
    internal readonly double del3;
    internal readonly double didt;
    internal readonly double dmdt;
    internal readonly double dnodt;
    internal readonly double domdt;
    internal readonly double e3;
    internal readonly double ee2;
    internal readonly double peo;
    internal readonly double pgho;
    internal readonly double pho;
    internal readonly double pinco;
    internal readonly double plo;
    internal readonly double se2;
    internal readonly double se3;
    internal readonly double sgh2;
    internal readonly double sgh3;
    internal readonly double sgh4;
    internal readonly double sh2;
    internal readonly double sh3;
    internal readonly double si2;
    internal readonly double si3;
    internal readonly double sl2;
    internal readonly double sl3;
    internal readonly double sl4;
    internal readonly double xfact;
    internal readonly double xgh2;
    internal readonly double xgh3;
    internal readonly double xgh4;
    internal readonly double xh2;
    internal readonly double xh3;
    internal readonly double xi2;
    internal readonly double xi3;
    internal readonly double xl2;
    internal readonly double zmol;
    internal readonly double zmos;
    internal readonly double xlamo;
    internal readonly double atime;
    internal readonly double xli;
    internal readonly double xni;
    internal readonly double xl4;
    internal readonly double bstar;
    internal readonly double argpo;
    internal readonly double inclo;
    internal readonly double mo;
    internal readonly double nodeo;
    internal readonly double xl3;

    /* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements. by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
    private DeepSpaceLongPeriodPeriodicContributions Dpper(double t, double ep, double inclp, double nodep, double argpp, double mp, bool init, char opsmode)
    {
        //  ---------------------- constants -----------------------------
        const double zns = 1.19459e-5;
        const double zes = 0.01675;
        const double znl = 1.5835218e-4;
        const double zel = 0.05490;

        //  --------------- calculate time varying periodics -----------
        var zm = zmos + (zns * t);

        // be sure that the initial call has time set to zero
        if (init)
        {
            zm = zmos;
        }

        var zf = zm + (2.0 * zes * Math.Sin(zm));
        var sinzf = Math.Sin(zf);
        var f2 = (0.5 * sinzf * sinzf) - 0.25;
        var f3 = -0.5 * sinzf * Math.Cos(zf);

        double ses = (se2 * f2) + (se3 * f3);
        double sis = (si2 * f2) + (si3 * f3);
        double sls = (sl2 * f2) + (sl3 * f3) + (sl4 * sinzf);
        double sghs = (sgh2 * f2) + (sgh3 * f3) + (sgh4 * sinzf);
        double shs = (sh2 * f2) + (sh3 * f3);

        zm = zmol + (znl * t);
        if (init)
        {
            zm = zmol;
        }

        zf = zm + (2.0 * zel * Math.Sin(zm));
        sinzf = Math.Sin(zf);
        f2 = (0.5 * sinzf * sinzf) - 0.25;
        f3 = -0.5 * sinzf * Math.Cos(zf);

        var sel = (ee2 * f2) + (e3 * f3);
        var sil = (xi2 * f2) + (xi3 * f3);
        var sll = (xl2 * f2) + (xl3 * f3) + (xl4 * sinzf);
        var sghl = (xgh2 * f2) + (xgh3 * f3) + (xgh4 * sinzf);
        var shll = (xh2 * f2) + (xh3 * f3);

        var pe = ses + sel;
        var pinc = sis + sil;
        var pl = sls + sll;
        var pgh = sghs + sghl;
        var ph = shs + shll;

        if (!init)
        {
            pe -= peo;
            pinc -= pinco;
            pl -= plo;
            pgh -= pgho;
            ph -= pho;
            inclp += pinc;
            ep += pe;
            var sinip = Math.Sin(inclp);
            var cosip = Math.Cos(inclp);

            /* ----------------- apply periodics directly ------------ */
            // sgp4fix for lyddane choice
            // strn3 used original inclination - this is technically feasible
            // gsfc used perturbed inclination - also technically feasible
            // probably best to readjust the 0.2 limit value and limit discontinuity
            // 0.2 rad = 11.45916 deg
            // use next line for original strn3 approach and original inclination
            // if (inclo >= 0.2)
            // use next line for gsfc version and perturbed inclination
            if (inclp >= 0.2)
            {
                ph /= sinip;
                pgh -= cosip * ph;
                argpp += pgh;
                nodep += ph;
                mp += pl;
            }
            else
            {
                //  ---- apply periodics with lyddane modification ----
                var sinop = Math.Sin(nodep);
                var cosop = Math.Cos(nodep);
                var alfdp = sinip * sinop;
                var betdp = sinip * cosop;
                var dalf = (ph * cosop) + (pinc * cosip * sinop);
                var dbet = (-ph * sinop) + (pinc * cosip * cosop);
                alfdp += dalf;
                betdp += dbet;
                nodep %= twoPi;

                //  sgp4fix for afspc written intrinsic functions
                //  nodep used without a trigonometric function ahead
                if (nodep < 0.0 && opsmode == 'a')
                {
                    nodep += twoPi;
                }
                var xls = mp + argpp + (cosip * nodep);
                var dls = pl + pgh - (pinc * nodep * sinip);
                xls += dls;
                var xnoh = nodep;
                nodep = Math.Atan2(alfdp, betdp);

                //  sgp4fix for afspc written intrinsic functions
                //  nodep used without a trigonometric function ahead
                if (nodep < 0.0 && opsmode == 'a')
                {
                    nodep += twoPi;
                }
                if (Math.Abs(xnoh - nodep) > pi)
                {
                    if (nodep < xnoh)
                    {
                        nodep += twoPi;
                    }
                    else
                    {
                        nodep -= twoPi;
                    }
                }
                mp += pl;
                argpp = xls - mp - (cosip * nodep);
            }
        }

        return new DeepSpaceLongPeriodPeriodicContributions(ep: ep, inclp: inclp, nodep: nodep, argpp: argpp, mp: mp);
    }



    /*-----------------------------------------------------------------------------
 *
 *                           procedure initl
 *
 *  this procedure initializes the sgp4 propagator. all the initialization is
 *    consolidated here instead of having multiple loops inside other routines.
 *
 *  author        : david vallado                  719-573-2600   28 jun 2005
 *
 *  inputs        :
 *    ecco        - eccentricity                           0.0 - 1.0
 *    epoch       - epoch time in days from jan 0, 1950. 0 hr
 *    inclo       - inclination of satellite
 *    no          - mean motion of satellite
 *    satn        - satellite number
 *
 *  outputs       :
 *    ainv        - 1.0 / a
 *    ao          - semi major axis
 *    con41       -
 *    con42       - 1.0 - 5.0 cos(i)
 *    cosio       - cosine of inclination
 *    cosio2      - cosio squared
 *    eccsq       - eccentricity squared
 *    method      - flag for deep space                    'd', 'n'
 *    omeosq      - 1.0 - ecco * ecco
 *    posq        - semi-parameter squared
 *    rp          - radius of perigee
 *    rteosq      - square root of (1.0 - ecco*ecco)
 *    sinio       - sine of inclination
 *    gsto        - gst at time of observation               rad
 *    no          - mean motion of satellite
 *
 *  locals        :
 *    ak          -
 *    d1          -
 *    del         -
 *    adel        -
 *    po          -
 *
 *  coupling      :
 *    getgravconst
 *    gstime      - find greenwich sidereal time from the julian date
 *
 *  references    :
 *    hoots, roehrich, norad spacetrack report #3 1980
 *    hoots, norad spacetrack report #6 1986
 *    hoots, schumacher and glover 2004
 *    vallado, crawford, hujsak, kelso  2006
 ----------------------------------------------------------------------------*/
    private static Spg4InitResult Initl(IPlanet planet, double ecco, double epoch, double inclo, double no, char opsmode)
    {
        var j2 = planet.J2;
        var xke = planet.xke();

        // sgp4fix use old way of finding gst
        // ----------------------- earth constants ---------------------
        // sgp4fix identify constants and allow alternate values

        // ------------- calculate auxillary epoch quantities ----------
        var eccsq = ecco * ecco;
        var omeosq = 1.0 - eccsq;
        var rteosq = Math.Sqrt(omeosq);
        var cosio = Math.Cos(inclo);
        var cosio2 = cosio * cosio;

        // ------------------ un-kozai the mean motion -----------------
        var ak = Math.Pow(xke / no, x2o3);
        var d1 = 0.75 * j2 * ((3.0 * cosio2) - 1.0) / (rteosq * omeosq);
        var delPrime = d1 / (ak * ak);
        var adel = ak * (1.0 - (delPrime * delPrime) - (delPrime * ((1.0 / 3.0) + (134.0 * delPrime * delPrime / 81.0))));
        delPrime = d1 / (adel * adel);
        no /= 1.0 + delPrime;

        var ao = Math.Pow(xke / no, x2o3);
        var sinio = Math.Sin(inclo);
        var po = ao * omeosq;
        var con42 = 1.0 - (5.0 * cosio2);
        var con41 = -con42 - cosio2 - cosio2;
        var ainv = 1.0 / ao;
        var posq = po * po;
        var rp = ao * (1.0 - ecco);
        var method = 'n';

        //  sgp4fix modern approach to finding sidereal time
        double gsto;
        if (opsmode == 'a')
        {
            //  sgp4fix use old way of finding gst
            //  count integer number of days from 0 jan 1970
            var ts70 = epoch - 7305.0;
            var ds70 = Math.Floor(ts70 + 1.0e-8);
            var tfrac = ts70 - ds70;

            //  find greenwich location at epoch
            const double c1 = 1.72027916940703639e-2;
            const double thgr70 = 1.7321343856509374;
            const double fk5r = 5.07551419432269442e-15;

            var c1p2p = c1 + twoPi;
            gsto = (thgr70 + (c1 * ds70) + (c1p2p * tfrac) + (ts70 * ts70 * fk5r)) % twoPi;
            if (gsto < 0.0)
            {
                gsto += twoPi;
            }
        }
        else
        {
            gsto = gstime(epoch + 2433281.5);
        }

        return new Spg4InitResult
        {
            no = no,
            method = method,
            ainv = ainv,
            ao = ao,
            con41 = con41,
            con42 = con42,
            cosio = cosio,

            cosio2 = cosio2,
            eccsq = eccsq,
            omeosq = omeosq,
            posq = posq,

            rp = rp,
            rteosq = rteosq,
            sinio = sinio,
            gsto = gsto,
        };
    }

    public class Spg4InitResult
    {
        internal char method;
        internal double no;
        internal double ainv;
        internal double ao;
        internal double con41;
        internal double con42;
        internal double cosio;
        internal double cosio2;
        internal double eccsq;
        internal double omeosq;
        internal double posq;
        internal double rp;
        internal double rteosq;
        internal double sinio;
        internal double gsto;
    }

    private static double gstime(double jdut1)
    {
        var tut1 = (jdut1 - 2451545.0) / 36525.0;

        var temp = (-6.2e-6 * tut1 * tut1 * tut1)
          + (0.093104 * tut1 * tut1)
          + (((876600.0 * 3600) + 8640184.812866) * tut1) + 67310.54841; // # sec
        temp = temp * deg2rad / 240.0 % twoPi; // 360/86400 = 1/240, to deg, to rad

        //  ------------------------ check quadrants ---------------------
        if (temp < 0.0)
        {
            temp += twoPi;
        }

        return temp;
    }


    /*----------------------------------------------------------------------------
     *
     *                             procedure sgp4
     *
     *  this procedure is the sgp4 prediction model from space command. this is an
     *    updated and combined version of sgp4 and sdp4, which were originally
     *    published separately in spacetrack report //3. this version follows the
     *    methodology from the aiaa paper (2006) describing the history and
     *    development of the code.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    satrec  - initialised structure from sgp4init() call.
     *    tsince  - time since epoch (minutes)
     *
     *  outputs       :
     *    r           - position vector                     km
     *    v           - velocity                            km/sec
     *  return code - non-zero on error.
     *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
     *                   2 - mean motion less than 0.0
     *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
     *                   4 - semi-latus rectum < 0.0
     *                   5 - epoch elements are sub-orbital
     *                   6 - satellite has decayed
     *
     *  locals        :
     *    am          -
     *    axnl, aynl        -
     *    betal       -
     *    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
     *    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
     *    cosisq  , cossu   , sinsu   , cosu    , sinu
     *    delm        -
     *    delomg      -
     *    dndt        -
     *    eccm        -
     *    emsq        -
     *    ecose       -
     *    el2         -
     *    eo1         -
     *    eccp        -
     *    esine       -
     *    argpm       -
     *    argpp       -
     *    omgadf      -
     *    pl          -
     *    r           -
     *    rtemsq      -
     *    rdotl       -
     *    rl          -
     *    rvdot       -
     *    rvdotl      -
     *    su          -
     *    t2  , t3   , t4    , tc
     *    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
     *    u   , ux   , uy    , uz     , vx     , vy     , vz
     *    inclm       - inclination
     *    mm          - mean anomaly
     *    nm          - mean motion
     *    nodem       - right asc of ascending node
     *    xinc        -
     *    xincp       -
     *    xl          -
     *    xlm         -
     *    mp          -
     *    xmdf        -
     *    xmx         -
     *    xmy         -
     *    nodedf      -
     *    xnode       -
     *    nodep       -
     *    np          -
     *
     *  coupling      :
     *    getgravconst-
     *    dpper
     *    dspace
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report //3 1980
     *    hoots, norad spacetrack report //6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     ----------------------------------------------------------------------------*/

    public OrbitalState<double> GetPosition(double tsince)
    {
        var satrec = this;
        var planet = satrec.planet;
        var earthRadius = planet.Radius;
        var xke = planet.xke();
        var j2 = planet.J2;
        var j3 = planet.J3;

        var j3oj2 = j3 / j2;
        var vkmpersec = earthRadius * xke / 60.0;

        /* ------------------ set mathematical constants --------------- */
        // sgp4fix divisor for divide by zero check on inclination
        // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
        // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency

        const double temp4 = 1.5e-12;

        // --------------------- clear sgp4 error flag -----------------
        var t = tsince;

        //  ------- update for secular gravity and atmospheric drag -----
        var xmdf = satrec.mo + (satrec.mdot * t);
        var argpdf = satrec.argpo + (satrec.argpdot * t);
        var nodedf = satrec.nodeo + (satrec.nodedot * t);
        var argpm = argpdf;
        var mm = xmdf;
        var t2 = t * t;
        var nodem = nodedf + (satrec.nodecf * t2);
        var tempa = 1.0 - (satrec.cc1 * t);
        var tempe = satrec.bstar * satrec.cc4 * t;
        var templ = satrec.t2cof * t2;

        if (satrec.isimp != 1)
        {
            var delomg = satrec.omgcof * t;
            //  sgp4fix use mutliply for speed instead of pow
            var delmtemp = 1.0 + (satrec.eta * Math.Cos(xmdf));
            var delm = satrec.xmcof * ((delmtemp * delmtemp * delmtemp) - satrec.delmo);
            var tempp = delomg + delm;
            mm = xmdf + tempp;
            argpm = argpdf - tempp;
            var t3 = t2 * t;
            var t4 = t3 * t;
            tempa = tempa - (satrec.d2 * t2) - (satrec.d3 * t3) - (satrec.d4 * t4);
            tempe += satrec.bstar * satrec.cc5 * (Math.Sin(mm) - satrec.sinmao);
            templ = templ + (satrec.t3cof * t3) + (t4 * (satrec.t4cof + (t * satrec.t5cof)));
        }
        var nm = satrec.no;
        var em = satrec.ecco;
        var inclm = satrec.inclo;

        if (satrec.method == 'd')
        {
            var tc = t;

            var dspaceResult = Dspace(satrec: satrec, t: t, tc: tc, em: em, argpm: argpm, inclm: inclm, xli: satrec.xli, mm: mm, xni: satrec.xni, nodem: nodem, nm: nm);

            em = dspaceResult.em;
            argpm = dspaceResult.argpm;
            inclm = dspaceResult.inclm;
            mm = dspaceResult.mm;
            nodem = dspaceResult.nodem;
            nm = dspaceResult.nm;
        }


        if (nm <= 0.0)
        {
            throw new PropagationException($"Error nm {nm}", 2);
        }

        var am = Math.Pow(xke / nm, x2o3) * tempa * tempa;
        nm = xke / Math.Pow(am, 1.5);
        em -= tempe;

        // fix tolerance for error recognition
        // sgp4fix am is fixed from the previous nm check
        if (em is >= 1.0 or < -0.001) // || (am < 0.95)
        {
            throw new PropagationException($"Error em {em}", 1);
        }

        //  sgp4fix fix tolerance to avoid a divide by zero
        if (em < 1.0e-6)
        {
            em = 1.0e-6;
        }
        mm += satrec.no * templ;
        var xlm = mm + argpm + nodem;

        nodem %= twoPi;
        argpm %= twoPi;
        xlm %= twoPi;
        mm = (xlm - argpm - nodem) % twoPi;

        // ----------------- compute extra mean quantities -------------
        var sinim = Math.Sin(inclm);
        var cosim = Math.Cos(inclm);

        // -------------------- add lunar-solar periodics --------------
        var ep = em;
        var xincp = inclm;
        var argpp = argpm;
        var nodep = nodem;
        var mp = mm;
        var sinip = sinim;
        var cosip = cosim;

        if (satrec.method == 'd')
        {
            var dpperResult = satrec.Dpper(t: t, init: false, ep: ep, inclp: xincp, nodep: nodep, argpp: argpp, mp: mp, opsmode: satrec.operationmode);

            xincp = dpperResult.inclp;

            if (xincp < 0.0)
            {
                xincp = -xincp;
                nodep += pi;
                argpp -= pi;
            }
            if (ep is < 0.0 or > 1.0)
            {
                throw new PropagationException($"Error ep {ep}", 3);
            }
        }

        var aycof = satrec.aycof;
        var xlcof = satrec.xlcof;

        //  -------------------- long period periodics ------------------
        if (satrec.method == 'd')
        {
            sinip = Math.Sin(xincp);
            cosip = Math.Cos(xincp);
            aycof = -0.5 * j3oj2 * sinip;

            //  sgp4fix for divide by zero for xincp = 180 deg
            xlcof = Math.Abs(cosip + 1.0) > 1.5e-12
                ? -0.25 * j3oj2 * sinip * (3.0 + (5.0 * cosip)) / (1.0 + cosip)
                : -0.25 * j3oj2 * sinip * (3.0 + (5.0 * cosip)) / temp4;
        }

        var axnl = ep * Math.Cos(argpp);
        var temp = 1.0 / (am * (1.0 - (ep * ep)));
        var aynl = (ep * Math.Sin(argpp)) + (temp * aycof);
        var xl = mp + argpp + nodep + (temp * xlcof * axnl);

        // --------------------- solve kepler's equation ---------------
        var u = (xl - nodep) % twoPi;
        var eo1 = u;
        var tem5 = 9999.9;
        var ktr = 1;

        var coseo1 = 0.0;
        var sineo1 = 0.0;

        //    sgp4fix for kepler iteration
        //    the following iteration needs better limits on corrections
        while (Math.Abs(tem5) >= 1.0e-12 && ktr <= 10)
        {
            sineo1 = Math.Sin(eo1);
            coseo1 = Math.Cos(eo1);
            tem5 = 1.0 - (coseo1 * axnl) - (sineo1 * aynl);
            tem5 = (u - (aynl * coseo1) + (axnl * sineo1) - eo1) / tem5;
            if (Math.Abs(tem5) >= 0.95)
            {
                tem5 = tem5 > 0.0 ? 0.95 : -0.95;
            }
            eo1 += tem5;
            ktr += 1;
        }

        //  ------------- short period preliminary quantities -----------
        var ecose = (axnl * coseo1) + (aynl * sineo1);
        var esine = (axnl * sineo1) - (aynl * coseo1);
        var el2 = (axnl * axnl) + (aynl * aynl);
        var pl = am * (1.0 - el2);
        if (pl < 0.0)
        {
            throw new PropagationException($"Error pl {pl}", 4);
        }

        var rl = am * (1.0 - ecose);
        var rdotl = Math.Sqrt(am) * esine / rl;
        var rvdotl = Math.Sqrt(pl) / rl;
        var betal = Math.Sqrt(1.0 - el2);
        temp = esine / (1.0 + betal);
        var sinu = am / rl * (sineo1 - aynl - (axnl * temp));
        var cosu = am / rl * (coseo1 - axnl + (aynl * temp));
        var su = Math.Atan2(sinu, cosu);
        var sin2u = (cosu + cosu) * sinu;
        var cos2u = 1.0 - (2.0 * sinu * sinu);
        temp = 1.0 / pl;
        var temp1 = 0.5 * j2 * temp;
        var temp2 = temp1 * temp;

        var con41 = satrec.con41;
        var x1mth2 = satrec.x1mth2;
        var x7thm1 = satrec.x7thm1;

        // -------------- update for short period periodics ------------
        if (satrec.method == 'd')
        {
            var cosisq = cosip * cosip;
            con41 = (3.0 * cosisq) - 1.0;
            x1mth2 = 1.0 - cosisq;
            x7thm1 = (7.0 * cosisq) - 1.0;
        }

        var mrt = (rl * (1.0 - (1.5 * temp2 * betal * con41)))
          + (0.5 * temp1 * x1mth2 * cos2u);

        // sgp4fix for decaying satellites
        if (mrt < 1.0)
        {
            throw new PropagationException($"decay condition {mrt}", 6);
        }

        su -= 0.25 * temp2 * x7thm1 * sin2u;
        var xnode = nodep + (1.5 * temp2 * cosip * sin2u);
        var xinc = xincp + (1.5 * temp2 * cosip * sinip * cos2u);
        var mvt = rdotl - (nm * temp1 * x1mth2 * sin2u / xke);
        var rvdot = rvdotl + (nm * temp1 * ((x1mth2 * cos2u) + (1.5 * con41)) / xke);

        // --------------------- orientation vectors -------------------
        var sinsu = Math.Sin(su);
        var cossu = Math.Cos(su);
        var snod = Math.Sin(xnode);
        var cnod = Math.Cos(xnode);
        var sini = Math.Sin(xinc);
        var cosi = Math.Cos(xinc);
        var xmx = -snod * cosi;
        var xmy = cnod * cosi;
        var ux = (xmx * sinsu) + (cnod * cossu);
        var uy = (xmy * sinsu) + (snod * cossu);
        var uz = sini * sinsu;
        var vx = (xmx * cossu) - (cnod * sinsu);
        var vy = (xmy * cossu) - (snod * sinsu);
        var vz = sini * cossu;

        // --------- position and velocity (in km and km/sec) ----------
        var r = new EarthCenteredInertial<double>(mrt * ux * earthRadius, mrt * uy * earthRadius, mrt * uz * earthRadius);
        var v = new EarthCenteredInertial<double>(((mvt * ux) + (rvdot * vx)) * vkmpersec, ((mvt * uy) + (rvdot * vy)) * vkmpersec, ((mvt * uz) + (rvdot * vz)) * vkmpersec);

        return new OrbitalState<double>(r, v);
    }

    public OrbitalState<double> GetPosition(DateTime utc)
    {
        var tsince = keplerianElements.GetMinutesPastEpoch(utc);

        return GetPosition(tsince);
    }

    /*-----------------------------------------------------------------------------
    *
    *                           procedure dsinit
    *
    *  this procedure provides deep space contributions to mean motion dot due
    *    to geopotential resonance with half day and one day orbits.
    *
    *  author        : david vallado                  719-573-2600   28 jun 2005
    *
    *  inputs        :
    *    cosim, sinim-
    *    emsq        - eccentricity squared
    *    argpo       - argument of perigee
    *    s1, s2, s3, s4, s5      -
    *    ss1, ss2, ss3, ss4, ss5 -
    *    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
    *    t           - time
    *    tc          -
    *    gsto        - greenwich sidereal time                   rad
    *    mo          - mean anomaly
    *    mdot        - mean anomaly dot (rate)
    *    no          - mean motion
    *    nodeo       - right ascension of ascending node
    *    nodedot     - right ascension of ascending node dot (rate)
    *    xpidot      -
    *    z1, z3, z11, z13, z21, z23, z31, z33 -
    *    eccm        - eccentricity
    *    argpm       - argument of perigee
    *    inclm       - inclination
    *    mm          - mean anomaly
    *    xn          - mean motion
    *    nodem       - right ascension of ascending node
    *
    *  outputs       :
    *    em          - eccentricity
    *    argpm       - argument of perigee
    *    inclm       - inclination
    *    mm          - mean anomaly
    *    nm          - mean motion
    *    nodem       - right ascension of ascending node
    *    irez        - flag for resonance           0-none, 1-one day, 2-half day
    *    atime       -
    *    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
    *    dedt        -
    *    didt        -
    *    dmdt        -
    *    dndt        -
    *    dnodt       -
    *    domdt       -
    *    del1, del2, del3        -
    *    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
    *    theta       -
    *    xfact       -
    *    xlamo       -
    *    xli         -
    *    xni
    *
    *  locals        :
    *    ainv2       -
    *    aonv        -
    *    cosisq      -
    *    eoc         -
    *    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
    *    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
    *    sini2       -
    *    temp        -
    *    temp1       -
    *    theta       -
    *    xno2        -
    *
    *  coupling      :
    *    getgravconst
    *
    *  references    :
    *    hoots, roehrich, norad spacetrack report #3 1980
    *    hoots, norad spacetrack report #6 1986
    *    hoots, schumacher and glover 2004
    *    vallado, crawford, hujsak, kelso  2006
    ----------------------------------------------------------------------------*/
    private static DsInitResult DsInit(IPlanet planet, double cosim, double argpo, double s1, double s2, double s3, double s4, double s5, double sinim, double ss1, double ss2, double ss3, double ss4, double ss5, double sz1, double sz3, double sz11, double sz13, double sz21, double sz23, double sz31, double sz33, double t, double tc, double gsto, double mo, double mdot, double no, double nodeo, double nodedot, double xpidot, double z1, double z3, double z11, double z13, double z21, double z23, double z31, double z33, double ecco, double eccsq, double emsq, double em, double argpm, double inclm, double mm, double nm, double nodem, int irez, double atime, double d2201, double d2211, double d3210, double d3222, double d4410, double d4422, double d5220, double d5232, double d5421, double d5433, double dedt, double didt, double dmdt, double dnodt, double domdt, double del1, double del2, double del3, double xfact, double xlamo, double xli, double xni)
    {
        var xke = planet.xke();

        const double q22 = 1.7891679e-6;
        const double q31 = 2.1460748e-6;
        const double q33 = 2.2123015e-7;
        const double root22 = 1.7891679e-6;
        const double root44 = 7.3636953e-9;
        const double root54 = 2.1765803e-9;
        const double rptim = 4.37526908801129966e-3; // equates to 7.29211514668855e-5 rad/sec
        const double root32 = 3.7393792e-7;
        const double root52 = 1.1428639e-7;
        const double znl = 1.5835218e-4;
        const double zns = 1.19459e-5;

        // -------------------- deep space initialization ------------
        irez = 0;
        if (nm is < 0.0052359877 and > 0.0034906585)
        {
            irez = 1;
        }
        if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5))
        {
            irez = 2;
        }

        // ------------------------ do solar terms -------------------
        var ses = ss1 * zns * ss5;
        var sis = ss2 * zns * (sz11 + sz13);
        var sls = -zns * ss3 * (sz1 + sz3 - 14.0 - (6.0 * emsq));
        var sghs = ss4 * zns * (sz31 + sz33 - 6.0);
        var shs = -zns * ss2 * (sz21 + sz23);

        // sgp4fix for 180 deg incl
        if (inclm is < 5.2359877e-2 or > (pi - 5.2359877e-2))
        {
            shs = 0.0;
        }
        if (sinim != 0.0)
        {
            shs /= sinim;
        }
        var sgs = sghs - (cosim * shs);

        // ------------------------- do lunar terms ------------------
        dedt = ses + (s1 * znl * s5);
        didt = sis + (s2 * znl * (z11 + z13));
        dmdt = sls - (znl * s3 * (z1 + z3 - 14.0 - (6.0 * emsq)));
        var sghl = s4 * znl * (z31 + z33 - 6.0);
        var shll = -znl * s2 * (z21 + z23);

        // sgp4fix for 180 deg incl
        if (inclm is < 5.2359877e-2 or > (pi - 5.2359877e-2))
        {
            shll = 0.0;
        }
        domdt = sgs + sghl;
        dnodt = shs;
        if (sinim != 0.0)
        {
            domdt -= cosim / sinim * shll;
            dnodt += shll / sinim;
        }

        // ----------- calculate deep space resonance effects --------
        const double dndt = 0.0;
        var theta = (gsto + (tc * rptim)) % twoPi;
        em += dedt * t;
        inclm += didt * t;
        argpm += domdt * t;
        nodem += dnodt * t;
        mm += dmdt * t;

        // sgp4fix for negative inclinations
        // the following if statement should be commented out
        // if (inclm < 0.0)
        // {
        //   inclm  = -inclm;
        //   argpm  = argpm - pi;
        //   nodem = nodem + pi;
        // }

        // -------------- initialize the resonance terms -------------
        if (irez != 0)
        {
            var aonv = nm / Math.Pow(xke, x2o3);

            // ---------- geopotential resonance for 12 hour orbits ------
            if (irez == 2)
            {
                var cosisq = cosim * cosim;
                var emo = em;
                em = ecco;
                var emsqo = emsq;
                emsq = eccsq;
                var eoc = em * emsq;
                var g201 = -0.306 - ((em - 0.64) * 0.440);

                double g211;
                double g310;
                double g322;
                double g410;
                double g422;
                double g520;
                double g533;
                double g521;
                double g532;

                if (em <= 0.65)
                {
                    g211 = 3.616 - (13.2470 * em) + (16.2900 * emsq);
                    g310 = -19.302 + (117.3900 * em) - (228.4190 * emsq) + (156.5910 * eoc);
                    g322 = -18.9068 + (109.7927 * em) - (214.6334 * emsq) + (146.5816 * eoc);
                    g410 = -41.122 + (242.6940 * em) - (471.0940 * emsq) + (313.9530 * eoc);
                    g422 = -146.407 + (841.8800 * em) - (1629.014 * emsq) + (1083.4350 * eoc);
                    g520 = -532.114 + (3017.977 * em) - (5740.032 * emsq) + (3708.2760 * eoc);
                }
                else
                {
                    g211 = -72.099 + (331.819 * em) - (508.738 * emsq) + (266.724 * eoc);
                    g310 = -346.844 + (1582.851 * em) - (2415.925 * emsq) + (1246.113 * eoc);
                    g322 = -342.585 + (1554.908 * em) - (2366.899 * emsq) + (1215.972 * eoc);
                    g410 = -1052.797 + (4758.686 * em) - (7193.992 * emsq) + (3651.957 * eoc);
                    g422 = -3581.690 + (16178.110 * em) - (24462.770 * emsq) + (12422.520 * eoc);
                    g520 = em > 0.715 ? -5149.66 + (29936.92 * em) - (54087.36 * emsq) + (31324.56 * eoc) : 1464.74 - (4664.75 * em) + (3763.64 * emsq);
                }
                if (em < 0.7)
                {
                    g533 = -919.22770 + (4988.6100 * em) - (9064.7700 * emsq) + (5542.21 * eoc);
                    g521 = -822.71072 + (4568.6173 * em) - (8491.4146 * emsq) + (5337.524 * eoc);
                    g532 = -853.66600 + (4690.2500 * em) - (8624.7700 * emsq) + (5341.4 * eoc);
                }
                else
                {
                    g533 = -37995.780 + (161616.52 * em) - (229838.20 * emsq) + (109377.94 * eoc);
                    g521 = -51752.104 + (218913.95 * em) - (309468.16 * emsq) + (146349.42 * eoc);
                    g532 = -40023.880 + (170470.89 * em) - (242699.48 * emsq) + (115605.82 * eoc);
                }
                var sini2 = sinim * sinim;
                var f220 = 0.75 * (1.0 + (2.0 * cosim) + cosisq);
                var f221 = 1.5 * sini2;
                var f321 = 1.875 * sinim * (1.0 - (2.0 * cosim) - (3.0 * cosisq));
                var f322 = -1.875 * sinim * (1.0 + (2.0 * cosim) - (3.0 * cosisq));
                var f441 = 35.0 * sini2 * f220;
                var f442 = 39.3750 * sini2 * sini2;

                var f522 = 9.84375 * sinim * ((sini2 * (1.0 - (2.0 * cosim) - (5.0 * cosisq)))
                     + (0.33333333 * (-2.0 + (4.0 * cosim) + (6.0 * cosisq))));
                var f523 = sinim * ((4.92187512 * sini2 * (-2.0 - (4.0 * cosim) + (10.0 * cosisq)))
                  + (6.56250012 * (1.0 + (2.0 * cosim) - (3.0 * cosisq))));
                var f542 = 29.53125 * sinim * (2.0 - (8.0 * cosim)
                  + (cosisq * (-12.0 + (8.0 * cosim) + (10.0 * cosisq))));
                var f543 = 29.53125 * sinim * (-2.0 - (8.0 * cosim)
                  + (cosisq * (12.0 + (8.0 * cosim) - (10.0 * cosisq))));

                var xno2 = nm * nm;
                var ainv2 = aonv * aonv;
                var temp1 = 3.0 * xno2 * ainv2;
                var temp = temp1 * root22;
                d2201 = temp * f220 * g201;
                d2211 = temp * f221 * g211;
                temp1 *= aonv;
                temp = temp1 * root32;
                d3210 = temp * f321 * g310;
                d3222 = temp * f322 * g322;
                temp1 *= aonv;
                temp = 2.0 * temp1 * root44;
                d4410 = temp * f441 * g410;
                d4422 = temp * f442 * g422;
                temp1 *= aonv;
                temp = temp1 * root52;
                d5220 = temp * f522 * g520;
                d5232 = temp * f523 * g532;
                temp = 2.0 * temp1 * root54;
                d5421 = temp * f542 * g521;
                d5433 = temp * f543 * g533;
                xlamo = (mo + nodeo + nodeo - (theta + theta)) % twoPi;
                xfact = mdot + dmdt + (2.0 * (nodedot + dnodt - rptim)) - no;
                em = emo;
                emsq = emsqo;
            }

            //  ---------------- synchronous resonance terms --------------
            if (irez == 1)
            {
                var g200 = 1.0 + (emsq * (-2.5 + (0.8125 * emsq)));
                var g310 = 1.0 + (2.0 * emsq);
                var g300 = 1.0 + (emsq * (-6.0 + (6.60937 * emsq)));
                var f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
                var f311 = (0.9375 * sinim * sinim * (1.0 + (3.0 * cosim))) - (0.75 * (1.0 + cosim));
                var f330 = 1.0 + cosim;

                f330 *= 1.875 * f330 * f330;
                del1 = 3.0 * nm * nm * aonv * aonv;
                del2 = 2.0 * del1 * f220 * g200 * q22;
                del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
                del1 = del1 * f311 * g310 * q31 * aonv;
                xlamo = (mo + nodeo + argpo - theta) % twoPi;
                xfact = mdot + xpidot + dmdt + domdt + dnodt - (no + rptim);
            }

            //  ------------ for sgp4, initialize the integrator ----------
            xli = xlamo;
            xni = no;
            atime = 0.0;
            nm = no + dndt;
        }

        return new DsInitResult(
            em: em,
      argpm: argpm,
      inclm: inclm,
      mm: mm,
      nm: nm,
      nodem: nodem,

      irez: irez,
      atime: atime,

      d2201: d2201,
      d2211: d2211,
      d3210: d3210,
      d3222: d3222,
      d4410: d4410,

      d4422: d4422,
      d5220: d5220,
      d5232: d5232,
      d5421: d5421,
      d5433: d5433,

      dedt: dedt,
      didt: didt,
      dmdt: dmdt,
      dndt: dndt,
      dnodt: dnodt,
      domdt: domdt,

      del1: del1,
      del2: del2,
      del3: del3,

      xfact: xfact,
      xlamo: xlamo,
      xli: xli,
      xni: xni
    );
    }


    internal class DsInitResult
    {
        public DsInitResult(int irez, double em, double argpm, double inclm, double mm, double nm, double nodem, double atime, double d2201, double d2211, double d3210, double d3222, double d4410, double d4422, double d5220, double d5232, double d5421, double d5433, double dedt, double didt, double dmdt, double dndt, double dnodt, double domdt, double del1, double del2, double del3, double xfact, double xlamo, double xli, double xni)
        {
            this.em = em;
            this.argpm = argpm;
            this.inclm = inclm;
            this.mm = mm;
            this.nm = nm;
            this.nodem = nodem;

            this.irez = irez;
            this.atime = atime;

            this.d2201 = d2201;
            this.d2211 = d2211;
            this.d3210 = d3210;
            this.d3222 = d3222;
            this.d4410 = d4410;

            this.d4422 = d4422;
            this.d5220 = d5220;
            this.d5232 = d5232;
            this.d5421 = d5421;
            this.d5433 = d5433;

            this.dedt = dedt;
            this.didt = didt;
            this.dmdt = dmdt;
            this.dndt = dndt;
            this.dnodt = dnodt;
            this.domdt = domdt;

            this.del1 = del1;
            this.del2 = del2;
            this.del3 = del3;

            this.xfact = xfact;
            this.xlamo = xlamo;
            this.xli = xli;
            this.xni = xni;
        }


        public readonly int irez;
        public readonly double em;
        public readonly double argpm;
        public readonly double inclm;
        public readonly double mm;
        public readonly double nm;
        public readonly double nodem;
        public readonly double atime;
        public readonly double d2201;
        public readonly double d2211;
        public readonly double d3210;
        public readonly double d3222;
        public readonly double d4410;
        public readonly double d4422;
        public readonly double d5220;
        public readonly double d5232;
        public readonly double d5421;
        public readonly double d5433;
        public readonly double dedt;
        public readonly double didt;
        public readonly double dmdt;
        public readonly double dndt;
        public readonly double dnodt;
        public readonly double domdt;
        public readonly double del1;
        public readonly double del2;
        public readonly double del3;
        public readonly double xfact;
        public readonly double xlamo;
        public readonly double xli;
        public readonly double xni;
    }

    /*-----------------------------------------------------------------------------
    *
    *                           procedure dspace
    *
    *  this procedure provides deep space contributions to mean elements for
    *    perturbing third body.  these effects have been averaged over one
    *    revolution of the sun and moon.  for earth resonance effects, the
    *    effects have been averaged over no revolutions of the satellite.
    *    (mean motion)
    *
    *  author        : david vallado                  719-573-2600   28 jun 2005
    *
    *  inputs        :
    *    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
    *    dedt        -
    *    del1, del2, del3  -
    *    didt        -
    *    dmdt        -
    *    dnodt       -
    *    domdt       -
    *    irez        - flag for resonance           0-none, 1-one day, 2-half day
    *    argpo       - argument of perigee
    *    argpdot     - argument of perigee dot (rate)
    *    t           - time
    *    tc          -
    *    gsto        - gst
    *    xfact       -
    *    xlamo       -
    *    no          - mean motion
    *    atime       -
    *    em          - eccentricity
    *    ft          -
    *    argpm       - argument of perigee
    *    inclm       - inclination
    *    xli         -
    *    mm          - mean anomaly
    *    xni         - mean motion
    *    nodem       - right ascension of ascending node
    *
    *  outputs       :
    *    atime       -
    *    em          - eccentricity
    *    argpm       - argument of perigee
    *    inclm       - inclination
    *    xli         -
    *    mm          - mean anomaly
    *    xni         -
    *    nodem       - right ascension of ascending node
    *    dndt        -
    *    nm          - mean motion
    *
    *  locals        :
    *    delt        -
    *    ft          -
    *    theta       -
    *    x2li        -
    *    x2omi       -
    *    xl          -
    *    xldot       -
    *    xnddt       -
    *    xndt        -
    *    xomi        -
    *
    *  coupling      :
    *    none        -
    *
    *  references    :
    *    hoots, roehrich, norad spacetrack report #3 1980
    *    hoots, norad spacetrack report #6 1986
    *    hoots, schumacher and glover 2004
    *    vallado, crawford, hujsak, kelso  2006
    ----------------------------------------------------------------------------*/
    private static DspaceResult Dspace(
        SGP4 satrec,
            double t, double tc,
             //
             double em, double argpm, double inclm, double xli, double mm, double xni, double nodem, double nm
        )
    {
        var irez = satrec.irez;
        var d2201 = satrec.d2201;
        var d2211 = satrec.d2211;
        var d3210 = satrec.d3210;
        var d3222 = satrec.d3222;
        var d4410 = satrec.d4410;
        var d4422 = satrec.d4422;
        var d5220 = satrec.d5220;
        var d5232 = satrec.d5232;
        var d5421 = satrec.d5421;
        var d5433 = satrec.d5433;
        var dedt = satrec.dedt;
        var del1 = satrec.del1;
        var del2 = satrec.del2;
        var del3 = satrec.del3;
        var didt = satrec.didt;
        var dmdt = satrec.dmdt;
        var dnodt = satrec.dnodt;
        var domdt = satrec.domdt;
        var argpo = satrec.argpo;
        var argpdot = satrec.argpdot;

        var gsto = satrec.gsto;
        var xfact = satrec.xfact;
        var xlamo = satrec.xlamo;
        var no = satrec.no;
        var atime = satrec.atime;


        const double fasx2 = 0.13130908;
        const double fasx4 = 2.8843198;
        const double fasx6 = 0.37448087;
        const double g22 = 5.7686396;
        const double g32 = 0.95240898;
        const double g44 = 1.8014998;
        const double g52 = 1.0508330;
        const double g54 = 4.4108898;
        const double rptim = 4.37526908801129966e-3; // equates to 7.29211514668855e-5 rad/sec
        const double stepp = 720.0;
        const double stepn = -720.0;
        const double step2 = 259200.0;

        double x2li;
        double x2omi;
        double xl;
        double xldot = 0;
        double xnddt = 0;
        double xndt = 0;
        double xomi;
        double dndt = 0.0;
        double ft = 0.0;

        //  ----------- calculate deep space resonance effects -----------
        double theta = (gsto + (tc * rptim)) % twoPi;
        em += dedt * t;
        inclm += didt * t;
        argpm += domdt * t;
        nodem += dnodt * t;
        mm += dmdt * t;

        // sgp4fix for negative inclinations
        // the following if statement should be commented out
        // if (inclm < 0.0)
        // {
        //   inclm = -inclm;
        //   argpm = argpm - pi;
        //   nodem = nodem + pi;
        // }

        /* - update resonances : numerical (euler-maclaurin) integration - */
        /* ------------------------- epoch restart ----------------------  */
        //   sgp4fix for propagator problems
        //   the following integration works for negative time steps and periods
        //   the specific changes are unknown because the original code was so convoluted

        // sgp4fix take out atime = 0.0 and fix for faster operation

        if (irez != 0)
        {
            //  sgp4fix streamline check
            if (atime == 0.0 || t * atime <= 0.0 || Math.Abs(t) < Math.Abs(atime))
            {
                atime = 0.0;
                xni = no;
                xli = xlamo;
            }

            // sgp4fix move check outside loop
            var delt = t > 0.0 ? stepp : stepn;

            var iretn = 381; // added for do loop
            while (iretn == 381)
            {
                //  ------------------- dot terms calculated -------------
                //  ----------- near - synchronous resonance terms -------
                if (irez != 2)
                {
                    xndt = (del1 * Math.Sin(xli - fasx2))
                      + (del2 * Math.Sin(2.0 * (xli - fasx4)))
                      + (del3 * Math.Sin(3.0 * (xli - fasx6)));
                    xldot = xni + xfact;
                    xnddt = (del1 * Math.Cos(xli - fasx2))
                                          + (2.0 * del2 * Math.Cos(2.0 * (xli - fasx4)))
                                          + (3.0 * del3 * Math.Cos(3.0 * (xli - fasx6)));
                    xnddt *= xldot;
                }
                else
                {
                    // --------- near - half-day resonance terms --------
                    xomi = argpo + (argpdot * atime);
                    x2omi = xomi + xomi;
                    x2li = xli + xli;
                    xndt = (d2201 * Math.Sin(x2omi + xli - g22))
                      + (d2211 * Math.Sin(xli - g22))
                      + (d3210 * Math.Sin(xomi + xli - g32))
                      + (d3222 * Math.Sin(-xomi + xli - g32))
                      + (d4410 * Math.Sin(x2omi + x2li - g44))
                      + (d4422 * Math.Sin(x2li - g44))
                      + (d5220 * Math.Sin(xomi + xli - g52))
                      + (d5232 * Math.Sin(-xomi + xli - g52))
                      + (d5421 * Math.Sin(xomi + x2li - g54))
                      + (d5433 * Math.Sin(-xomi + x2li - g54));
                    xldot = xni + xfact;
                    xnddt = (d2201 * Math.Cos(x2omi + xli - g22))
                      + (d2211 * Math.Cos(xli - g22))
                      + (d3210 * Math.Cos(xomi + xli - g32))
                      + (d3222 * Math.Cos(-xomi + xli - g32))
                      + (d5220 * Math.Cos(xomi + xli - g52))
                      + (d5232 * Math.Cos(-xomi + xli - g52))
                      + (2.0 * d4410 * Math.Cos(x2omi + x2li - g44))
                      + (d4422 * Math.Cos(x2li - g44))
                      + (d5421 * Math.Cos(xomi + x2li - g54))
                      + (d5433 * Math.Cos(-xomi + x2li - g54));
                    xnddt *= xldot;
                }

                //  ----------------------- integrator -------------------
                //  sgp4fix move end checks to end of routine
                if (Math.Abs(t - atime) >= stepp)
                {
                    iretn = 381;
                }
                else
                {
                    ft = t - atime;
                    iretn = 0;
                }

                if (iretn == 381)
                {
                    xli += (xldot * delt) + (xndt * step2);
                    xni += (xndt * delt) + (xnddt * step2);
                    atime += delt;
                }
            }

            nm = xni + (xndt * ft) + (xnddt * ft * ft * 0.5);
            xl = xli + (xldot * ft) + (xndt * ft * ft * 0.5);
            if (irez != 1)
            {
                mm = xl - (2.0 * nodem) + (2.0 * theta);
                dndt = nm - no;
            }
            else
            {
                mm = xl - nodem - argpm + theta;
                dndt = nm - no;
            }
            nm = no + dndt;
        }

        var ret = new DspaceResult(atime: atime, em: em, argpm: argpm, inclm: inclm, xli: xli, mm: mm, xni: xni, nodem: nodem, dndt: dndt, nm: nm);

        return ret;
    }

    private class DspaceResult
    {
        public DspaceResult(double atime, double em, double argpm, double inclm, double xli, double mm, double xni, double nodem, double nm, double dndt)
        {
            this.atime = atime;
            this.em = em;
            this.argpm = argpm;
            this.inclm = inclm;
            this.xli = xli;
            this.mm = mm;
            this.xni = xni;
            this.nodem = nodem;
            this.nm = nm;
            this.dndt = dndt;
        }


        public readonly double atime;
        public readonly double em;
        public readonly double argpm;
        public readonly double inclm;
        public readonly double xli;
        public readonly double mm;
        public readonly double xni;
        public readonly double nodem;
        public readonly double nm;
        public readonly double dndt;
    }
}


