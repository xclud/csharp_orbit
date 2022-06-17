using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace System;

internal sealed class DeepSpaceLongPeriodPeriodicContributions
{
    public DeepSpaceLongPeriodPeriodicContributions(double ep, double inclp, double nodep, double argpp, double mp)
    {
        this.ep = ep;
        this.inclp = inclp;
        this.argpp = argpp;
        this.mp = mp;
        this.nodep = nodep;
    }

    public readonly double ep;
    public readonly double inclp;
    public readonly double nodep;
    public readonly double argpp;
    public readonly double mp;
}
