package orc.consoperators;

import org.apache.commons.math.MathException;
import org.apache.commons.math3.util.FastMath;

import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.LogNormalDistributionModel;


public class ConsOperatorUtils {
    final private static double SQRT2 = Math.sqrt(2.0);
    final static private double EPSILON = 1e-8;

    public static double getHRForLN (double rNew, double qOld, ParametricDistribution distribution) {
        double stdev = ((LogNormalDistributionModel) distribution).SParameterInput.get().getArrayValue();
        double miu = - 0.5 * stdev * stdev;

        double b = FastMath.log(rNew);
        double c = 2 * stdev * stdev;
        double x = b - miu;
        double x_sq = x * x / c;
        double rateHR = b + x_sq;

        double a = erfInv(2 * qOld - 1);
        double quantileHR = miu + SQRT2 * stdev * a + a * a;
        return quantileHR - rateHR;
    }

    public static double getHRForPieceWise (double rNew, double qOld, double qNew, ParametricDistribution distribution) {
        PiecewiseLinearDistribution pld = (PiecewiseLinearDistribution) distribution;
        double logHR = Math.log(pld.getDerivativeAtQuantile(qOld));
        logHR +=  Math.log(pld.getDerivativeAtQuantileInverse(rNew, qNew));
        return logHR;
    }

    public static double getHRUseNumericApproximation (double rNew, double qOld, ParametricDistribution distribution) {
        double logHR = 0;
        try {
            double r0 = distribution.inverseCumulativeProbability(qOld);
            double r0h = distribution.inverseCumulativeProbability(qOld + EPSILON);
            logHR += FastMath.log((r0h - r0) / EPSILON);

            double q0 = distribution.cumulativeProbability(rNew);
            double q0h = distribution.cumulativeProbability(rNew + EPSILON);
            if (q0h != q0) {
            	logHR += FastMath.log((q0h - q0) / EPSILON);
            } else {
            	logHR = Double.POSITIVE_INFINITY;
            }
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability!");
        }
        return -logHR;
    }

    public static double erfInv(final double x) {
        double w = - FastMath.log((1.0 - x) * (1.0 + x));
        double p;

        if (w < 6.25) {
            w -= 3.125;
            p = -3.6444120640178196996e-21;
            p = -1.685059138182016589e-19 + p * w;
            p = 1.2858480715256400167e-18 + p * w;
            p = 1.115787767802518096e-17 + p * w;
            p = -1.333171662854620906e-16 + p * w;
            p = 2.0972767875968561637e-17 + p * w;
            p = 6.6376381343583238325e-15 + p * w;
            p = -4.0545662729752068639e-14 + p * w;
            p = -8.1519341976054721522e-14 + p * w;
            p = 2.6335093153082322977e-12 + p * w;
            p = -1.2975133253453532498e-11 + p * w;
            p = -5.4154120542946279317e-11 + p * w;
            p = 1.051212273321532285e-09 + p * w;
            p = -4.1126339803469836976e-09 + p * w;
            p = -2.9070369957882005086e-08 + p * w;
            p = 4.2347877827932403518e-07 + p * w;
            p = -1.3654692000834678645e-06 + p * w;
            p = -1.3882523362786468719e-05 + p * w;
            p = 0.0001867342080340571352 + p * w;
            p = -0.00074070253416626697512 + p * w;
            p = -0.0060336708714301490533 + p * w;
            p = 0.24015818242558961693 + p * w;
            p = 1.6536545626831027356 + p * w;
        } else if (w < 16.0) {
            w = FastMath.sqrt(w) - 3.25;
            p = 2.2137376921775787049e-09;
            p = 9.0756561938885390979e-08 + p * w;
            p = -2.7517406297064545428e-07 + p * w;
            p = 1.8239629214389227755e-08 + p * w;
            p = 1.5027403968909827627e-06 + p * w;
            p = -4.013867526981545969e-06 + p * w;
            p = 2.9234449089955446044e-06 + p * w;
            p = 1.2475304481671778723e-05 + p * w;
            p = -4.7318229009055733981e-05 + p * w;
            p = 6.8284851459573175448e-05 + p * w;
            p = 2.4031110387097893999e-05 + p * w;
            p = -0.0003550375203628474796 + p * w;
            p = 0.00095328937973738049703 + p * w;
            p = -0.0016882755560235047313 + p * w;
            p = 0.0024914420961078508066 + p * w;
            p = -0.0037512085075692412107 + p * w;
            p = 0.005370914553590063617 + p * w;
            p = 1.0052589676941592334 + p * w;
            p = 3.0838856104922207635 + p * w;
        } else if (!Double.isInfinite(w)) {
            w = FastMath.sqrt(w) - 5.0;
            p = -2.7109920616438573243e-11;
            p = -2.5556418169965252055e-10 + p * w;
            p = 1.5076572693500548083e-09 + p * w;
            p = -3.7894654401267369937e-09 + p * w;
            p = 7.6157012080783393804e-09 + p * w;
            p = -1.4960026627149240478e-08 + p * w;
            p = 2.9147953450901080826e-08 + p * w;
            p = -6.7711997758452339498e-08 + p * w;
            p = 2.2900482228026654717e-07 + p * w;
            p = -9.9298272942317002539e-07 + p * w;
            p = 4.5260625972231537039e-06 + p * w;
            p = -1.9681778105531670567e-05 + p * w;
            p = 7.5995277030017761139e-05 + p * w;
            p = -0.00021503011930044477347 + p * w;
            p = -0.00013871931833623122026 + p * w;
            p = 1.0103004648645343977 + p * w;
            p = 4.8499064014085844221 + p * w;
        } else {
            p = Double.POSITIVE_INFINITY;
        }
        return p * x;
    }

    // To test the hastings ratio under lognormal distribution
    public static double calculateHastingsRatio(double r, double q, double stdev) {
        double miu = - 0.5 * stdev * stdev; // miu of lognormal
        double a = ConsOperatorUtils.erfInv(2 * q - 1);
        double b = FastMath.log(r);
        double c = 2 * stdev * stdev;
        double x = b - miu;
        double x_sq = x * x / c;
        double d = Math.sqrt(c);
        return -b - x_sq + miu + (d * a) + (a * a);
    }
}
