package orc.consoperators;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;





@Description("Approximates parametric distribution by piecewise linear approximation.")
public class PiecewiseLinearDistribution extends ParametricDistribution {
    final public Input<ParametricDistribution> distrInput = new Input<>("distr", "Underlying parametric distribution that is approximated.", Validate.REQUIRED);
    final public Input<Integer> numberOfDiscreteRates = new Input<>("bins", "number of bins used to approximate the distribution piecewise linearly. (default 100)", 100);
    
    final public Input<Double> limitInput = new Input<>("limit", "fraction of bins at end of distribution to be cut off -- should be between 0 and 1", 0.1);
    final public Input<Boolean> cutOffEndInput = new Input<>("cutOffEnd", "approximate values below for quantiles close to 0 and 1 with rate at limit/bins and (1-limit/bins) respectively."
    		+ "If false, for quantiles below limit/bins and above (1-limit/bins) the inverseCumulativeProbability methods of the underlying parametric distribution is called"
    		+ "(which is slower, but should not happen very often).", true);

    LinearPiecewiseImpl dist = new LinearPiecewiseImpl();
    ContinuousDistribution underlyingDistr = new NormalDistributionImpl();
    
    // ParametricDistribution distribution;
    
    protected double[] rates; //the output rates
    protected double[] storedRates; //
    
    protected double limit_0 = 0.1; // limit_0 / (numberOfDiscreteRates-1) is the minimum quantile
    protected double limitLow, limitUp;
    protected boolean cutOffEnd = true;
    protected int rateCountMin1;
    protected double intervalSize;

    Map<Double, Double> cache = new HashMap<>(); 
    Map<Double, Double> storedcache = new HashMap<>(); 

    public PiecewiseLinearDistribution() {}
    
    /** construct log-normal distribution with mean = 1 **/ 
    public PiecewiseLinearDistribution(RealParameter S) {
    	LogNormalDistributionModel lnormal = new LogNormalDistributionModel();
    	lnormal.initByName("S", S, "M", "1.0", "meanInRealSpace", true);
    	initByName("distr", lnormal);
    }
    
    
    @Override
	public void initAndValidate() {
    	int dim = numberOfDiscreteRates.get();
    	rates = new double[dim];
    	storedRates = new double[dim];
    	rateCountMin1 = rates.length-1;
    	intervalSize = 1.0 / rateCountMin1;

    	ParametricDistribution distribution = distrInput.get();
    	org.apache.commons.math.distribution.Distribution d = distribution.getDistribution();
    	if (d instanceof ContinuousDistribution) {
    		underlyingDistr = (ContinuousDistribution) d;
    	} else {
    		throw new IllegalArgumentException("Expected parametric distribution that is continuous");
    	}
    	limit_0 = limitInput.get();
    	if (limit_0 <= 0 || limit_0 >= 1) {
    		throw new IllegalArgumentException("limit should be between 0 and 1");
    	}
    	cutOffEnd = cutOffEndInput.get();
    	if (cutOffEnd) {
	    	limitLow = limit_0 / dim;
	    	limitUp = 1.0 - limit_0 / dim;
    	} else {
    		// mark the lowest and highest bin for
    		// using the underlying distribution
	    	limitLow = 1.0 / dim;
	    	limitUp = 1.0 - 1.0 / dim;
    	}
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    void refresh() {
    }

    @Override
    public org.apache.commons.math.distribution.Distribution getDistribution() {
        refresh();
        return dist;
    }

    public class LinearPiecewiseImpl implements ContinuousDistribution {

        public LinearPiecewiseImpl() {
        }


        @Override
        public double cumulativeProbability(double x) throws MathException {
        	underlyingDistr = getUnderlyingDistr();
        	// Return exact cdf using piecewise linear approximation
            int i = getIntervalFor(x);

            if (!cutOffEnd) {
            	// TODO: needs testing
            	if (i == 0) {
//            		if (x < rates[0]) {
                		double cdf = underlyingDistr.cumulativeProbability(x);
                		cdf = Math.max(cdf,  1e-256);
                		return cdf;
//            		}
//                	double cdf = (limitLow + (1-limit_0) * (x - rates[i]) / (rates[i+1] - rates[i])) * intervalSize;
//            		return cdf;
            	} else if (i >= rateCountMin1 - 1) {
            		double cdf = underlyingDistr.cumulativeProbability(x);
            		cdf = Math.min(cdf,  0.9999999999999999);
            		return cdf;
//            	} else if (i == rateCountMin1-1) {
//                	double cdf = (i + (1-limit_0) * (x - rates[i]) / (rates[i+1] - rates[i])) * intervalSize;
//            		return cdf;
            	}
            }
            
            if (i < rateCountMin1) {
            	double cdf = (i + (x - rates[i]) / (rates[i+1] - rates[i])) * intervalSize;
            	cdf = Math.max(0, cdf);
            	return cdf;
            }
            return 1.0;
        }
        
		@Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public double inverseCumulativeProbability(double q) throws MathException {
        	underlyingDistr = getUnderlyingDistr();
        	if (!cutOffEnd && (q <= limitLow  || q >= limitUp)) {        		
            	final Double qD = q;
            	if (cache.containsKey(qD)) {
            		return cache.get(qD);
            	}
            	double x = underlyingDistr.inverseCumulativeProbability(q);
            	cache.put(qD, x);
        		return x;
        	}
            double v = q * rateCountMin1;
            int i = (int) v;
            
            // make sure cached rates are calculated
            if (rates[i] == 0.0) {
    	        try {
    	        	if (i > 0) {
    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) * intervalSize);
    	        	} else {
    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 * intervalSize);
    	        	}
    	        } catch (MathException e) {
    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
    	        }
            }
            if (i < rateCountMin1 && rates[i + 1] == 0.0) {
    	        try {
    	        	if (i < rateCountMin1 - 1) {
    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) * intervalSize);
    	        	} else {
    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rateCountMin1 - limit_0) * intervalSize);
    	        	}
    	        } catch (MathException e) {
    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
    	        }
            }


            sanitycheck();
            
            if (cutOffEnd) {
            	if (i == 0) {
            		v = (v-limit_0) / (1-limit_0);
            		if (v < 0) {
            			v = 0;
            		}
            	} else if (i == rateCountMin1) {
            		v = i + (v - i) / (1-limit_0);
            	}
            }
            // return piecewise linear approximation
            double r = rates[i];
            if (i < rateCountMin1) {
            	r += (rates[i+1] - rates[i]) * (v - i);
            }
            return r;
        }

        @Override
        public double density(double x) {
        	throw new IllegalArgumentException("not implemented yet");
        }

        @Override
        public double logDensity(double x) {
        	return Math.log(density(x));
        }


        protected int getIntervalFor(double x) throws MathException {
            double q = underlyingDistr.cumulativeProbability(x);
            int i = getIntervalFor(x, q);
			return i;
		}

        protected int getIntervalFor(double r, double qNew) {
	        double v = qNew * rateCountMin1;
	        int i = (int) v;
	        if (rates[i] == 0.0) {
		        try {
		        	if (i > 0) {
		        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) * intervalSize);
		        	} else {
		        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 * intervalSize);
		        	}
		        } catch (MathException e) {
		            throw new RuntimeException("Failed to compute inverse cumulative probability!");
		        }
	        }
	        if (i < rateCountMin1 && rates[i + 1] == 0.0) {
		        try {
		        	if (i < rateCountMin1 - 1) {
		        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) * intervalSize);
		        	} else {
		        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rateCountMin1 - limit_0) * intervalSize);
		        	}
		        } catch (MathException e) {
		            throw new RuntimeException("Failed to compute inverse cumulative probability!");
		        }
	        }
	        
	        // test boundary: r should be between rates[i] and rates[i+1]
	        // but due to numerical errors in the piecewise linear approximation could
	        // fall just outside, so make sure the condition is met, and rates are calculated
	        while (i > 0 && r < rates[i]) {
	        	i--;
	            if (i > 0 && rates[i] == 0.0) {
	    	        try {
	    	        	if (i > 0) {
	    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) * intervalSize);
	    	        	} else {
	    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 * intervalSize);
	    	        	}
	    	        } catch (MathException e) {
	    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	    	        }
	            }
	        }
	        while (i < rateCountMin1 && r > rates[i+1]) {
	        	i++;
	            if (i < rateCountMin1 && rates[i + 1] == 0.0) {
	    	        try {
	    	        	if (i < rateCountMin1 - 1) {
	    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) * intervalSize);
	    	        	} else {
	    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rateCountMin1 - limit_0) * intervalSize);
	    	        	}
	    	        } catch (MathException e) {
	    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	    	        }
	            }
	        }
            sanitycheck();
	        return i;
		}
    } // class LinearPiecewiseImpl

    @Override
    protected double getMeanWithoutOffset() {
    	throw new IllegalArgumentException("not implemented yet");
    }
    
    @Override
    protected void store() {
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        underlyingDistr = getUnderlyingDistr();
    	storedcache.clear();
    	for (Double d : cache.values()) {
    		storedcache.put(d, cache.get(d));
    	}
//System.err.println("PLD  store" + ((BEASTInterface)underlyingDistr).getInput("mode").get() + " " + Arrays.toString(storedRates));
        super.store();

        sanitycheck();
    }

    private void sanitycheck() {
//    	if (underlyingDistr != getUnderlyingDistr()) {
//			int h = 43;
// 			h--;    		
//    	}
//
//        for (int i = 1; i < rates.length; i++) {
//			if (rates[i-1] > 0) {
//				int j = i;
//				while (j < rates.length && rates[j] == 0) {
//					j++;
//				}
//				if (j < rates.length && rates[i-1] > rates[j]) {
//					int h = 43;
//		 			h--;
//				}
//			}
//        }
    }
    
    @Override
    protected void restore() {
        double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;

        Map<Double,Double> tmp2 = storedcache;
    	storedcache = cache;
    	cache = tmp2;
        
        underlyingDistr = getUnderlyingDistr();
        
        // System.err.println("PLDrestore " + ((BEASTInterface)underlyingDistr).getInput("mode").get() + " " + Arrays.toString(rates));
    	super.restore();
    }

    protected ContinuousDistribution getUnderlyingDistr() {
    	return (ContinuousDistribution) distrInput.get().getDistribution();
	}

	@Override
    protected boolean requiresRecalculation() {
    	
    	if (distrInput.get() != null && distrInput.get().isDirtyCalculation()) {
    		Arrays.fill(rates, 0.0);
    		cache.clear();
    		underlyingDistr = getUnderlyingDistr();
    		return true;
    	}
    	
    	return super.requiresRecalculation();
    }
    
	/**
	 * Assumes quantile parameterisation of clockModel, so clockModel must be specified
	 * @param q quantile
	 * @return derivative of rate distribution at quantile q
	 */
	final static private double EPSILON = 1e-12;
	public double getDerivativeAtQuantile(double q) {
    	if (!cutOffEnd && (q < limitLow  || q > limitUp)) {
			try {
				if (q + EPSILON >= 1.0) {
					return Double.POSITIVE_INFINITY;
				}
	    		double r = underlyingDistr.inverseCumulativeProbability(q);
	    		double rPlusH = underlyingDistr.inverseCumulativeProbability(q + EPSILON);
	    		double dR = (rPlusH - r) / EPSILON;
	    		return dR;
			} catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!" + e.getMessage());
			}
    	}

		// use cached rates
        double v = q * rateCountMin1;
        int i = (int) v;
        if (rates[i] == 0.0) {
	        try {
	        	if (i > 0) {
	        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) * intervalSize);
	        	} else {
	        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 * intervalSize);
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
        }
        if (i < rateCountMin1 && rates[i + 1] == 0.0) {
	        try {
	        	if (i < rateCountMin1 - 1) {
	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) * intervalSize);
	        	} else {
	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rateCountMin1 - limit_0) * intervalSize);
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
        }
        if (!cutOffEnd) {
        	if (i == 0 || i == rateCountMin1 - 1) {
            	double r = (rates[i+1] - rates[i]) / ((1.0 - limit_0) * intervalSize);
            	return r;        		
            }
        }
        if (i < rateCountMin1) {
        	double r = (rates[i+1] - rates[i]) / intervalSize;
        	return r;
        }
        return 0;
    }
    
	/**
	 * Assumes quantile parameterisation of clockModel, so clockModel must be specified
	 * @param r rate
	 * @param qNew: quantile associated with rate r
	 * @return derivative of quantile distribution at rate r
	 */
    public double getDerivativeAtQuantileInverse(double r, double qNew) {
    	// TODO: verify the following is correct
    	return 1.0 / getDerivativeAtQuantile(qNew);
    	
//    	if (!cutOffEnd && (qNew < limitLow  || qNew > limitUp)) {
//			try {
//				double r0 = underlyingDistr.inverseCumulativeProbability(qNew);
//	    		double rPlusH = underlyingDistr.inverseCumulativeProbability(qNew + EPSILON);
//	    		double dR = 1.0/((rPlusH - r0) / EPSILON);
//	    		return dR;
//			} catch (MathException e) {
//				e.printStackTrace();
//			}
//    	}
//
//    	int i = dist.getIntervalFor(r, qNew);
//        if (i < rateCountMin1) {
//            double derivative = (1.0/rateCountMin1)/(rates[i+1] - rates[i]);
//            return derivative;        	
//        }
//        return 0;
    }
    
    
    
    /**
     * @return The minimum rate supported by the function 
     * @throws MathException 
     */
    public double getRangeMin() throws MathException {
    	if (!cutOffEnd) {
    		return underlyingDistr.inverseCumulativeProbability(0.0);
    	}
    	if (rates[0] == 0.0) {
    		rates[0] = underlyingDistr.inverseCumulativeProbability(limit_0 * intervalSize);
    	}
    	return rates[0];
    	//return getDerivativeAtQuantile(limit_0 / rateCountMin1);
    }
    
    
    /**
     * @return The maximum rate supported by the function 
     * @throws MathException 
     */
    public double getRangeMax() throws MathException {
    	if (!cutOffEnd) {
    		return underlyingDistr.inverseCumulativeProbability(1.0);
    	}
    	if (rates[rateCountMin1] == 0.0) {
    		rates[rateCountMin1] = underlyingDistr.inverseCumulativeProbability((rateCountMin1 - limit_0) * intervalSize);
    	}
    	return rates[rateCountMin1];
    	//return getDerivativeAtQuantile((rateCountMin1 - limit_0) / rateCountMin1);
    }
    
    
    
    
 }
