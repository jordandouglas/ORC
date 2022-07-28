package orc.consoperators;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.math.distributions.ParametricDistribution;

@Description("Caches parametric distribution for inverscumulative methods.")
public class CachedDistribution extends ParametricDistribution {
    final public Input<ParametricDistribution> distrInput = new Input<>("distr", "Underlying parametric distribution that is cached.", Validate.REQUIRED);
    
    ContinuousDistribution dist;
	ContinuousDistribution underlyingDistr;
    
    Map<Double, Double> cache = new HashMap<>(); 
    Map<Double, Double> storedcache = new HashMap<>(); 

    @Override
	public void initAndValidate() {
    	dist = new CachedImpl();
        underlyingDistr = getUnderlyingDistr();
    }
    
    
    /**
     * make sure internal state is up to date *
     */
    void refresh() {
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }

    public class CachedImpl implements ContinuousDistribution {
    	
        public CachedImpl() {
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
        	// System.err.println(x);
        	if (x <= 0) {
        		return 0;
        	}
        	if (Double.isInfinite(x)) {
        		return 1;
        	}
        	return underlyingDistr.cumulativeProbability(x);
        }
        
		@Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }
		
        @Override
        public double inverseCumulativeProbability(double q) throws MathException {
        	if (q >= 1) {
        		return Double.POSITIVE_INFINITY;
        	}
        	if (Double.isNaN(q)) {
        		return 0;
        	}
        	final Double qD = q;
        	if (cache.containsKey(qD)) {
        		return cache.get(qD);
        	}
        	double x = underlyingDistr.inverseCumulativeProbability(q);
        	cache.put(qD, x);
        	return x;
//        	return underlyingDistr.inverseCumulativeProbability(q);
        }

        @Override
        public double density(double x) {
        	return underlyingDistr.density(x);
        }

        @Override
        public double logDensity(double x) {
        	return Math.log(density(x));
        }
    } // class CachedImpl

    @Override//        	if (cache.containsKey(qD)) {
//	return cache.get(qD);
//}

    protected double getMeanWithoutOffset() {
    	throw new IllegalArgumentException("not implemented yet");
    }
    
    @Override
    protected void store() {
//		cache.clear();
    	storedcache.clear();
    	for (Double d : cache.values()) {
    		storedcache.put(d, cache.get(d));
    	}
        super.store();
        underlyingDistr = getUnderlyingDistr();
    }
    
    @Override
    protected void restore() {
//		cache.clear();
    	Map<Double,Double> tmp = storedcache;
    	storedcache = cache;
    	cache = tmp;
        underlyingDistr = getUnderlyingDistr();
    	super.restore();
    }

    protected ContinuousDistribution getUnderlyingDistr() {
    	return (ContinuousDistribution) distrInput.get().getDistribution();
	}

	@Override
    protected boolean requiresRecalculation() {

		
    	if (distrInput.get().isDirtyCalculation()) {
    		cache.clear();
    		underlyingDistr = getUnderlyingDistr();
    		return true;
    	}
    	
    	return super.requiresRecalculation();
    }
    
 }
