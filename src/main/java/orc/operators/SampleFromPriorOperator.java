package orc.operators;

import java.util.ArrayList;
import java.util.List;



import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.spec.inference.distribution.IID;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.Tensor;
import beast.base.core.Log;
import beast.base.util.Randomizer;
import orc.distribution.BranchRatePrior;




/**
 * Samples a parameter from its prior distribution
 * If no prior is specified, then assumes a uniform distribution between parameter min and max
 * If there is no min/max then throws error
 * If there is more than one element in the parameter, then the number of elements sampled 
 * 		per iteration is sample from a Binomial(n = ndimensions, p) distribution where np is tunable
 * @author Jordan Douglas
 *
 */
public class SampleFromPriorOperator extends Operator {
	
	
	
    final public Input<Tensor<?,?>> paramInput = new Input<>("parameter", "the parameter to sample", Input.Validate.REQUIRED);
    final public Input<TensorDistribution<?,?>> priorInput = new Input<>("prior", "the prior distribution of the parameter", Input.Validate.REQUIRED);
    //final public Input<ScalarDistribution<RealScalar<Real>, Double>> prior2Input = new Input<>("prior2", "the prior distribution of the parameter", Validate.XOR, priorInput);
    final public Input<Double> npInput = new Input<>("np", "tunable parameter describing the mean number of elements in the parameter vector to sample", 1.0);
    
    
    Tensor<?,?> parameter;
    TensorDistribution<?,?> prior;
    //ScalarDistribution<RealScalar<Real>, Double> prior;
    double np;

    

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
		parameter = paramInput.get();
		prior = priorInput.get();
		
		// Check parameters
    	if (!(parameter instanceof RealScalarParam || parameter instanceof RealVectorParam)) {
     		throw new IllegalArgumentException("Only Real/Int Scalar/Vector Parameers are allowed as parameter inputs");
    	}
    	
    	
    	if (parameter instanceof RealScalarParam) {
    		if (!(prior instanceof ScalarDistribution<?,?>)) {
    			throw new IllegalArgumentException("Prior must be a Scalar distribution to match the parameter!");
    		}
    	}
    	
    	else if (parameter instanceof RealVectorParam) {
    		if (!(prior instanceof IID) && !(prior instanceof BranchRatePrior)) {
    			throw new IllegalArgumentException("Parameters must have an iid prior");
    		}
    		if (prior instanceof IID) {
	    		IID<?,?,?> iid = (IID<?,?,?>) prior;
	    		if (!(iid.distInput.get() instanceof ScalarDistribution<?,?>)) {
	    			throw new IllegalArgumentException("IID's prior must be a Scalar distribution to match the parameter!");
	    		}
    		}
    	}
		
		
		
		if (prior == null) {
			throw new IllegalArgumentException("prior2");
			//prior = null; // TODO prior2Input.get().distInput.get();
		}
		np = npInput.get();
		this.validateNP(true);
		
		
	}
	
	
	
	private double getParamValue(int index) {
		if (parameter instanceof RealScalarParam) {
			RealScalarParam<?> param = (RealScalarParam<?>) parameter;
			return param.get();
		}else {
			RealVectorParam<?> param = (RealVectorParam<?>) parameter;
			return param.get(index);
		}
	}
	
	private void setParamValue(int index, double value) {
		if (parameter instanceof RealScalarParam) {
			RealScalarParam<?> param = (RealScalarParam<?>) parameter;
			param.set(value);
		}else {
			RealVectorParam<?> param = (RealVectorParam<?>) parameter;
			param.set(index, value);
		}
	}
	
	

	@Override
	public double proposal() {
		
		//System.out.println("---------- prior=" + prior.density(1.0) + " ----------");
		
		double logHR = 0;
		List<Integer> paramsToSample = this.sampleParamsToSample();

		ScalarDistribution<?,?> priorDist = null;
		
		if (parameter instanceof RealScalarParam) {
			RealScalarParam<?> param = (RealScalarParam<?>) parameter;
			priorDist = (ScalarDistribution<?,?>) prior;
		}else if (prior instanceof IID) {
			RealVectorParam<?> param = (RealVectorParam<?>) parameter;
			IID<?,?,?> iid = (IID<?,?,?>) prior;
			priorDist = (ScalarDistribution<?,?>) iid.distInput.get();
		}else {
			RealVectorParam<?> param = (RealVectorParam<?>) parameter;
			BranchRatePrior d = (BranchRatePrior) prior;
			priorDist = d.getDist();
		}
		
		
		
		
		
		
		for (Integer p : paramsToSample) {


			// After proposal
			try {
				
				
				double oldLogP = 0;
				double newLogP = 0;
			
				
				//System.out.println("Sampling parameter " + p);
				
				
				// Before proposal
				double oldX = this.getParamValue(p);
				oldLogP = prior == null ? 0 : priorDist.logDensity(oldX);
				
				
				// Sample x from the prior
				double u = Randomizer.nextFloat();
				double newX = (double) priorDist.inverseCumulativeProbability(u);
				
				// Calculate log-density of new x
				newLogP = prior == null ? 0 : priorDist.logDensity(newX);
				
				// Set new value
				this.setParamValue(p, newX);
				//parameter.set(p, newX);
				
				
				//System.out.println("From " + oldX + " to " + newX);
				
				
				// Hastings ratio
				logHR += oldLogP - newLogP;
				
				
			} catch (Exception e) {
				// TODO Auto-generated catch block
				//e.printStackTrace();
				return Double.NEGATIVE_INFINITY;
			}
			
		}
		

		return logHR;
	}
	
	
	
	/**
	 * @return The elements in 'parameter' to sample (without replacement)
	 */
	private List<Integer> sampleParamsToSample() {
		
		List<Integer> sampledIndices = new ArrayList<Integer>();
		if (parameter instanceof RealScalarParam) {
			sampledIndices.add(0);
		}else {
			RealVectorParam<?> p = (RealVectorParam<?>) parameter;
			 double prob = Math.min(1.0, this.np / this.parameter.size());
			for (int i = 0; i < this.parameter.size(); i ++) {
				boolean toSample = Randomizer.nextFloat() < prob;
				if (toSample) {
					sampledIndices.add(i);
				}
			}
			
			
		}
		
		
	    
	    
	    // Check that at least 1 thing is being changed. If not, then select one uniformly at random
	    if (sampledIndices.size() == 0) {
	    	sampledIndices.add(Randomizer.nextInt(this.parameter.size()));
	    }
	    
    	
		return sampledIndices;
	}
	
	
    @Override
    public double getCoercableParameterValue() {
        return this.np;
    }


    @Override
    public void setCoercableParameterValue(double value) {
    	this.np = value;
    	this.validateNP(false);
    }
    
    
    /**
     * Ensures that np is no larger than n
     * If np is larger than n, then the tunable parameter may wander 
     * off into infinity and never return if needed
     */
    private void validateNP(boolean verbose) {
    	if (this.np > this.parameter.size()) {
    		if (verbose) Log.warning("Setting np to " + this.parameter.size() + " so that it is no larger than n");
    		this.np = this.parameter.size();
    	}
    	if (this.np < 0) {
    		if (verbose) Log.warning("Setting np to 0 so that it is non-negative");
    		this.np = 0;
    	}
    }
    
    @Override
    public void optimize(double logAlpha) {
    	
        double delta = calcDelta(logAlpha);
        delta += Math.log(this.np);
        this.np = Math.exp(delta);
        this.validateNP(false);
        
    }
    
    
    


    
    

}
