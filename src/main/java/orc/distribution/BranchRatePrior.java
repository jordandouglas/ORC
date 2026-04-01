package orc.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;


import org.apache.commons.statistics.distribution.LogNormalDistribution;


@Description("Log normal prior distribution on branch rates, with a mean branch rate of 1")
public class BranchRatePrior extends TensorDistribution<RealVector<NonNegativeReal>, Double>  {
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree that contains the branch rates", Input.Validate.REQUIRED); 
	final public Input<RealScalar<PositiveReal>> sigmaInput = new Input<>("sigma", "relaxed clock standard deviation (log normal).", Input.Validate.REQUIRED); 
	
	
	@Override
	public double calculateLogP() {
	
		RealVector<NonNegativeReal> branchRates = paramInput.get();
		logP = calcLogP(branchRates.getElements());
		return logP;
	
	}



	@Override
	public void refresh() {
		// TODO Auto-generated method stub
		
	}



	@Override
	public Double getLowerBoundOfParameter() {
		return 0.0;
	}



	@Override
	public Double getUpperBoundOfParameter() {
		// TODO Auto-generated method stub
		return Double.POSITIVE_INFINITY;
	}



	@Override
	protected double calcLogP(Double... value) {
		return this.calcLogP(Arrays.asList(value));
	}
	
	
    private double calcLogP(List<Double> branchRates) {
    	
    	// Check sigma is positive
        double s = sigmaInput.get().get();
        final double mean = 1;
        if (s <= 0) {
    		return Double.NEGATIVE_INFINITY;
    	}
        double m = Math.log(mean) - (0.5 * s * s);
        
        
        LogNormalDistribution dist = LogNormalDistribution.of(m, s);
        
        
        Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount()-1;
        
		
		double logp = 0;
        for (int nodeNr = 0; nodeNr < dimension; nodeNr ++) {
        	
        	
        	double rate = branchRates.get(nodeNr);
        	if (rate < 0) {
        		logp = Double.NEGATIVE_INFINITY;
        		return logp;
        	}

        	
        	logp += dist.logDensity(rate);
        	
        }
    	
        return logp;
        
    }

    
    public ScalarDistribution<?,?> getDist(){
    	
    	final double mean = 1;
    	double s = sigmaInput.get().get();;
    	
    	
    	RealScalar<Real> mParam = new RealScalarParam<>(mean, Real.INSTANCE);
    	RealScalar<PositiveReal> sParam = new RealScalarParam<>(s, PositiveReal.INSTANCE);
    	
    	
    	LogNormal dist = new LogNormal();
    	dist.initByName("M", mParam, "S", sParam, "meanInRealSpace", true);
    	
    	
    	return dist;
    	
    }


	@Override
	public List<Double> sample() {
		
		
		List<Double> rates = new ArrayList<Double>();
		
		Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount();
		
		
		 
        // Check sigma is positive
        double s = sigmaInput.get().get();
        final double mean = 1;
        if (s <= 0) {
        	throw new IllegalArgumentException("Cannot sample branch rates because sigma is non-positive " + s);
    	}
        double m = Math.log(mean) - (0.5 * s * s);
        LogNormalDistribution dist = LogNormalDistribution.of(m, s);
        
        
        for (int nodeNr = 0; nodeNr < dimension; nodeNr ++) {
        	
        	
        	try {
				double rateOfBranch = dist.inverseCumulativeProbability(Randomizer.nextFloat());
				rates.add(rateOfBranch);
				
			} catch (Exception e) {
				e.printStackTrace();
				throw new IllegalArgumentException("Unexpected error when sampling from LN(" + m + ", " + s + ")");
			}
        	
        }
        
        
        return rates;
		
		
	}
	


}
