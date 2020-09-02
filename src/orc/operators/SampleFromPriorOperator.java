package orc.operators;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;

import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;




/**
 * Samples a parameter from its prior distribution
 * If no prior is specified, then assumes a uniform distribution between parameter min and max
 * If there is no min/max then throws error
 * If there is more than one element in the parameter, then the number of elements sampled 
 * 		per iteration is sample from a Binomial(n = ndimensions, p) distribution where np is tunable
 * If a tree is also provided, then the parameter (eg. rates, population sizes) is assumed to correspond to tree
 * 		branches. In this case, parameters within the same neighbourhood on the tree are sampled together
 * @author Jordan Douglas
 *
 */
public class SampleFromPriorOperator extends Operator {
	
	
    final public Input<Parameter> paramInput = new Input<>("parameter", "the parameter sample", Input.Validate.REQUIRED);
    final public Input<ParametricDistribution> priorInput = new Input<>("prior", "the prior distribution of the parameter");
    final public Input<Tree> treeInput = new Input<>("tree", "the tree that the parameter belong to (if applicable)", Input.Validate.OPTIONAL);
    final public Input<Double> npInput = new Input<>("np", "tunable parameter describing the mean number of elements in the parameter vector to sample", 1.0);
    
    
    boolean realParam;
    Parameter parameter;
    ParametricDistribution prior;
    double np;

    Tree tree;
    boolean useTree;
    boolean useRootBranch;
    
    

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
		parameter = paramInput.get();
		prior = priorInput.get();
		tree = treeInput.get();
		np = npInput.get();
		this.validateNP(true);
		
		
		if (! (parameter instanceof RealParameter) && !(parameter instanceof IntegerParameter)) {
			throw new IllegalArgumentException("Parameters must be Real or Integer parameters!");
		}
		
		this.realParam = this.parameter instanceof RealParameter;
	
		// Check that parameter dimensions match tree dimensions
		if (tree != null) {
			
			
			// Each element tree node is associated with an element in the parameter vector
			if (parameter.getDimension() == tree.getNodeCount()) {
				useTree = true;
				useRootBranch = true;
			}
			
			
			// Each element tree node is associated with an element in the parameter vector, EXCEPT the root node
			else if (parameter.getDimension() == tree.getNodeCount() - 1) {
				useTree = true;
				useRootBranch = false;
			}
			
			
			// Illegal input
			else {
				throw new IllegalArgumentException("Please ensure that the dimension of the parameter " + parameter.getID() + " matches the node count in the tree. "
						+ "Otherwise, do not specify the tree.");
			}
			
		}
		
		// No tree
		else {
			useTree = false;
			useRootBranch = false;
		}
		
		
		this.validatePrior();
		
		
		
		
		
	}
	
	
	private void validatePrior() {
		
		// Check that there is a proper prior
		if (prior == null) {
			if (this.realParam) {
				RealParameter param = (RealParameter) this.parameter;
				if  (param.getLower() == Double.NEGATIVE_INFINITY || param.getUpper() == Double.POSITIVE_INFINITY) {
					throw new IllegalArgumentException("Please either provide a prior or a lower and upper limit for parameter " + parameter.getID());
				}
			}else {
				IntegerParameter param = (IntegerParameter) this.parameter;
				if  (param.getLower() == Integer.MIN_VALUE || param.getUpper() == Integer.MAX_VALUE) {
					throw new IllegalArgumentException("Please either provide a prior or a lower and upper limit for parameter " + parameter.getID());
				}
				
			}
			

		}
		
		
	}

	@Override
	public double proposal() {
		
		//System.out.println("---------- prior=" + prior.density(1.0) + " ----------");
		
		double logHR = 0;
		List<Integer> paramsToSample = this.sampleParamsToSample();
		
		this.validatePrior();
		
		
		for (Integer p : paramsToSample) {
			
			//System.out.println("Sampling parameter " + p);


			// After proposal
			try {
				
				
				double oldLogP = 0;
				double newLogP = 0;
				
				// Real / double
				if (this.realParam) {
					
					RealParameter param = (RealParameter) parameter;
					
					// Before proposal
					double oldX = param.getValue(p);
					oldLogP = prior == null ? 0 : prior.logDensity(oldX);
					
					
					// Sample x from the prior
					double u = Randomizer.nextFloat();
					double newX = prior == null ? Randomizer.nextFloat() * (param.getUpper() - param.getLower()) + param.getLower() : prior.inverseCumulativeProbability(u);
					
					// Calculate log-density of new x
					newLogP = prior == null ? 0 : prior.logDensity(newX);
					
					// Set new value
					param.setValue(p, newX);
					
					
					//System.out.println("From " + oldX + " to " + newX);
					
					
				}
				
				// Integer / int
				else {
					
					
					IntegerParameter param = (IntegerParameter) parameter;
					
					// Before proposal
					int oldX = param.getValue(p);
					oldLogP = prior == null ? 0 : prior.logDensity(oldX);
					
					
					// Sample x from the prior
					int newX = (int) (prior == null ? Randomizer.nextInt(param.getUpper() - param.getLower() + 1) + param.getLower() : prior.inverseCumulativeProbability(Randomizer.nextFloat()));
					
					// Calculate log-density of new x
					newLogP = prior == null ? 0 : prior.logDensity(newX);
					
					// Set new value
					param.setValue(p, newX);
					
					
					//System.out.println("From " + oldX + " to " + newX);
					
				}
					

				
				// Hastings ratio
				logHR += oldLogP - newLogP;
				
				
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		

		return logHR;
	}
	
	
	
	/**
	 * @return The elements in 'parameter' to sample (without replacement)
	 */
	private List<Integer> sampleParamsToSample() {
		
		
		List<Integer> sampledIndices = new ArrayList<Integer>();
	    double prob = Math.min(1.0, this.np / this.parameter.getDimension());
		
	    
	    // Incorporate tree structure
	    if (useTree) {
	    	
	    	
	    	
	    	
	    	// Sample the baseline node
	    	int baseNodeNum = Randomizer.nextInt(this.parameter.getDimension());
	    	Node node = tree.getNode(baseNodeNum);
	    	sampledIndices.add(baseNodeNum);
	    	
	    	//System.out.print("" + baseNodeNum + " -> ");
	    	
	    	// Follow the lineage up to the root and decide whether to include each other child on the way
	    	do {
	    		
	    		
	    		

	    		// Sample this node (and its other children) ??
	    		if (!node.isRoot() || this.useRootBranch) {
	    			
	    			
	    			// Sample the parent?
	    			Node parent = node.getParent();
	    			boolean sampleParent = parent != null && Randomizer.nextFloat() < prob;
	    			if (parent.isRoot() && !this.useRootBranch) sampleParent = false;
	    			if (sampleParent) sampledIndices.add(parent.getNr());
	    			
	    			
	    			// Sample any siblings of node?
	    			for (int c = 0; c < parent.getChildCount(); c++) {
	    				Node sibling = parent.getChild(c);
	    				if (sibling.equals(node)) continue;

	    				boolean sampleSibling = Randomizer.nextFloat() < prob;
		    			if (sampleSibling) sampledIndices.add(sibling.getNr());
	    				
	    			}
	    			
	    			
	    			//if (parent != null) System.out.print(parent.getNr() + (sampleParent ? "y" : "n") + " -> ");
	    			
	    			
	    			node = parent;
	    			
	    			
	    		}
	    		

	    		
	    	} while (!node.isRoot());
	    	

	    	
	    }
	    
	    
	    // No tree
	    else {
			for (int i = 0; i < this.parameter.getDimension(); i ++) {
				boolean toSample = Randomizer.nextFloat() < prob;
				if (toSample) {
					sampledIndices.add(i);
				}
			}
	    }
		
	    
	    
	    // Check that at least 1 thing is being changed. If not, then select one uniformly at random
	    if (sampledIndices.size() == 0) {
	    	sampledIndices.add(Randomizer.nextInt(this.parameter.getDimension()));
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
    	if (this.np > this.parameter.getDimension()) {
    		if (verbose) Log.warning("Setting np to " + this.parameter.getDimension() + " so that it is no larger than n");
    		this.np = this.parameter.getDimension();
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
    
    
    
    
    @Override
    public List<StateNode> listStateNodes() {
    	
    	// The tree (if specified) is not operated on
    	List<StateNode> stateNodes = new ArrayList<StateNode>();
    	if (this.realParam) {
			RealParameter param = (RealParameter) parameter;
			stateNodes.add(param);
    	}else {
			IntegerParameter param = (IntegerParameter) parameter;
			stateNodes.add(param);
    	}
    	
    	return stateNodes;
    }
    
    
    /*
    @Override
    public double getTargetAcceptanceProbability() {
    	
    	// Since each parameter is resampled from its prior, the larger the acceptance rate, the better
        return 1.0;
    }
    */

    
    
    

}
