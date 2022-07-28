package orc.consoperators;


import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math.MathException;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.CompoundRealParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("For internal nodes: propose a new node time")
@Citation(value =
        "Zhang, R., Drummond, A. (2020) Improving the performance of Bayesian phylogenetic inference\n" +
                "  under relaxed clock models. BMC Evol Biol 20, 54", DOI = "https://doi.org/10.1186/s12862-020-01609-4",
        year = 2020, firstAuthorSurname = "Zhang")
public class InConstantDistanceOperator extends TreeOperator {
	
	
	// Ensure this is always included among the two child branches
	enum ChildInclusions {
		none, 		// Operate on all any internal node (default)
		leaf,		// Make sure that the internal node always has one or more leaf children
		internal	// Make sure that the internal node always has one or more internal node children
	}
	
	
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<RealParameter> quantileInput = new Input<>("quantiles", "the quantiles of each branch rate.", Input.Validate.XOR,rateInput);
    final public Input<UCRelaxedClockModel> clockModelInput = new Input<>("clockModel", "relaxed clock model used to deal with quantiles", Input.Validate.REQUIRED);
    final public Input<KernelDistribution> proposalKernelInput = new Input<>("kernel", "Proposal kernel for a random walk on the internal node height.");
    
    final public Input<ChildInclusions> includeInput = new Input<>("include", "Option flag to ensure that this operator always targets nodes which have a 'leaf' or 'internal' child. "
    		+ "Default: no constraints 'none'.", ChildInclusions.none, ChildInclusions.values());
	
    
    
    // Proposal kernel
    private KernelDistribution kernel;
    
    private ChildInclusions include;
    private double twindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private enum rateMode {
        quantiles,
        rates
    }
    private rateMode mode = rateMode.rates;

    @Override
    public void initAndValidate() {
    	this.twindowSize = twindowSizeInput.get();
        if (rateInput.get() == null) {
        	this.quantiles = quantileInput.get();
            this.mode = rateMode.quantiles;
        } else {
        	this.rates = rateInput.get();
            this.mode = rateMode.rates;
        }
        this.kernel = proposalKernelInput.get();
        this.include = includeInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        ParametricDistribution rateDistribution = clockModelInput.get().rateDistInput.get();
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree

        // the chosen node to work on
        Node node;

        // the original node times
        double t_x, t_j, t_k;

        // the original rates
        double r_x, r_k, r_j;

        // the original quantiles
        double q_x = 0.5; double q_k = 0.5; double q_j = 0.5;

        // the proposed quantiles
        double q_x_ = 0.5; double q_k_ = 0.5; double q_j_ = 0.5;

        double hastingsRatio = 0.0;

        // Step 1: randomly select an internal node, denoted by node x.
        // avoid fake nodes used to connect direct ancestors into tree.
        node = this.sampleNode(tree);
        if (node == null) return Double.NEGATIVE_INFINITY;

       // the number of this node
        int nodeNr = node.getNr();
        // if this node has max number, then use the free index stored in root node to get rate.
        if (nodeNr == branchCount) {
            nodeNr = node.getTree().getRoot().getNr();
        }

       // time for this node
       t_x = node.getHeight();


       // Step 2: Access to the child nodes of this node
       // son
       Node son = node.getChild(0); // get the left child of this node, i.e. son
       t_j = son.getHeight(); // node time of son
       int sonNr = son.getNr(); // node number of son
       if (sonNr == branchCount) {
           sonNr = son.getTree().getRoot().getNr();
        }

       // daughter
       Node daughter = node.getChild(1); //get the right child of this node, i.e. daughter
       t_k = daughter.getHeight(); // node time of daughter
       int dauNr = daughter.getNr(); // node number of daughter
       if (dauNr == branchCount) {
            dauNr = daughter.getTree().getRoot().getNr();
       }

       switch (mode) {
           case rates: {
               r_x = rates.getValue(nodeNr); // rate of branch above this node
               r_j = rates.getValue(sonNr); // rate of branch above son
               r_k = rates.getValue(dauNr); // rate of branch above daughter
               break;
           }

           case quantiles: {
               q_x = quantiles.getValue(nodeNr);
               q_j = quantiles.getValue(sonNr);
               q_k = quantiles.getValue(dauNr);
               try {
                   r_x = rateDistribution.inverseCumulativeProbability(q_x);
                   r_j = rateDistribution.inverseCumulativeProbability(q_j);
                   r_k = rateDistribution.inverseCumulativeProbability(q_k);
               } catch (MathException e) {
                   e.printStackTrace();
                   return Double.NEGATIVE_INFINITY;
               }
               break;
           }

           default: {
               return Double.NEGATIVE_INFINITY;
           }
       }


       // Step3: to propose a new node time for this node
       double a;
       if (kernel != null) a = kernel.getRandomDelta(1, twindowSize);
       else a = Randomizer.uniform(-twindowSize, twindowSize);
       double t_x_ = t_x + a;

       // deal with the boundary cases
       double upper = node.getParent().getHeight();
       double lower = Math.max(t_j, t_k);

        if (t_x_ == lower || t_x_ == upper) {
            return Double.NEGATIVE_INFINITY;
        }
        // fold the proposed node time
        double err; double n; double r;
        if (t_x_ > upper) {
            err = t_x_ - upper;
            n = Math.floor(err / (upper - lower));
            r = err - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = upper - r;
            } else {
                t_x_ = lower + r;
            }
        } else if (t_x_ < lower) {
            err = lower - t_x_;
            n = Math.floor(err / (upper - lower));
            r = err - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = lower + r;
            } else {
                t_x_ = upper - r;
            }
        }

        // set the proposed node time
        node.setHeight(t_x_);


       // Step4: propose the new rates
       // there are three rates in total
       // r_x_, r_j_, r_k_
       double r_x_ = r_x * (upper - t_x) / (upper - t_x_);
       double r_j_ = r_j * (t_x - t_j) / (t_x_ - t_j);
       double r_k_ = r_k * (t_x - t_k) / (t_x_ - t_k);


       // set the proposed new rates or quantiles
        switch (mode) {
            case rates: {
                // set rates directly
                rates.setValue(nodeNr, r_x_);
                rates.setValue(sonNr, r_j_);
                rates.setValue(dauNr, r_k_);
                break;
            }

            case quantiles: {
                try {
                    // reject rates if exceeding piecewise approximation's range
                    if (rateDistribution instanceof PiecewiseLinearDistribution) {
                        PiecewiseLinearDistribution piecewise = (PiecewiseLinearDistribution) rateDistribution;
                        double rmin = piecewise.getRangeMin();
                        double rmax = piecewise.getRangeMax();
                        if (r_x_ <= rmin || r_x_ >= rmax) return Double.NEGATIVE_INFINITY;
                        if (r_j_ <= rmin || r_j_ >= rmax) return Double.NEGATIVE_INFINITY;
                        if (r_k_ <= rmin || r_k_ >= rmax) return Double.NEGATIVE_INFINITY;
                    }

                    // new quantiles of proposed rates
                    q_x_ = rateDistribution.cumulativeProbability(r_x_);
                    q_j_ = rateDistribution.cumulativeProbability(r_j_);
                    q_k_ = rateDistribution.cumulativeProbability(r_k_);
                    
                    
                    if (q_x_ <= 0 || q_x_ >= 1) return Double.NEGATIVE_INFINITY;
                    if (q_j_ <= 0 || q_j_ >= 1) return Double.NEGATIVE_INFINITY;
                    if (q_k_ <= 0 || q_k_ >= 1) return Double.NEGATIVE_INFINITY;

                    // set quantiles
                    quantiles.setValue(nodeNr, q_x_);
                    quantiles.setValue(sonNr, q_j_);
                    quantiles.setValue(dauNr, q_k_);


                } catch (MathException e) {
                    e.printStackTrace();
                    return Double.NEGATIVE_INFINITY;
                }
                break;
            }

            default: {

            }

        }

        

        // Step4: calculate the Hastings ratio
        double nu =(upper - t_x) * (t_x - t_j) * (t_x - t_k) ;
        double de = (upper - t_x_) * (t_x_ - t_j) * (t_x_ - t_k);
        hastingsRatio = Math.log(nu / de);
        switch (mode) {
            case rates: {

                break;
            }

            case quantiles: {
                if (rateDistribution instanceof CachedDistribution && ((CachedDistribution)rateDistribution).distrInput.get() instanceof LogNormalDistributionModel) {
                    hastingsRatio += ConsOperatorUtils.getHRForLN(r_x_, q_x, ((CachedDistribution)rateDistribution).distrInput.get())
                            + ConsOperatorUtils.getHRForLN(r_j_, q_j, ((CachedDistribution)rateDistribution).distrInput.get())
                            + ConsOperatorUtils.getHRForLN(r_k_, q_k, ((CachedDistribution)rateDistribution).distrInput.get());
                } else if (rateDistribution instanceof LogNormalDistributionModel) {
                    hastingsRatio += ConsOperatorUtils.getHRForLN(r_x_, q_x, rateDistribution)
                                  + ConsOperatorUtils.getHRForLN(r_j_, q_j, rateDistribution)
                                  + ConsOperatorUtils.getHRForLN(r_k_, q_k, rateDistribution);
                }

                else if (rateDistribution instanceof PiecewiseLinearDistribution) {
                    if (((PiecewiseLinearDistribution)rateDistribution).distrInput.get() instanceof LogNormalDistributionModel) {
                        hastingsRatio += ConsOperatorUtils.getHRForLN(r_x_, q_x, ((PiecewiseLinearDistribution)rateDistribution).distrInput.get())
                                + ConsOperatorUtils.getHRForLN(r_j_, q_j, ((PiecewiseLinearDistribution)rateDistribution).distrInput.get())
                                + ConsOperatorUtils.getHRForLN(r_k_, q_k, ((PiecewiseLinearDistribution)rateDistribution).distrInput.get());
                    } else {
                    	hastingsRatio += ConsOperatorUtils.getHRForPieceWise(r_x_, q_x, q_x_, rateDistribution)
                                  + ConsOperatorUtils.getHRForPieceWise(r_j_, q_j, q_j_, rateDistribution)
                                  + ConsOperatorUtils.getHRForPieceWise(r_k_, q_k, q_k_, rateDistribution);
                    }
                }

                else {
                    hastingsRatio += ConsOperatorUtils.getHRUseNumericApproximation(r_x_, q_x, rateDistribution)
                                  + ConsOperatorUtils.getHRUseNumericApproximation(r_j_, q_j, rateDistribution)
                                  + ConsOperatorUtils.getHRUseNumericApproximation(r_k_, q_k, rateDistribution);
                }

                break;
            }

            default: {

            }

        }
        return hastingsRatio;
    }



    protected Node sampleNode(Tree tree) {
    	Node node = null;
    	int nodeCount = tree.getNodeCount();
    	
    	
    	if (tree.getLeafNodeCount() < 3) return null;
    	
    	switch (this.include) {
    	
    	
    		// Sample an internal node uniformly at random
	    	case none: {
	    		
	    		do {
		             final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
		             node = tree.getNode(nodeNr);
		        } while (node.isRoot() || node.isLeaf() || node.isFake());
	    		
	    		break;
	    	}
	    	
	    	// Sample an internal node such that at least one of its children are leaves
	    	case leaf:{
	
	    		boolean hasLeafChild;
	    		do {
		             final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
		             node = tree.getNode(nodeNr);
		             hasLeafChild = node.getChild(0).isLeaf() || node.getChild(1).isLeaf();
		        } while (node.isRoot() || node.isLeaf() || node.isFake() || !hasLeafChild);
	    		
	    		break;
	    	}
	    	
	    	
	    	// Sample an internal node such that at least one of its children are internal nodes
	    	// As long as there are at least 5 leaves (in the binary tree), such a node will always exist
	    	case internal:{
	    		
	    		if (tree.getLeafNodeCount() < 5) return null;
	    		
	    		boolean hasInternalChild;
	    		do {
		             final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
		             node = tree.getNode(nodeNr);
		             hasInternalChild = !node.getChild(0).isLeaf() || !node.getChild(1).isLeaf();
	    		} while (node.isRoot() || node.isLeaf() || node.isFake() || !hasInternalChild);
	    		
	    		break;
	    	}
    	
    	
    	}
    	
    	
    	 return node;
    	 
	}


	// Tuning the parameter: twindowsize represents the range of Uniform distribution
    @Override
    public double getCoercableParameterValue() {
        return twindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        twindowSize = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(twindowSize);
        twindowSize = Math.exp(delta);
    }
    
    
    @Override
    public double getTargetAcceptanceProbability() {
    	if (this.kernel == null) return super.getTargetAcceptanceProbability();
    	else if (kernel instanceof KernelDistribution.Bactrian || kernel instanceof KernelDistribution.Mirror) {
    		return 0.3;
    	}
    	return super.getTargetAcceptanceProbability();
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = twindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }


    
    @Override
    public List<StateNode> listStateNodes() {
    	List<StateNode> stateNodes = super.listStateNodes();
    	boolean hasCompoundRealParameter = true;
    	while (hasCompoundRealParameter) {
    		hasCompoundRealParameter = false;
        	for (int i = 0; i < stateNodes.size(); i++) {
        		StateNode s = stateNodes.get(i);
        		if (s instanceof CompoundRealParameter) {
        			CompoundRealParameter c = (CompoundRealParameter) s;
        			stateNodes.remove(i);
        			stateNodes.addAll(c.parameterListInput.get());
        			i--;
        			hasCompoundRealParameter = true;
        		}
        	}
    	}
    	return stateNodes;
    }
    
}

