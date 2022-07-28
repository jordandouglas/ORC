package orc.consoperators;

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


import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math.MathException;

@Description("Small pulley: Propose a new genetic distance")
@Citation(value =
        "Zhang, R., Drummond, A. (2020) Improving the performance of Bayesian phylogenetic inference\n" +
                "  under relaxed clock models. BMC Evol Biol 20, 54", DOI = "https://doi.org/10.1186/s12862-020-01609-4",
        year = 2020, firstAuthorSurname = "Zhang")
public class SmallPulley extends TreeOperator {
    final public  Input<Double> dwindowSizeInput =
            new Input<>("dwindowSize", "the size of the window in Big Pulley");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<RealParameter> quantileInput = new Input<>("quantiles", "the quantiles of each branch rate.", Input.Validate.XOR,rateInput);
    final public Input<UCRelaxedClockModel> clockModelInput = new Input<>("clockModel", "relaxed clock model used to deal with quantiles", Input.Validate.REQUIRED);
    final public Input<KernelDistribution> proposalKernelInput = new Input<>("kernel", "Proposal kernel for a random walk on the genetic distances.");
	
    
    // Proposal kernel
    private KernelDistribution kernel;

    private double dwindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private enum rateMode {
        quantiles,
        rates
    }
    private rateMode mode = rateMode.rates;


    @Override
    public void initAndValidate() {
        dwindowSize = dwindowSizeInput.get();
        if (rateInput.get() == null) {
            quantiles = quantileInput.get();
            mode = rateMode.quantiles;
        } else {
            rates = rateInput.get();
            mode = rateMode.rates;
        }
        kernel = proposalKernelInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        ParametricDistribution rateDistribution = clockModelInput.get().rateDistInput.get();
        int branchCount = tree.getNodeCount() - 1; //the number of branches of the tree

        // original rates
        double r_j; double r_k;

        // the original node time
        double t_x, t_j, t_k;

        // the original quantiles
        double q_j = 0.5; double q_k = 0.5;

        // the proposed quantiles
        double q_j_ = 0.5; double q_k_ = 0.5;

        double hastingsRatio = 0.0;

        // Step1: get the root of the tree
        Node node = tree.getRoot();
        // root time
        t_x = node.getHeight();//get the time of this node

        // Step2: access to the child nodes of the root
        // son
        Node son = node.getChild(0); // get the left child of this node, i.e. son
        t_j = son.getHeight(); // node time of son

        // daughter
        Node daughter = node.getChild(1); // get the right child of this node, i.e. daughter
        t_k = daughter.getHeight(); // node time of daughter

        int sonNr = son.getNr(); // node number of son
        if (sonNr == branchCount) {
            sonNr = son.getTree().getRoot().getNr();
        }

        int dauNr = daughter.getNr(); // node time of daughter
        if (dauNr == branchCount) {
            dauNr = daughter.getTree().getRoot().getNr();
        }

        // get the rates on branches above son and daughter
        switch (mode) {
            case rates: {
                r_j = rates.getValue(sonNr); // rate of branch above son
                r_k = rates.getValue(dauNr); // rate of branch above daughter
                break;
            }

            case quantiles: {
                q_j = quantiles.getValue(sonNr);
                q_k = quantiles.getValue(dauNr);
                try {
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

        // d: the distance to be proposed
        // i.e. distance on the branch above son
        double d = r_j * (t_x - t_j);
        // D: total distances of branches above son and daughter
        double D = r_k * (t_x - t_k) + d;

        // Step3: propose new genetic distance
        double b;
        if (kernel != null) b = kernel.getRandomDelta(1, dwindowSize);
        else b = Randomizer.uniform(-dwindowSize, dwindowSize);
        double d_ = d + b;
        if (d_ == 0.0 || d_ == D) {
        return Double.NEGATIVE_INFINITY;
        }
        // reflect the proposed distance
        double err; double n; double r;
        if (d_ > D) {
            err = d_ - D;
            n = Math.floor(err / D);
            r = err - n * D;
            if (n % 2 == 0) {
                d_ = D - r;
            } else {
                d_ = r;
            }
        } else if (d_ < 0.0) {
            err = - d_;
            n = Math.floor(err / D);
            r = err - n * D;
            if (n % 2 == 0) {
                d_ = r;
            } else {
                d_ = - r;
            }
        }

        //Step 3: propose new rates
        double r_j_ = d_ / (t_x - t_j);
        double r_k_ = (D - d_) / (t_x - t_k);

        //Step 4: set the proposed new rates
        switch (mode) {
            case rates: {
                // set rates directly
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
                        if (r_j_ <= rmin || r_j_ >= rmax) return Double.NEGATIVE_INFINITY;
                        if (r_k_ <= rmin || r_k_ >= rmax) return Double.NEGATIVE_INFINITY;
                    }

                    // new quantiles of proposed rates
                    q_j_ = rateDistribution.cumulativeProbability(r_j_);
                    q_k_ = rateDistribution.cumulativeProbability(r_k_);
                    
                    if (q_j_ <= 0 || q_j_ >= 1) return Double.NEGATIVE_INFINITY;
                    if (q_k_ <= 0 || q_k_ >= 1) return Double.NEGATIVE_INFINITY;

                    // set quantiles
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

        //Step5: calculate the Hastings ratio
        switch (mode) {
            case rates: {
                break;
            }

            case quantiles: {
                if (rateDistribution instanceof CachedDistribution && ((CachedDistribution)rateDistribution).distrInput.get() instanceof LogNormalDistributionModel) {
                    hastingsRatio = ConsOperatorUtils.getHRForLN(r_j_, q_j, ((CachedDistribution)rateDistribution).distrInput.get())
                            + ConsOperatorUtils.getHRForLN(r_k_, q_k, ((CachedDistribution)rateDistribution).distrInput.get());
                } else if (rateDistribution instanceof LogNormalDistributionModel) {
                    hastingsRatio = ConsOperatorUtils.getHRForLN(r_j_, q_j, rateDistribution)
                            + ConsOperatorUtils.getHRForLN(r_k_, q_k, rateDistribution);
                }

                else if (rateDistribution instanceof PiecewiseLinearDistribution) {
                    if (((PiecewiseLinearDistribution)rateDistribution).distrInput.get() instanceof LogNormalDistributionModel) {
                        hastingsRatio = ConsOperatorUtils.getHRForLN(r_j_, q_j, ((PiecewiseLinearDistribution)rateDistribution).distrInput.get())
                                + ConsOperatorUtils.getHRForLN(r_k_, q_k, ((PiecewiseLinearDistribution)rateDistribution).distrInput.get());
                    } else {
                    	hastingsRatio = ConsOperatorUtils.getHRForPieceWise(r_j_, q_j, q_j_, rateDistribution)
                            + ConsOperatorUtils.getHRForPieceWise(r_k_, q_k, q_k_, rateDistribution);
                    }
                }

                else {
                    hastingsRatio = ConsOperatorUtils.getHRUseNumericApproximation(r_j_, q_j, rateDistribution)
                            + ConsOperatorUtils.getHRUseNumericApproximation(r_k_, q_k, rateDistribution);
                }

                break;
            }

            default: {

            }

        }
        return hastingsRatio;
    }


    @Override
    public double getCoercableParameterValue() {
        return dwindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        dwindowSize = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    
    
    @Override
    public double getTargetAcceptanceProbability() {
    	if (this.kernel == null) return super.getTargetAcceptanceProbability();
    	else if (kernel instanceof KernelDistribution.Bactrian || kernel instanceof KernelDistribution.Mirror) {
    		return 0.3;
    	}
    	return super.getTargetAcceptanceProbability();
    }

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(dwindowSize);
        dwindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = dwindowSize * ratio;

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



