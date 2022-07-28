package orc.operators;


import java.util.ArrayList;
import java.util.List;



import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.operator.Exchange;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.BeautiDoc;
import orc.consoperators.ConsOperatorUtils;
import orc.consoperators.InConstantDistanceOperator;
import orc.consoperators.PiecewiseLinearDistribution;


/**
 * A meta-operator which proposes a tree using Narrow Exchange and then recomputes branch rates using one of the NER variants in the orc.NER subpackage
 * Only applicable under a relaxed clock model
 * @author Jordan Douglas
 */
public class MetaNEROperator extends InConstantDistanceOperator {


	public final Input<ParametricDistribution> distributionInput = new Input<>("distr", "Rate distribution. Required if using the quantile parameterisation.");
	public final Input<Boolean> optimiseInput = new Input<>("optimize", "Whether to optimise the tunable parameter (default true)", true);




    // Time proposal
    protected double tdp;

    // Rate proposals
    protected double rap, rbp, rcp, rdp;

    // Quantile proposals
    protected double qap, qbp, qcp, qdp;


    // Clock model
    protected enum ClockMode {
        quantiles,
        rates
    }

    // Random walk size optimisation
    boolean optimise;


    // Inputs
    protected ClockMode clockMode;
    protected RealParameter rates;
    protected double twindowSize = 0;
    protected ParametricDistribution rateDistribution;
    protected Tree tree;
    protected KernelDistribution proposalKernel;



    public MetaNEROperator() {
    	twindowSizeInput.setRule(Validate.OPTIONAL);
    	clockModelInput.setRule(Validate.FORBIDDEN);
    }


	@Override
	public void initAndValidate() {


		// Real rates or quantile parameterisation?
		if (rateInput.get() != null) {
			rates = rateInput.get();
			clockMode = ClockMode.rates;
		}else {
			rates = quantileInput.get();
			clockMode = ClockMode.quantiles;
			rateDistribution = distributionInput.get();
			if (rateDistribution == null) {
				throw new IllegalArgumentException("Please specify the rate distribution 'distr' when using the quantile parameterisation");
			}
		}


		// Window size, tree, and sub-operator (use this operator if unspecified)
		this.twindowSize = twindowSizeInput.get() == null ? 0 : twindowSizeInput.get();
		this.tree =  treeInput.get();
		this.proposalKernel = proposalKernelInput.get();
		this.optimise = optimiseInput.get();


	}

	@Override
	public double proposal() {


		this.tree =  treeInput.get();

		// Get nodes which operator may apply to
        final List<Node> applicableNodesBeforeOperation = getApplicableNodes(this.tree);
        if (applicableNodesBeforeOperation.size() == 0) return Double.NEGATIVE_INFINITY;



        // Sample a node E uniformly at random
        int Eindex = Randomizer.nextInt(applicableNodesBeforeOperation.size());
        Node E = applicableNodesBeforeOperation.get(Eindex);


		// Access to the child nodes of E: C and D. D must not be a leaf.
        int nodeDNum = Randomizer.nextInt(2);
        Node D = E.getChild(nodeDNum);
        Node C = E.getChild(1 - nodeDNum);
        if (D.getHeight() < C.getHeight()) {
        	nodeDNum = 1 - nodeDNum;
            D = E.getChild(nodeDNum);
            C = E.getChild(1 - nodeDNum);
        }


        // Tree with dated tips
        if (D.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }


        // Randomly identify the two children of D as A and B
        Node A, B = null;
        if (Randomizer.nextBoolean()) {
        	A = D.getChild(0);
        	B = D.getChild(1);
        } else {
        	A = D.getChild(1);
        	B = D.getChild(0);
        }


        // Get original node times
        double ta = A.getHeight(); // Fixed
        double tb = B.getHeight(); // Fixed
        double tc = C.getHeight(); // Fixed
        double td = D.getHeight(); // Free
        double te = E.getHeight(); // Fixed

        // Make the proposal
        double logJD = 0;
        switch (clockMode) {

	        case rates: {

	            // Get node rates
	        	try {
		            double ra = rates.getValue(A.getNr()); // Free
		            double rb = rates.getValue(B.getNr()); // Free
		            double rc = rates.getValue(C.getNr()); // Free
		            double rd = rates.getValue(D.getNr()); // Free
		            
		            // Can be overriden by an NER variant
		            logJD = this.proposalRates(twindowSize, ta, tb, tc, td, te, ra, rb, rc, rd);
	        	} catch (ArrayIndexOutOfBoundsException e) {
	        		// can happen when third party operator does not put the root as highest numbered node
	        		return Double.NEGATIVE_INFINITY;
	        	}


	            // Ensure that the constraints have not been broken
	            if (!this.validateProposalRates(ta, tb, tc, td, te)) return Double.NEGATIVE_INFINITY;
	            break;

	        }

	        case quantiles: {

	        	// Get node quantiles
	            double qa = rates.getValue(A.getNr()); // Free
	            double qb = rates.getValue(B.getNr()); // Free
	            double qc = rates.getValue(C.getNr()); // Free
	            double qd = rates.getValue(D.getNr()); // Free


	            logJD = this.proposalQuantiles(twindowSize, ta, tb, tc, td, te, qa, qb, qc, qd);

	            // Ensure that the constraints have not been broken
	            if (!this.validateProposalQuantiles(ta, tb, tc, td, te)) return Double.NEGATIVE_INFINITY;
	            break;
	        }

        }

        // Irreversible proposal
        if (logJD == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;



        // Rearrange the tree
        exchangeNodes(A, C, D, E);





        // Set the new times + rates/quantiles
        D.setHeight(this.tdp);
        switch (clockMode) {
	        case rates: {

	            // Set the new rates
	        	rates.setValue(A.getNr(), this.rap);
	        	rates.setValue(B.getNr(), this.rbp);
	        	rates.setValue(C.getNr(), this.rcp);
	        	rates.setValue(D.getNr(), this.rdp);
	        	break;
	        }
	        case quantiles: {

	        	// Set the new quantiles
	        	rates.setValue(A.getNr(), this.qap);
	        	rates.setValue(B.getNr(), this.qbp);
	        	rates.setValue(C.getNr(), this.qcp);
	        	rates.setValue(D.getNr(), this.qdp);
	        	break;

	        }

        }




        // Hastings ratio
        double proposalForward = 0;
        double proposalBackward = 0;
        final List<Node> applicableNodesAfterOperation = getApplicableNodes(this.tree);
    	proposalForward = -Math.log(applicableNodesBeforeOperation.size());
    	proposalBackward = -Math.log(applicableNodesAfterOperation.size());
        double logHastingsRatioExchange = proposalBackward - proposalForward;





		//return Double.POSITIVE_INFINITY;
        return logHastingsRatioExchange + logJD;

	}



    // p1 becomes the parent of c2
    // p2 becomes the parent of c1
	private void exchangeNodes(Node c1, Node c2, Node p1, Node p2) {
        replace(p1, c1, c2);
        replace(p2, c2, c1);
    }


	/**
	 *
	 * @return A list of all nodes which this operator may apply to
	 */
    private List<Node> getApplicableNodes(Tree tree){


    	// Get all internal/root nodes which have grandchildren
    	List<Node> nodes = new ArrayList<Node>();


        // Iterate through all internal nodes
        for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {

        	Node E = tree.getNode(i);

        	// Ensure that at least 1 child is not a leaf
        	if (E.getChildCount() == 2 && (!E.getChild(0).isLeaf() || !E.getChild(1).isLeaf())) {
        		nodes.add(E);
        	}

        }

    	return nodes;

    }


    /**
     * @param windowSize - tunable step size
     * @return Sample the step size made during the random walk.
     * If there is no proposal kernel distribution then return 0
     */
    protected double getRandomWalkStepSize(double windowSize) {

    	// Sample a random walk step size
    	if (proposalKernel != null) return proposalKernel.getRandomDelta(1, windowSize);

    	// No random walk
    	return 0;


    }




	/**
	 * Makes proposals: td', ra', rb', rc', and rd' and returns the hastings ratio
	 * Assumes real rate parameterisation
	 * The proposal variables tdp, rap, rbp, rcp, and rdp contain the new values
	 * @param rWindowSize - proposal width size
	 * @param ta - original time of node A
	 * @param tb - original time of node B
	 * @param tc - original time of node C
	 * @param td - original time of node D
	 * @param te - original time of node E
	 * @param ra - original rate of node A
	 * @param rb - original rate of node B
	 * @param rc - original rate of node C
	 * @param rd - original rate of node D
	 * @return logJD - the natural logarithm of the determinant of the Jacobian
	 */
	protected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te,
												  double ra, double rb, double rc, double rd) {


		double r_td = this.getRandomWalkStepSize(rWindowSize);


		// Null proposal: keep all rates constant. Keep tD constant unless a random walk is applied.
		this.tdp = td;
		this.rap = ra;
		this.rbp = rb;
		this.rcp = rc;
		this.rdp = rd + r_td;
		return 0;

	}

	/**
	 * Makes proposals: td', qa', qb', qc', and qd' and returns the hastings ratio
	 * Assumes quantile parameterisation
	 * The proposal variables tdp, qap, qbp, qcp, and qdp contain the new values
	 * @param rWindowSize - proposal width size
	 * @param ta - original time of node A
	 * @param tb - original time of node B
	 * @param tc - original time of node C
	 * @param td - original time of node D
	 * @param te - original time of node E
	 * @param qa - original quantile of node A
	 * @param qb - original quantile of node B
	 * @param qc - original quantile of node C
	 * @param qd - original quantile of node D
	 * @return logJD - the natural logarithm of the determinant of the Jacobian
	 */
	protected double proposalQuantiles(double rWindowSize, double ta, double tb, double tc, double td, double te,
			  double qa, double qb, double qc, double qd) {


		double logJD = 0;
		try {


			// Convert quantiles to rates using the i-cdf
			double ra = rateDistribution.inverseCumulativeProbability(qa);
			double rb = rateDistribution.inverseCumulativeProbability(qb);
			double rc = rateDistribution.inverseCumulativeProbability(qc);
			double rd = rateDistribution.inverseCumulativeProbability(qd);


			// Propose new rates + times
			logJD = this.proposalRates(rWindowSize, ta,  tb,  tc, td, te, ra, rb, rc, rd);
			if (logJD == Double.NEGATIVE_INFINITY) return logJD;


			// Convert proposed rates into to proposed quantiles using the cdf
			this.qap = rateDistribution.cumulativeProbability(this.rap);
			this.qbp = rateDistribution.cumulativeProbability(this.rbp);
			this.qcp = rateDistribution.cumulativeProbability(this.rcp);
			this.qdp = rateDistribution.cumulativeProbability(this.rdp);


			// Contribution of the cdf to the hastings ratio
			if (this.qap != qa) logJD += getQuantileHastingsRatioContribution(this.rap, qa, this.qap);
			if (this.qbp != qb) logJD += getQuantileHastingsRatioContribution(this.rbp, qb, this.qbp);
			if (this.qcp != qc) logJD += getQuantileHastingsRatioContribution(this.rcp, qc, this.qcp);
			if (this.qdp != qd) logJD += getQuantileHastingsRatioContribution(this.rdp, qd, this.qdp);


		} catch (Exception e) {
			return Double.NEGATIVE_INFINITY;
		}

		return logJD;


	}



	/**
	 * Get the contribution of a q -> r -> r' -> q' transformation to the hastings ratio
	 * Only required in the quantiles parameterisation
	 * @param rNew - the branch rate after the proposal
	 * @param qOld - the branch rate quantile before the proposal
	 * @param qNew - the branch rate quantile after the proposal
	 * @return
	 */
	protected double getQuantileHastingsRatioContribution(double rNew, double qOld, double qNew) {


		double logHR = 0;
		if (rateDistribution instanceof LogNormalDistributionModel) {
			logHR = ConsOperatorUtils.getHRForLN(rNew, qOld, rateDistribution);
        }

        else if (rateDistribution instanceof PiecewiseLinearDistribution) {
        	logHR = ConsOperatorUtils.getHRForPieceWise(rNew, qOld, qNew, rateDistribution);
        }

        else {
        	logHR = ConsOperatorUtils.getHRUseNumericApproximation(rNew, qOld, rateDistribution);
        }

		return logHR;


	}



    @Override
    public void optimize(double logAlpha) {
    	if (!this.optimise) return;
    	double delta = calcDelta(logAlpha);
        delta += Math.log(twindowSize);
        twindowSize = Math.exp(delta);
    }


    @Override
    public double getCoercableParameterValue() {
        return twindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        twindowSize = value;
    }



    /**
     * Ensure that the new state is value
	 * @param ta - Height of node A
	 * @param tb - Height of node B
	 * @param tc - Height of node C
	 * @param td - Height of node D
	 * @param te - Height of node E
     * @return True if valid, false if invalid
     */
	protected boolean validateProposalRates(double ta, double tb, double tc, double td, double te) {

		if (this.tdp > te || this.tdp < tc || this.tdp < tb || this.tdp < ta) return false;
		if (this.rap <= 0) return false;
		if (this.rbp <= 0) return false;
		if (this.rcp <= 0) return false;
		if (this.rdp <= 0) return false;
		return true;

	}


	/**
	 * Ensure that the new state is value
	 * @param ta - Height of node A
	 * @param tb - Height of node B
	 * @param tc - Height of node C
	 * @param td - Height of node D
	 * @param te - Height of node E
	 * @return True if valid, false if invalid
	 */
	protected boolean validateProposalQuantiles(double ta, double tb, double tc, double td, double te) {

		if (this.tdp > te || this.tdp < tc || this.tdp < tb || this.tdp < ta) return false;
		if (this.qap <= 0 || this.qap >= 1) return false;
		if (this.qbp <= 0 || this.qbp >= 1) return false;
		if (this.qcp <= 0 || this.qcp >= 1) return false;
		if (this.qdp <= 0 || this.qdp >= 1) return false;


		// Ensure that proposed quantiles are not associated with rates which go outside the rate boundaries
		if (rateDistribution instanceof PiecewiseLinearDistribution) {
            PiecewiseLinearDistribution piecewise = (PiecewiseLinearDistribution) rateDistribution;

            try {

	            double rmin = piecewise.getRangeMin();
	            double rmax = piecewise.getRangeMax();
	            if (this.rap < rmin || this.rap > rmax) return false;
	            if (this.rbp < rmin || this.rbp > rmax) return false;
	            if (this.rcp < rmin || this.rcp > rmax) return false;
	            if (this.rdp < rmin || this.rdp > rmax) return false;

            } catch (Exception e) {
            	e.printStackTrace();
            	return false;
            }
        }


		return true;

	}


	/**
	 * Remove all narrow exchange operators from the beauti doc
	 * @param doc
	 * @return
	 */
    public static void customConnector(BeautiDoc doc) {

    	try {

	    	// Find and destroy the narrow exchanges
	    	for (String str : doc.pluginmap.keySet()) {

	    		BEASTInterface obj = doc.pluginmap.get(str);
	    		if (obj instanceof Exchange) {

	    			Exchange operator = (Exchange) obj;
	    			if (operator.isNarrowInput.get()) {
	    				operator.m_pWeight.set(0.0);
	    			}


	    		}

	    	}

    	} catch (Exception e) {

		}

    }




}
