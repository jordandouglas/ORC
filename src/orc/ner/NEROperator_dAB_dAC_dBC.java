package orc.ner;

import beast.base.util.Randomizer;
import orc.operators.MetaNEROperator;

public class NEROperator_dAB_dAC_dBC extends MetaNEROperator {

	@Override
	protected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te, 
										double ra, double rb, double rc, double rd) {

		// Random proposals
		double r_alpha = 0;
		double r_beta = 0;
		double r_gamma = 0;
		double r_delta = 0;
		double r_epsilon = this.getRandomWalkStepSize(rWindowSize);


		// Propose new rates + times
		this.rap = (te*(r_delta + rd) + ra*ta - ra*td - (r_delta + rd)*(r_epsilon + td))/(ta - te);
		this.rbp = -(rb*(tb - td))/(r_epsilon - tb + td);
		this.rcp = -(rc*tc - rc*te + rd*td - rd*te)/(r_epsilon - tc + td);
		this.rdp = r_delta + rd;
		this.tdp = r_epsilon + td;


		// Jacobian determinant
		double JD = ((ta - td)*(tb - td)*(tc - te))/((ta - te)*(r_epsilon - tb + td)*(r_epsilon - tc + td));
		if (JD <= 0) return Double.NEGATIVE_INFINITY;
		return Math.log(JD);

	}

}
