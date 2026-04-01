open module orc {
    requires beast.pkgmgmt;
    requires beast.base;
    requires org.apache.commons.statistics.distribution;


    exports orc.consoperators;
    exports orc.inference;
    exports orc.ner;
    exports orc.operators;
    exports orc.distribution;

    provides beast.base.core.BEASTInterface with
        
        orc.consoperators.InConstantDistanceOperator,
        orc.consoperators.SmallPulley,
        orc.consoperators.SimpleDistance,
        orc.inference.TipRateLogger,
        orc.distribution.BranchRatePrior,
        orc.ner.NEROperator_dAB,
        orc.ner.NEROperator_dAB_dAC,
        orc.ner.NEROperator_dAB_dAC_dAE,
        orc.ner.NEROperator_dAB_dAC_dAE_dBC,
        orc.ner.NEROperator_dAB_dAC_dAE_dBC_dCE,
        orc.ner.NEROperator_dAB_dAC_dAE_dCE,
        orc.ner.NEROperator_dAB_dAC_dBC,
        orc.ner.NEROperator_dAB_dAC_dBC_dBE,
        orc.ner.NEROperator_dAB_dAC_dBC_dCE,
        orc.ner.NEROperator_dAB_dAC_dBE,
        orc.ner.NEROperator_dAB_dAC_dCE,
        orc.ner.NEROperator_dAB_dAE,
        orc.ner.NEROperator_dAB_dAE_dBC,
        orc.ner.NEROperator_dAB_dAE_dBC_dCE,
        orc.ner.NEROperator_dAB_dAE_dCE,
        orc.ner.NEROperator_dAB_dBC,
        orc.ner.NEROperator_dAB_dBC_dBE,
        orc.ner.NEROperator_dAB_dBC_dCE,
        orc.ner.NEROperator_dAB_dBE,
        orc.ner.NEROperator_dAB_dBE_dCE,
        orc.ner.NEROperator_dAB_dCE,
        orc.ner.NEROperator_dAC,
        orc.ner.NEROperator_dAC_dAE,
        orc.ner.NEROperator_dAC_dAE_dBC,
        orc.ner.NEROperator_dAC_dAE_dBC_dCE,
        orc.ner.NEROperator_dAC_dAE_dBE,
        orc.ner.NEROperator_dAC_dAE_dBE_dCE,
        orc.ner.NEROperator_dAC_dAE_dCE,
        orc.ner.NEROperator_dAC_dBC,
        orc.ner.NEROperator_dAC_dBC_dBE,
        orc.ner.NEROperator_dAC_dBC_dCE,
        orc.ner.NEROperator_dAC_dBE,
        orc.ner.NEROperator_dAC_dBE_dCE,
        orc.ner.NEROperator_dAC_dCE,
        orc.ner.NEROperator_dAE,
        orc.ner.NEROperator_dAE_dBC,
        orc.ner.NEROperator_dAE_dBC_dBE,
        orc.ner.NEROperator_dAE_dBC_dCE,
        orc.ner.NEROperator_dAE_dBE,
        orc.ner.NEROperator_dAE_dBE_dCE,
        orc.ner.NEROperator_dAE_dCE,
        orc.ner.NEROperator_dBC,
        orc.ner.NEROperator_dBC_dBE,
        orc.ner.NEROperator_dBC_dCE,
        orc.ner.NEROperator_dBE,
        orc.ner.NEROperator_dBE_dCE,
        orc.ner.NEROperator_dCE,
        orc.operators.MetaNEROperator,
        orc.operators.SampleFromPriorOperator;
}


