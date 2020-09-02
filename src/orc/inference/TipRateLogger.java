package orc.inference;


import java.io.PrintStream;
import beast.evolution.tree.Tree;
import beast.core.parameter.RealParameter;
import beast.core.parameter.IntegerParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.core.Input;
import beast.core.Description;
import beast.core.Function;
import beast.core.Loggable;
import beast.core.BEASTObject;

@Description("Logs all rates associated with the leaves of a tree. Number of dimensions = number of taxa. If clock is specified then logger willlog actual rates instead of the parameter itself")
public class TipRateLogger extends BEASTObject implements Loggable, Function
{
    public final Input<BranchRateModel> clockModelInput = new Input<>("clock", "Clock model for computing rates");
    public final Input<IntegerParameter> ratesCategoriesInput = new Input<>("categories", "Vector of rate categories");
    public final Input<RealParameter> ratesInput = new Input<>("rates", "Vector of rates/quantiles");
    public final Input<Tree> treeInput = new Input<>("tree", "tree for which the rates apply", Input.Validate.REQUIRED);
    
    
    
    boolean usingCategories;
    private BranchRateModel clockModel;
    private Tree tree;
    private IntegerParameter catRates;
    private RealParameter realRates;
    private int ntaxa;
    

    
    public void initAndValidate() {
        this.catRates = (IntegerParameter)this.ratesCategoriesInput.get();
        this.realRates = (RealParameter)this.ratesInput.get();
        this.clockModel = (BranchRateModel)this.clockModelInput.get();
        if (this.catRates == null && this.realRates == null && this.clockModel == null) {
            throw new IllegalArgumentException("Must specify categories or real/quantile rates or the clock model");
        }
        this.usingCategories = (this.catRates != null);
        this.tree = this.treeInput.get();
        this.ntaxa = this.treeInput.get().getLeafNodeCount();
    }
    
    public double[] calcValues() {
        final double[] array = new double[this.ntaxa];
        for (int i = 0; i < this.ntaxa; ++i) {
            if (this.clockModel != null) {
                array[i] = this.clockModel.getRateForBranch(this.tree.getNode(i));
            }
            else if (this.usingCategories) {
                array[i] = this.catRates.getArrayValue(i);
            }
            else {
                array[i] = this.realRates.getArrayValue(i);
            }
        }
        return array;
    }
    
    public int getDimension() {
        return this.ntaxa;
    }
    
    public double getArrayValue() {
        return this.calcValues()[0];
    }
    
    public double getArrayValue(final int n) {
        if (n > this.ntaxa) {
            throw new IllegalArgumentException();
        }
        return this.calcValues()[n];
    }
    
    public void init(final PrintStream printStream) {
        String id = this.getID();
        if (id == null) {
            id = "";
        }
        for (int i = 0; i < this.ntaxa; ++i) {
            printStream.print(id + ".taxon." + this.tree.getNode(i).getID() + "\t");
        }
    }
    
    public void log(final long n, final PrintStream printStream) {
        for (int i = 0; i < this.ntaxa; ++i) {
            double d;
            if (this.clockModel != null) {
                d = this.clockModel.getRateForBranch(this.tree.getNode(i));
            }
            else if (this.usingCategories) {
                d = this.catRates.getArrayValue(i);
            }
            else {
                d = this.realRates.getArrayValue(i);
            }
            printStream.print(d + "\t");
        }
    }
    
    public void close(final PrintStream printStream) {
    }
}