package depict;

import depict.math.IntegerDoubleObject;
import java.util.concurrent.Callable;

/**
 *
 * @author DEPICTdevelopers
 */
public class PerformPathwayAnalysisSusceptibilityLociPermutationTask implements Callable<depict.math.IntegerDoubleObject> {

    private int perm;
    private int pathway;
    private int nrUniqueLoci;
    private int[] nrGenesPerUniqueLoci;
    private int[][] geneStartIndexPermutedLoci;
    private int[] geneIDToMatrixID;
    private int nrGenes;
    depict.matrix.ExpressionDataset dataset;
    
    public PerformPathwayAnalysisSusceptibilityLociPermutationTask(int perm, int pathway, int nrUniqueLoci, int[] nrGenesPerUniqueLoci, int[][] geneStartIndexPermutedLoci, int[] geneIDToMatrixID, int nrGenes, depict.matrix.ExpressionDataset dataset) {
        this.perm = perm;
        this.pathway = pathway;
        this.nrUniqueLoci = nrUniqueLoci;
        this.nrGenesPerUniqueLoci = nrGenesPerUniqueLoci;
        this.geneStartIndexPermutedLoci = geneStartIndexPermutedLoci;
        this.geneIDToMatrixID = geneIDToMatrixID;
        this.nrGenes = nrGenes;
        this.dataset = dataset;
    }

    public IntegerDoubleObject call() throws Exception {
        double[] zScores = new double[nrUniqueLoci];
        for (int m = 0; m < nrUniqueLoci; m++) {
            for (int vm = 0; vm < nrGenesPerUniqueLoci[m]; vm++) {
                int geneM = geneIDToMatrixID[geneStartIndexPermutedLoci[perm][m] + vm];
                double zScore = dataset.rawData[geneM][pathway];
                zScores[m]+=zScore;
            }
            zScores[m]/=nrGenesPerUniqueLoci[m];
            if (Double.isNaN(zScores[m])) {
                System.out.println(perm + "\t" + pathway + "\t" + m + "\t" + zScores[m] + "\t" + nrGenesPerUniqueLoci[m]);
            }
        }
        
        double mean = mean(zScores);
        double variance = variance(zScores, mean);
        double zScore = mean / Math.sqrt(variance);
        return new depict.math.IntegerDoubleObject(perm, zScore);
    }

    private double mean(double[] v) {
        double sum = 0;
        for (int k = 0; k < v.length; k++) sum += v[k];
        return (sum / (double) v.length);
    }
    
    private double variance(double[] v, double mean) {
        double ans = 0.0;
        for (int i = 0; i < v.length; i++) ans += (v[i] - mean) * (v[i] - mean);
        return ans / (v.length - 1);
    }
    
}
