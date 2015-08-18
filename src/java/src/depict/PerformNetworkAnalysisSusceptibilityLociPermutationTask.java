package depict;

import depict.math.IntegerDoubleObject;
import java.util.concurrent.Callable;

/**
 *
 * @author DEPICTdevelopers
 */
public class PerformNetworkAnalysisSusceptibilityLociPermutationTask implements Callable<depict.math.IntegerDoubleObject> {

    private int l;
    private int vm;
    private int geneL;
    private int perm;
    private int nrUniqueLoci;
    private int[] nrGenesPerUniqueLoci;
    private int[][] geneStartIndexPermutedLoci;
    private int[] geneIDToMatrixID;
    private int nrGenes;
    private depict.matrix.SymmetricShortDistanceMatrix matrix;
    private double[] matrixValToPValue;
    private short[][] matrixUniform;
    
    public PerformNetworkAnalysisSusceptibilityLociPermutationTask(int perm, int l, int vm, int geneL, int nrUniqueLoci, int[] nrGenesPerUniqueLoci, int[][] geneStartIndexPermutedLoci, int[] geneIDToMatrixID, int nrGenes, depict.matrix.SymmetricShortDistanceMatrix matrix, double[] matrixValToPValue, short[][] matrixUniform) {
        this.perm = perm;
        this.l = l;
        this.vm = vm;
        this.geneL = geneL;
        this.nrUniqueLoci = nrUniqueLoci;
        this.nrGenesPerUniqueLoci = nrGenesPerUniqueLoci;
        this.geneStartIndexPermutedLoci = geneStartIndexPermutedLoci;
        this.geneIDToMatrixID = geneIDToMatrixID;
        this.nrGenes = nrGenes;
        this.matrix = matrix;
        this.matrixValToPValue = matrixValToPValue;
        this.matrixUniform = matrixUniform;
    }

    @Override
    public IntegerDoubleObject call() throws Exception {
        double minLogPValueSumPerm = 0;
        for (int m = 0; m < nrUniqueLoci; m++) {
            if (l != m) {
                double pValueMin = 1;
                for (int vm = 0; vm < nrGenesPerUniqueLoci[m]; vm++) {

                    int geneM = geneIDToMatrixID[geneStartIndexPermutedLoci[perm][m] + vm];
                    while (geneM == geneL) {
                        geneM = (int) (Math.random() * (double) nrGenes);
                    }

                    double pValue = 1;
                    if (matrixUniform == null) {
                        pValue = matrixValToPValue[matrix.get(geneL, geneM)];
                    } else {
                        pValue = matrixValToPValue[(int) matrixUniform[geneL][geneM]];
                    }
                    if (pValue < pValueMin) {
                        pValueMin = pValue;
                    }


                }

                pValueMin = 1d - Math.pow(1d - pValueMin, (double) nrGenesPerUniqueLoci[m]);
                if (pValueMin < 1E-16) {
                    pValueMin = 1E-16;
                }
                minLogPValueSumPerm += -Math.log(pValueMin);


            }
        }
        return new depict.math.IntegerDoubleObject(perm, minLogPValueSumPerm);
    }
}
