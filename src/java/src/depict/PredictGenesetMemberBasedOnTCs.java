package depict;
import java.util.*;
import java.io.*;


/**
 *
 * @author DEPICTdevelopers
 */
public class PredictGenesetMemberBasedOnTCs {

    public PredictGenesetMemberBasedOnTCs(String[] args) throws InterruptedException {

        cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.DRand();
        try {            

            String genenetworkFilename = "/Users/tp/Documents/Work/gene_network_components/affymetrix/GPL570-GPL96-GPL1261-GPL1355CombinedEigenvectorsENSGOnlyCollapsedAndNormalized"; 
            //"/Users/tp/Documents/Work/gene_network_components/RNA-seq/expression_table.genelevel.v71.htseq.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.PCAOverSamplesPrincipalComponents.txt"
            String genesetPath =  "/Users/tp/Drive/Work/Schizophrenia/Results/Reconstitution/v2/"; // There needs to be at least two gene sets, otherwise NaNs in cofunc.
            String cofuncPath = genesetPath;   
            String flag = "fmrp_presynaps";
            int genesetStart = 0; 
            int genesetEnd = 100;

/*
            String genenetworkFilename = args[0];  
            String genesetPath = args[1]; 
            String cofuncPath = args[2]; 
            String flag = args[3]; 
            int genesetStart = Integer.parseInt(args[4]); 
            int genesetEnd = Integer.parseInt(args[5]);
*/            
            String genesetFilename = genesetPath + "/" + flag + ".tab";
            String cofuncFilenamePrefix = cofuncPath + "/" + flag;
            int minGenesetSize = 10;
            int maxGenesetSize = 1000;
            
            depict.matrix.ExpressionDataset dataset = new depict.matrix.ExpressionDataset(genenetworkFilename);
            depict.matrix.ExpressionDataset datasetTransposed = new depict.matrix.ExpressionDataset(dataset.fileName);
            datasetTransposed.transposeDataset();
            
            //Identify strata in the dataset (mix of different species, that cause missing values per species):
            int nrStrata = 0;
            int[] nrTCsPerStrata = new int[dataset.nrSamples];
            nrTCsPerStrata[nrStrata]++;
            for (int s=1; s<dataset.nrSamples; s++) {
                for (int p=0; p<dataset.nrProbes; p++) {
                    if (Double.isNaN(dataset.rawData[p][s])!=Double.isNaN(dataset.rawData[p][s - 1])) {
                        nrStrata++;
                        break;
                    }
                }
                nrTCsPerStrata[nrStrata]++;
            }            
            nrStrata++;
            int[] strataStartColumn = new int[nrStrata];
            if (nrStrata>1) {
                for (int s=0; s<nrStrata; s++) {
                    if (s>0) strataStartColumn[s] = strataStartColumn[s-1] + nrTCsPerStrata[s-1];
                    System.out.println("Number of species:\t" + s + "\t" + nrTCsPerStrata[s] + "\t" + strataStartColumn[s]);
                }
            }

            // Read in gene sets
            // Determine size of gene set
            List<String> geneSets = new ArrayList<String>();
            try{
                Scanner scan = new Scanner(new File(genesetFilename));
                String line;
                while (scan.hasNextLine()) {
                    line = scan.nextLine();
                    
                    if (line.charAt(0) == "#".charAt(0)) { continue; } // Discard header

                    //Check whether there are at least 10 genes for each TC with data for this geneset:
                    int minNrGenesPerStrata = dataset.nrProbes;
                    for (int s=0; s<nrStrata; s++) {
                        int nrGenesAvailable = 0;
                        for (int gene=0; gene<dataset.nrProbes; gene++) {
                            if (!Double.isNaN(datasetTransposed.rawData[strataStartColumn[s]][gene])) {
                                //if (datasetGeneset.rawData[gene][geneset]==1) {
                                if (line.contains(dataset.probeNames[gene])) {
                                    nrGenesAvailable++;
                                }
                            }
                        }
                        minNrGenesPerStrata = Math.min(nrGenesAvailable, minNrGenesPerStrata);
                    }
                    if (minNrGenesPerStrata >= minGenesetSize && minNrGenesPerStrata <= maxGenesetSize) {
                        geneSets.add(line);
                    }                
                }
            } catch (Exception e) {
                e.printStackTrace();
                System.out.println(e.getMessage());
            }
            
            // Construct matrix with gene sets (binary representation)
            depict.matrix.ExpressionDataset datasetGeneset = new depict.matrix.ExpressionDataset(dataset.nrProbes, geneSets.size());
            datasetGeneset.probeNames = dataset.probeNames;
            datasetGeneset.recalculateHashMaps();            
            int genesetCounter = 0;
            for (String s : geneSets) {
                String[] split = s.split("\t"); 
                datasetGeneset.sampleNames[genesetCounter] = split[0] + "_" + split[1];
                String[] genesetGenes = Arrays.copyOfRange(split, 3, split.length);;

                datasetGeneset.recalculateHashMaps();
                HashMap hashGenesetGenes = new HashMap();
                for (int h=0; h<genesetGenes.length; h++) {
                    hashGenesetGenes.put(genesetGenes[h], null);
                }
                for (int p=0; p <datasetGeneset.nrProbes; p++) {
                    if (hashGenesetGenes.containsKey(datasetGeneset.probeNames[p])) {
                        datasetGeneset.rawData[p][genesetCounter] = 1;
                    } else {
                        datasetGeneset.rawData[p][genesetCounter] = 0;
                    }
                }
                genesetCounter++;
            }
            System.gc();
            
            // Construct confunctionality matrix
            depict.matrix.ExpressionDataset datasetCofunc = new depict.matrix.ExpressionDataset(dataset.nrProbes, geneSets.size());
            datasetCofunc.probeNames = dataset.probeNames;  // Ensemble gene IDs
            datasetCofunc.sampleNames = datasetGeneset.sampleNames; // Gene sets
            datasetCofunc.recalculateHashMaps();
            
            String[] queryGenesets = datasetCofunc.sampleNames;
            for (int q = genesetStart; q < Math.min(queryGenesets.length,genesetEnd); q++) {
                int geneset = ((Integer) datasetGeneset.hashSamples.get(queryGenesets[q])).intValue();

                // Count number of genes available
                int[] nrGenesInGenesetPerStrata = new int[nrStrata];    
                int[] nrGenesNotInGenesetPerStrata = new int[nrStrata];
                for (int s=0; s<nrStrata; s++) {
                    for (int gene=0; gene<dataset.nrProbes; gene++) { // Loop over genes
                        if (!Double.isNaN(datasetTransposed.rawData[strataStartColumn[s]][gene])) { // Make sure not NA in expression data
                            if (datasetGeneset.rawData[gene][geneset]==1) {
                                nrGenesInGenesetPerStrata[s]++;
                            } else {
                                nrGenesNotInGenesetPerStrata[s]++;
                            }
                        }
                    }
                }

                // Score each gene set across TCs
                double[] zScoreProfile = new double[dataset.nrSamples]; // Gene sets
                for (int s=0; s<nrStrata; s++) {

                    double[] vals1 = new double[nrGenesInGenesetPerStrata[s]];
                    double[] vals2 = new double[nrGenesNotInGenesetPerStrata[s]];
                    for (int tc=0; tc<nrTCsPerStrata[s]; tc++) {
                        int vals1Itr = 0;
                        int vals2Itr = 0;
                        for (int gene=0; gene<dataset.nrProbes; gene++) { // Genes
                            if (!Double.isNaN(datasetTransposed.rawData[tc + strataStartColumn[s]][gene])) { // Look up in expression
                                if (datasetGeneset.rawData[gene][geneset]==1) { // Gene is in gene set
                                    vals1[vals1Itr] = datasetTransposed.rawData[tc + strataStartColumn[s]][gene];
                                    vals1Itr++;
                                } else { // Gene is not in gene set 
                                    vals2[vals2Itr] = datasetTransposed.rawData[tc + strataStartColumn[s]][gene];
                                    vals2Itr++;
                                }
                            }
                        }
                        double n1 = vals1.length;
                        double n2 = vals2.length;
                        double mean1 = mean(vals1);
                        double mean2 = mean(vals2);
                        double var1 = variance(vals1, mean1);
                        double var2 = variance(vals2, mean2);
                        double t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);
                        double df = ((var1/n1+var2/n2) * (var1/n1+var2/n2)) / ( ((var1/n1) * (var1/n1)) / (n1-1) + ((var2/n2) * (var2/n2)) / (n2-1));
                        cern.jet.random.StudentT tDist = new cern.jet.random.StudentT(df, randomEngine);
                        double pValue = 1;
                        double zScore = 0;
                        if (t < 0) { // mean1 < mean2 results in negative z-score
                            pValue = tDist.cdf(t);
                            if (pValue < 2.0E-323) pValue = 2.0E-323;
                            zScore = cern.jet.stat.Probability.normalInverse(pValue);
                        } else { // mean1 > mean2 results in positive z-score
                            pValue = tDist.cdf(-t);
                            if (pValue < 2.0E-323) pValue = 2.0E-323;
                            zScore = -cern.jet.stat.Probability.normalInverse(pValue);
                        }
                        pValue*=2;
                        zScoreProfile[tc + strataStartColumn[s]] = zScore;
                    }

                }

                // Score each gene across TC's
                double[] zScoresIndividualGenes = new double[dataset.nrProbes];
                for (int p=0; p<dataset.nrProbes; p++) { // Loop through genes

                    double[] zScoreProfileToUse = null;

                    // Prevent overfitting, don't use the above z-score vector
                    if (datasetGeneset.rawData[p][geneset]==1) {

                        zScoreProfileToUse = new double[dataset.nrSamples]; // Gene sets
                        int[] nrGenesInGenesetPerStrataToUse = new int[nrStrata];
                        int[] nrGenesNotInGenesetPerStrataToUse = new int[nrStrata];
                        for (int s=0; s<nrStrata; s++) {
                            for (int gene=0; gene<dataset.nrProbes; gene++) { // loop through genes
                                if (p!=gene) {
                                    if (!Double.isNaN(datasetTransposed.rawData[strataStartColumn[s]][gene])) {
                                        if (datasetGeneset.rawData[gene][geneset]==1) {
                                            nrGenesInGenesetPerStrataToUse[s]++;
                                        } else {
                                            nrGenesNotInGenesetPerStrataToUse[s]++;
                                        }
                                    }
                                }
                            }
                        }                            
                        for (int s=0; s<nrStrata; s++) {
                            double[] vals1 = new double[nrGenesInGenesetPerStrataToUse[s]];
                            double[] vals2 = new double[nrGenesNotInGenesetPerStrataToUse[s]];
                            for (int tc=0; tc<nrTCsPerStrata[s]; tc++) {
                                int vals1Itr = 0;
                                int vals2Itr = 0;
                                for (int gene=0; gene<dataset.nrProbes; gene++) {
                                    if (p!=gene) {
                                        if (!Double.isNaN(datasetTransposed.rawData[tc + strataStartColumn[s]][gene])) {
                                            if (datasetGeneset.rawData[gene][geneset]==1) {
                                                vals1[vals1Itr] = datasetTransposed.rawData[tc + strataStartColumn[s]][gene];
                                                vals1Itr++;
                                            } else {
                                                vals2[vals2Itr] = datasetTransposed.rawData[tc + strataStartColumn[s]][gene];
                                                vals2Itr++;
                                            }
                                        }
                                    }
                                }
                                double n1 = vals1.length;
                                double n2 = vals2.length;
                                double mean1 = mean(vals1);
                                double mean2 = mean(vals2);
                                double var1 = variance(vals1, mean1);
                                double var2 = variance(vals2, mean2);
                                double t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);
                                double df = ((var1/n1+var2/n2) * (var1/n1+var2/n2)) / ( ((var1/n1) * (var1/n1)) / (n1-1) + ((var2/n2) * (var2/n2)) / (n2-1));
                                cern.jet.random.StudentT tDist = new cern.jet.random.StudentT(df, randomEngine);
                                double pValue = 1;
                                double zScore = 0;
                                if (t < 0) { // mean1 < mean2 results in negative z-score
                                    pValue = tDist.cdf(t);
                                    if (pValue < 2.0E-323) pValue = 2.0E-323;
                                    zScore = cern.jet.stat.Probability.normalInverse(pValue);
                                } else { // mean1 > mean2 results in positive z-score
                                    pValue = tDist.cdf(-t);
                                    if (pValue < 2.0E-323) pValue = 2.0E-323;
                                    zScore = -cern.jet.stat.Probability.normalInverse(pValue);
                                }
                                zScoreProfileToUse[tc + strataStartColumn[s]] = zScore;
                            }

                        }                            
                    } else {
                        zScoreProfileToUse = zScoreProfile;
                    }

                    // Compute correlation between gene and gene set
                    Vector vec1 = new Vector();
                    Vector vec2 = new Vector();
                    for (int tc=0; tc<dataset.nrSamples; tc++) {
                        if (!Double.isNaN(dataset.rawData[p][tc])) {
                            vec1.add(dataset.rawData[p][tc]);
                            vec2.add(zScoreProfileToUse[tc]);
                        }
                    }
                    double[] vals1 = new double[vec1.size() * 2];
                    double[] vals2 = new double[vec2.size() * 2];
                    for (int v=0; v<vec1.size(); v++) {
                        vals1[v * 2] = ((Double) vec1.get(v)).doubleValue();
                        vals2[v * 2] = ((Double) vec2.get(v)).doubleValue();
                        vals1[v * 2 + 1] = - ((Double) vec1.get(v)).doubleValue();
                        vals2[v * 2 + 1] = - ((Double) vec2.get(v)).doubleValue();
                    }
                    double correlation = JSci.maths.ArrayMath.correlation(vals1, vals2);

                    cern.jet.random.StudentT tDistColt = new cern.jet.random.StudentT(vals1.length/2 - 2, randomEngine);
                    double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (vals1.length/2 - 2)));
                    double pValue = 1;
                    double zScore = 0;
                    if (t < 0) { // Means negative correlation between gene and gene, results in negative z-score 
                        pValue = tDistColt.cdf(t);
                        if (pValue < 2.0E-323) pValue = 2.0E-323;
                        zScore = cern.jet.stat.Probability.normalInverse(pValue);
                    } else { // Positive correlation between gene and gene set, results in positive z-score (used in DEPICT output for top genes)
                        pValue = tDistColt.cdf(-t);
                        if (pValue < 2.0E-323) pValue = 2.0E-323;
                        zScore = -cern.jet.stat.Probability.normalInverse(pValue);
                    }
                    pValue*=2;

                    datasetCofunc.rawData[p][q] = zScore;
                    zScoresIndividualGenes[p] = zScore; // Used for AUC
                }

                // Compute AUC
                Vector vecInGeneset = new Vector();
                Vector vecNotInGeneset = new Vector();
                for (int p=0; p<dataset.nrProbes; p++) { // Loop through genes
                    if (datasetGeneset.rawData[p][geneset]==1) {
                        vecInGeneset.add(zScoresIndividualGenes[p]);
                    } else {
                        vecNotInGeneset.add(zScoresIndividualGenes[p]);
                    }
                }
                double[] vals1 = new double[vecInGeneset.size()];
                for (int v=0; v<vals1.length; v++) {
                    vals1[v] = ((Double) vecInGeneset.get(v)).doubleValue();
                }
                double[] vals2 = new double[vecNotInGeneset.size()];
                for (int v=0; v<vals2.length; v++) {
                    vals2[v] = ((Double) vecNotInGeneset.get(v)).doubleValue();
                }

                depict.math.WilcoxonMannWhitney wmw = new depict.math.WilcoxonMannWhitney();
                double pValue = wmw.returnWilcoxonMannWhitneyPValue(vals1, vals2);
                System.out.println("PredictionPerformanceGeneset:\t" + queryGenesets[q] + "\tGenesetNr:\t" + geneset + "\t\tNrGenesInGeneset\t" + vals1.length + "\tNrGenesOutsideGeneset:\t" + vals2.length + "\tP-Value:\t" + pValue + "\tAUC:\t" + wmw.getAUC());
            }            
            datasetCofunc.save(cofuncFilenamePrefix + ".txt");
            datasetCofunc.standardNormalizeData();
            datasetCofunc.transposeDataset();
            datasetCofunc.standardNormalizeData();
            datasetCofunc.transposeDataset();
            datasetCofunc.save(cofuncFilenamePrefix + "_normalized.binary");
            //datasetCofunc.save(cofuncFilenamePrefix + "_" + genesetStart + "-" + genesetEnd + ".txt");
            //datasetCofunc.save(cofuncFilenamePrefix + ".txt");
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }

        System.exit(0);

    }


    
    private double mean(double[] v) {
        double sum = 0;
        for (int k = 0; k < v.length; k++) {
            sum += v[k];
        }
        return (sum / (double) v.length);
    }

    private double variance(double[] v, double mean) {
        double ans = 0.0;
        for (int i = 0; i < v.length; i++) {
            ans += (v[i] - mean) * (v[i] - mean);
        }
        return ans / (v.length - 1);
    }
    
}
