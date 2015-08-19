package depict;

import depict.math.StringDoubleObject;
import depict.math.StringIntegerObjectSorterOnIntegerValues;
import depict.math.StringDoubleObjectSorter;
import depict.math.IntegerDoubleObject;
import depict.math.StringIntegerObject;
import depict.matrix.SymmetricShortDistanceMatrix;
import depict.matrix.ExpressionDataset;
import java.io.*;
import java.util.*;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author DEPICTdevelopers
 */
public class PerformPathwayAnalysisSusceptibilityLoci {

    public PerformPathwayAnalysisSusceptibilityLoci(String[] args) throws InterruptedException {

        Calendar cal = Calendar.getInstance();
        Date currentTime = cal.getTime();
        System.out.println("DEPICT run from " + currentTime + "\n");

        String directoryPermutedLociDefinedBySignificantSNPssAndLDInformation = args[0];
        String filenameLociDefinedBySignificantSNPssAndLDInformation = args[1];        
        String outputFileLabel = args[2];
        boolean conductNetworkAnalysis = Integer.parseInt(args[3]) == 1;
        boolean conductPathwayAnalysis = Integer.parseInt(args[4]) == 1;
        boolean conductTissueAnalysis = Integer.parseInt(args[5]) == 1;
        int nrCores = Integer.parseInt(args[6]);
        String resultsDirectory = args[7];
        String cofuncMatrixFile = args[8];
        String filenameGeneAnnotation = args[9];  
        String tissueMatrixFile = args[10];
        int maxTopGenesPerGeneSet = Integer.parseInt(args[11]);
        int nrReps = Integer.parseInt(args[12]);
        int nrPerms = Integer.parseInt(args[13]);
        int HLAstart = Integer.parseInt(args[14]);
        int HLAend = Integer.parseInt(args[15]);
        String goMappingFile = args[16];
        String mgiMappingFile = args[17];
        String inwebMappingFile = args[18];
        String tissueMappingFile = args[19];
        String filenameGenericIlluminaProbeIDToEnsembl = args[20]; //Additional annotation on SNPs that affect Ensembl genes. This is only used as additional annotation, not used in any of the calculations:
        String filenameGenericIlluminaProbeIDEQTLs = args[21];
        boolean calculateGenePrioritizationPValueForGenesOutsideLoci = Integer.parseInt(args[22]) == 1;
        String chrToBeLeftOut = args[23]; // Usefull for benchmarking purposes
        
        //Which GWAS gene P-Value file to use:
        String filenameGenePValues = null; 
        double genePValueThresholdToUse = 0.00;

        //Parameters to use:
        boolean enforceUniformPValueDistributionOnCofunctionalityEveryGene = false; // Do we want to enforce a uniform P-Value distribution for every gene: If we do this, for every gene the co-functionality distribution with all other genes is uniformly distributed
        boolean includeSexChromosomesWhenPermuting = false; //Do we want to include genes on chr. X and Y in the analysis? Most GWAS only assess autosomes, so usually we want to put this flag up
        boolean includeHLA = false; //Do we want to include genes within the HLA? This is conversatively defined as Chr 6, 20,000,000 - 40,000,000bp
        boolean correctForGeneSize = false; // This is taken care of by using background loci
        String correctInPermutationsForEffectiveNrTestedSNPsPerGene = null; //filenameGenePValues;
        String gzip_encoding = "US-ASCII";
        String gtex_identifier = "gtex";
        
        String filenameDatabaseToUse = cofuncMatrixFile;
        if (conductTissueAnalysis) {
            filenameDatabaseToUse = tissueMatrixFile;
        }

        // Maximum allowed top genes to be printed for a given gene set or tissue
        int limitMaxTopGenesPerGeneSet = 101;

        // Locus file columns
        int locusFileChrCol = 1;
        int locusFileNearestCol = 5;
        int locusFileGenesCol = 6;
        int locusFileGWASPvalueCol = 4;

        // Output files
        String outputGeneFileName = resultsDirectory + "/" + outputFileLabel + "_geneprioritization.txt";
        String outputGeneGenomeFileName = resultsDirectory + "/" + outputFileLabel + "_geneprioritization_outside_input_loci.txt";
        String outputGenesetFileName = resultsDirectory + "/" + outputFileLabel + "_genesetenrichment.txt";        
        String outputTissueFileName = resultsDirectory + "/" + outputFileLabel + "_tissueenrichment.txt";

        //Load the predicted Z-Scores for a certain database (e.g. GO_BP, MGI, Reactome, KEGG etc.):
        ExpressionDataset dataset = new ExpressionDataset(filenameDatabaseToUse, "\t", null, null);
        
        //Save cofunc as plain text or binary
        //dataset.save("/tmp/cofunc.binary");
        //dataset.save("/tmp/cofunc.txt");
        //System.exit(0);
        
        // Create mapping from gene set identifier to column in dataset (used to list top genes for a given gene set)
        Map<String, Integer> geneSetIDToMatrixID = new HashMap<String, Integer>();
        for (int gs = 0; gs < dataset.nrSamples; gs++) {
            geneSetIDToMatrixID.put(dataset.sampleNames[gs], gs);
        }

        // Print analyses that will be run
        if (conductNetworkAnalysis) {
            System.out.println("Performing gene prioritization\n");
        }
        if (conductTissueAnalysis) {
            System.out.println("Performing tissue/cell type enrichment analysis\n");
        } else {
            System.out.println("Performing reconstituted gene set enrichment analysis\n");
        }
        
        //Now load the gene annotation (HG 19, NCBI 37) from file. This file is sorted, such that the first gene on chr. 1 is first in the file, and the last gene on chr. Y is last in the file:
        HashMap hashEnsemblIndex = new HashMap(); //Hash: What is the index of every Ensembl gene?
        HashMap hashEnsemblToHGNC = new HashMap(); //Hash: Holding Ensembl to HGNC gene symbol conversion
        Vector vecEnsemblIndex = new Vector(); //Vector: Holding all Ensembl genes
        int[] ensemblChr = new int[70000]; //Chromosome 
        int[] ensemblChrPosStart = new int[70000]; //Chr position start of gene
        int[] ensemblChrPosEnd = new int[70000]; //Chr position end of gene
        String[] ensemblNames = new String[70000]; //Array holding all Ensembl names
        Vector vecEnsemblGeneLengths = new Vector(); //Vector: Contains objects that describe the gene and the length of every gene:
        HashMap hashEnsemblBiotype = new HashMap(); //HashMap: Contains information what type of gene we are dealing with (protein coding, lincRNA, antisense, etc.)
        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(filenameGeneAnnotation)));
            String str = in.readLine(); // Skip header line
            int itr = 0;
            int locus_col = 0;
            int start_pos_col = 2;
            int end_pos_col = 3;
            int biotype_col = 4;
            int hgnc_col = 5;
            int chr_col = 6;
            while ((str = in.readLine()) != null) {
                String[] data = str.split("\t");
                String probeset = data[locus_col];
                if (!hashEnsemblIndex.containsKey(probeset)) {
                    int chr; 
                    if (data[chr_col].equals("X")) {
                        chr = 23;
                    } else if (data[chr_col].equals("Y")) {
                        chr = 24;
                    } else { 
                        chr = Integer.parseInt(data[chr_col]);
                    }
                    int chrPosStart = Integer.parseInt(data[start_pos_col]);
                    int chrPosEnd = Integer.parseInt(data[end_pos_col]);
                    if (Arrays.asList(dataset.probeNames).contains(probeset) ) { // Only save if gene is in cofunctionality matrix
                        if (includeSexChromosomesWhenPermuting || chr <= 22) {
                            if (includeHLA || !(chr == 6 && chrPosEnd > HLAstart && chrPosStart < HLAend)) {
                                ensemblChr[itr] = chr;
                                ensemblChrPosStart[itr] = chrPosStart;
                                ensemblChrPosEnd[itr] = chrPosEnd;
                                ensemblNames[itr] = probeset; // Ensembl gene ID

                                //Determine the size of the gene. If we want to correct for the gene size in the analyses, we need this information:
                                int ensemblGeneLength = (int) Math.round((double) (chrPosEnd - chrPosStart) / 10d);
                                if (!correctForGeneSize) {
                                    ensemblGeneLength = 1;
                                }

                                StringIntegerObject object = new StringIntegerObject(probeset, ensemblGeneLength);
                                vecEnsemblGeneLengths.add(object);

                                hashEnsemblBiotype.put(probeset, data[biotype_col]);
                                String hgnc = data[hgnc_col];
                                hashEnsemblToHGNC.put(probeset, hgnc);
                                vecEnsemblIndex.add(probeset);
                                hashEnsemblIndex.put(probeset, itr);

                                itr++;
                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }
        int nrGenes = hashEnsemblIndex.size();
                
        System.out.println("Nr. of unique genes to be used in the analysis:\t" + vecEnsemblGeneLengths.size());

        //Check for NaN values, if found bail out:
        boolean nanValuesPresent = false;
        HashMap hashValidColumns = new HashMap();
        for (int s = 0; s < dataset.nrSamples; s++) {
            boolean columnNoNaNs = true;
            for (int p = 0; p < dataset.nrProbes; p++) {
                if (Double.isNaN(dataset.rawData[p][s])) {
                    System.out.println("NaN value found in network matrix!:\t" + dataset.probeNames[p] + "\t" + dataset.sampleNames[s] + "\t" + dataset.rawData[p][s]);
                    nanValuesPresent = true;
                    columnNoNaNs = false;
                }
            }
            if (columnNoNaNs) {
                hashValidColumns.put(dataset.sampleNames[s], null);
            }
        }
        if (nanValuesPresent) {
            //dataset = new ExpressionDataset(filenameDatabaseToUse, "\t", null, hashValidColumns);
            //dataset.save(dataset.fileName);
            System.out.println("ERROR!!! There are NaN values in the Z-Score matrix!!!");
            System.exit(0);
        }

        //Determine whether the correlation matrix has already been calculated for this particular database:
        String coexpressionSuffix = ".Coexpression.dat";
        File fileCoexpressionNetwork = new File(dataset.fileName + coexpressionSuffix);
        if (!fileCoexpressionNetwork.canRead()) {

            System.out.println("Correlation matrix has not been generated yet: Calculating this now. This will take a while!");
            dataset.standardNormalizeData();

            //Initialize symmetric matrix, in short format (2 bytes per gene-pair), to save space:
            depict.matrix.SymmetricShortDistanceMatrix matrix = new depict.matrix.SymmetricShortDistanceMatrix(dataset.nrProbes);

            final int nrSamples = dataset.nrSamples;
            double sampleCountMinusOne = nrSamples - 1;

            cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.DRand();
            cern.jet.random.StudentT tDistColt = new cern.jet.random.StudentT(dataset.nrSamples - 2, randomEngine);

            //Precalculate correlation to Z-Score index, to speed up things:
            int[] correlationToZScore = new int[2000001];
            for (int corrInt = 0; corrInt < 2000001; corrInt++) {
                double correlation = (double) (corrInt - 1000000) / 1000000d;
                double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (dataset.nrSamples - 2)));
                double pValue = 0;
                double zScore = 0;
                if (t < 0) { //Avoid complementation issues:
                    pValue = tDistColt.cdf(t);
                    if (pValue < 2.0E-323) {
                        pValue = 2.0E-323;
                    }
                    zScore = cern.jet.stat.Probability.normalInverse(pValue);
                } else {
                    pValue = tDistColt.cdf(-t); //Take two sided P-Value
                    if (pValue < 2.0E-323) {
                        pValue = 2.0E-323;
                    }
                    zScore = -cern.jet.stat.Probability.normalInverse(pValue);
                }
                int zScoreInt = 32768 + (int) Math.round((zScore * 100d)); //Symmetric matrix is defined as short, so start in the middle of range (32768)
                if (zScoreInt > matrix.MAX_VALUE) {
                    zScoreInt = matrix.MAX_VALUE - 1;
                }
                if (zScoreInt < 0) {
                    zScoreInt = 0;
                }
                correlationToZScore[corrInt] = zScoreInt;
            }

            //Calculate pairwise correlations and Z-Scores:
            final int nrProbes = dataset.nrProbes;
            for (int p = 0; p < nrProbes; p++) {
                double[] rawDataP = dataset.rawData[p];
                for (int q = p + 1; q < nrProbes; q++) {
                    double[] rawDataQ = dataset.rawData[q];
                    double covarianceInterim = 0;
                    for (int s = 0; s < nrSamples; s++) {
                        covarianceInterim += rawDataP[s] * rawDataQ[s];
                    }
                    double correlation = covarianceInterim / sampleCountMinusOne;
                    int corrInt = (int) Math.round(correlation * 1000000d + 1000000d);
                    //Store index:
                    matrix.set(p, q, correlationToZScore[corrInt]);
                }
                if (p % 100 == 0) {
                    System.out.println(p);
                }
            }

            System.out.println("Correlation matrix has been generated. Saving to file:\t" + dataset.fileName + coexpressionSuffix);
            matrix.save(new File(dataset.fileName + coexpressionSuffix));
        }

        //Now load the correlation matrix:
        depict.matrix.SymmetricShortDistanceMatrix matrix = new depict.matrix.SymmetricShortDistanceMatrix(dataset.nrProbes);
        matrix.load(new File(dataset.fileName + coexpressionSuffix));

        //Define array, allows quick look-up from matrix short index value to P-Value:
        double[] matrixValToPValue = new double[matrix.maxValue() + 1];
        for (int v = 0; v < matrix.maxValue(); v++) {
            double z = ((double) v - 32768d) / 100d;
            double pValue = cern.jet.stat.Probability.normal(-z);
            if (pValue < 2.0E-323) {
                pValue = 2.0E-323;
            }
            matrixValToPValue[v] = pValue;
        }

        //Do we want to enforce a uniform P-Value distribution?
        short[][] matrixUniform = null;
        if (enforceUniformPValueDistributionOnCofunctionalityEveryGene) {

            System.out.println("Enforcing uniform P-Value distribution on the co-functionality correlation matrix:");

            //Enforce normal distribution on every gene.
            matrixUniform = new short[dataset.nrProbes][dataset.nrProbes];
            for (int p = 0; p < dataset.nrProbes; p++) {


                //For every gene, store the Z-Score index as absolute value:
                double[] valsAbs = new double[dataset.nrProbes - 1];
                int itr = 0;
                for (int q = 0; q < dataset.nrProbes; q++) {
                    if (p != q) {
                        valsAbs[itr] = Math.abs(matrix.get(p, q) - 32768);
                        itr++;
                    }
                }
                //Now get ranks, instead of Z-Scores:
                jsc.util.Rank rank = new jsc.util.Rank(valsAbs, 0d);
                valsAbs = rank.getRanks();
                itr = 0;
                for (int q = 0; q < dataset.nrProbes; q++) {
                    if (p != q) {
                        matrixUniform[p][q] = (short) (Math.round(valsAbs[itr]) - 1);
                        itr++;
                    }
                }
                if (p % 1000 == 0) {
                    System.out.println(p);
                }
            }
            //Update the array for getting a P-Value:
            for (int v = 0; v < matrix.maxValue(); v++) {
                matrixValToPValue[v] = 0;
            }
            for (int p = 0; p < dataset.nrProbes; p++) {
                double pValue = 1 - (0.5d + (double) p) / (double) (dataset.nrProbes + 1);
                matrixValToPValue[p] = pValue;
            }
            System.out.println("");
        }

        //Now initialize the vector, indicating for every gene that will be included for analysis what the index in the correlation matrix is:
        //This is important, as we sometimes do not want to include genes on Chr. X or Chr. Y, or only want to do analysis on a subset of the data:
        int[] geneIDToMatrixID = new int[nrGenes];
        for (int g = 0; g < nrGenes; g++) {
            geneIDToMatrixID[g] = ((Integer) dataset.hashProbes.get((String) vecEnsemblIndex.get(g))).intValue();
        }


        //Load eQTL annotation information, not used at all for any calculations in the method, but only for output purposes!!!
        HashMap hashGenericIlluminaProbeIDToEnsembl = new HashMap();
        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(filenameGenericIlluminaProbeIDToEnsembl)));
            String str = in.readLine();
            int itr = 0;
            while ((str = in.readLine()) != null) {
                String[] data = str.split("\t");
                if (data[4].trim().length() > 1) {
                    if (dataset.hashProbes.containsKey(data[4].trim())) {
                        hashGenericIlluminaProbeIDToEnsembl.put(data[0], data[4].trim());
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }
        HashMap hashEnsemblGeneCisEQTLSNP = new HashMap();
        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(filenameGenericIlluminaProbeIDEQTLs)));
            String str = in.readLine();
            int itr = 0;
            while ((str = in.readLine()) != null) {
                String[] data = str.split("\t");
                if (hashGenericIlluminaProbeIDToEnsembl.containsKey(data[4])) {
                    String ensemblGene = (String) hashGenericIlluminaProbeIDToEnsembl.get(data[4]);
                    if (!hashEnsemblGeneCisEQTLSNP.containsKey(ensemblGene)) {
                        hashEnsemblGeneCisEQTLSNP.put(ensemblGene, data[1]);
                    } else {
                        hashEnsemblGeneCisEQTLSNP.put(ensemblGene, (String) hashEnsemblGeneCisEQTLSNP.get(ensemblGene) + ";" + data[1]);
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }



        //Gene definitions have been loaded, now load the locus information!!!
        if (correctInPermutationsForEffectiveNrTestedSNPsPerGene != null) {
            vecEnsemblGeneLengths = new Vector();
            try {
                java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(correctInPermutationsForEffectiveNrTestedSNPsPerGene)));
                String str = in.readLine();
                while ((str = in.readLine()) != null) {
                    String[] data = str.split("\t");
                    if (hashEnsemblIndex.containsKey(data[0])) {
                        try {
                            int nrEffectiveTestsTimes100 = (int) Math.round(Double.parseDouble(data[5]) * 100d);
                            StringIntegerObject object = new StringIntegerObject(data[0], nrEffectiveTestsTimes100);
                            //System.out.println(data[0] + "\t" + nrEffectiveTestsTimes100);
                            vecEnsemblGeneLengths.add(object);
                        } catch (Exception e) {
                        }
                    }
                }
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }
            System.out.println("In permutations there will be a correction for the effective number of tested SNPs per gene");
            System.out.println("Number of genes with effective number of tested SNPs available:\t" + vecEnsemblGeneLengths.size());
            if (vecEnsemblGeneLengths.size() == 0) {
                System.out.println("Error! Could not parse the number of effective tests per gene!!!!");
                System.exit(0);
            }
        }

        //Make an array that we can use to randomly pick genes, while adhering to differences in size they have or the effective number of tested SNPs per gene:
        StringIntegerObjectSorterOnIntegerValues sorter = new StringIntegerObjectSorterOnIntegerValues();
        sorter.sort(vecEnsemblGeneLengths);
        int totalLengthGenes = 0;
        for (int g = 0; g < vecEnsemblGeneLengths.size(); g++) {
            StringIntegerObject object = (StringIntegerObject) vecEnsemblGeneLengths.get(g);
            int geneLength = object.intValue;
            totalLengthGenes += geneLength;
        }
        //Array randomGene will hold references to gene identifiers. When randomly picking an element from this array, we will randomly pick a gene, but in such a way that we more often pick bigger genes than smaller genes.
        short[] randomGene = new short[totalLengthGenes];
        totalLengthGenes = 0;
        for (int g = 0; g < vecEnsemblGeneLengths.size(); g++) {
            StringIntegerObject object = (StringIntegerObject) vecEnsemblGeneLengths.get(g);
            int geneLength = object.intValue;
            int geneIndex = ((Integer) hashEnsemblIndex.get(object.stringValue)).intValue();
            for (int l = 0; l < geneLength; l++) {
                randomGene[totalLengthGenes + l] = (short) geneIndex;
            }
            totalLengthGenes += geneLength;
        }

        
        Vector[] vecLociPerGene = new Vector[nrGenes];
        for (int e = 0; e < nrGenes; e++) {
            vecLociPerGene[e] = new Vector();
        }
        Vector vecLoci = new Vector();
        Vector vecLociNames = new Vector();
        int nrNonUniqueGenesInLoci = 0;
        HashMap hashGWASLocusPValue = null;
        HashMap hashClosestGeneToSNP = null;
        if (filenameLociDefinedBySignificantSNPssAndLDInformation != null) {
            if (hashClosestGeneToSNP == null) {
                hashClosestGeneToSNP = new HashMap();
            }
            try {
                java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(filenameLociDefinedBySignificantSNPssAndLDInformation)));                
                String str = in.readLine(); // Skip header
                while ((str = in.readLine()) != null) {
                    String[] data = str.replace(" ", "").split("\t");
                    if ( str.startsWith("#") || data[locusFileChrCol].equals(chrToBeLeftOut) ) {
                        // Skip this if
                        //   - comment line in locus definition file
                        //   - This locus is on the chr that was chosen to be left out
                    } else {
                        hashClosestGeneToSNP.put(data[0], data[locusFileGWASPvalueCol]);
                        Vector vecGenes = new Vector();
                        for (int d = locusFileNearestCol; d < data.length; d++) {
                            if (d == locusFileNearestCol) {
                                String[] data2 = data[locusFileNearestCol].split(";");
                                for (int e = 0; e < data2.length; e++) {                               
                                    if (hashEnsemblIndex.containsKey(data2[e])) {
                                        hashClosestGeneToSNP.put(data2[e], null);
                                        if (data.length == locusFileNearestCol+1 && !vecGenes.contains(data2[e])) {
                                            int ensemblIndex = ((Integer) hashEnsemblIndex.get(data2[e])).intValue();
                                            vecLociPerGene[ensemblIndex].add(vecLoci.size());
                                            vecGenes.add(data2[e]);
                                        }
                                    }
                                }
                            }
                            if (d == locusFileGenesCol) {
                                String[] data2 = data[locusFileGenesCol].split(";");
                                for (int e = 0; e < data2.length; e++) {
                                    if (hashEnsemblIndex.containsKey(data2[e])) {
                                        if (!vecGenes.contains(data2[e])) {
                                            int ensemblIndex = ((Integer) hashEnsemblIndex.get(data2[e])).intValue();
                                            vecLociPerGene[ensemblIndex].add(vecLoci.size());
                                            vecGenes.add(data2[e]);
                                        }
                                    }
                                }
                            }
                        }
                        if (vecGenes.size() > 0) {
                            vecLociNames.add(data[0]);
                            vecLoci.add(vecGenes);
                        }
                    }
                }
                //System.out.println("Initial number of (potentially overlapping) loci:\t" + vecLoci.size());
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
                System.exit(-1);
            }
        }

        //Initial genes for prioritization have been selected, now ascertain whether any of these genes overlap.
        //Now determine whether loci overlap. If so, combine these loci:
        // Allready assessed at the locus construction stage.  There should be no overlapping loci.
        depict.matrix.SymmetricShortDistanceMatrix matrixLoci = new depict.matrix.SymmetricShortDistanceMatrix(vecLoci.size());
        matrixLoci.setAllElements(matrixLoci.MAX_VALUE);
        for (int e = 0; e < nrGenes; e++) {
            if (vecLociPerGene[e].size() > 1) {
                //At least two loci use this gene, combine these loci:
                for (int l = 0; l < vecLociPerGene[e].size(); l++) {
                    int locusL = ((Integer) vecLociPerGene[e].get(l)).intValue();
                    for (int m = 0; m < vecLociPerGene[e].size(); m++) {
                        int locusM = ((Integer) vecLociPerGene[e].get(m)).intValue();
                        matrixLoci.set(locusL, locusM, 1);
                    }
                }
            }
        }
        for (int p = 0; p < matrixLoci.size(); p++) {
            matrixLoci.set(p, p, matrixLoci.MAX_VALUE);
        }
        matrixLoci.getAllPairsShortestPath();
        Vector vecUniqueLoci = new Vector();
        Vector vecUniqueLociNames = new Vector();
        boolean[] lociAlreadyCombined = new boolean[vecLoci.size()];
        for (int g = 0; g < vecLoci.size(); g++) {
            if (!lociAlreadyCombined[g]) {
                Vector vecGenes = (Vector) vecLoci.get(g);
                String uniqueLocusName = (String) vecLociNames.get(g);
                for (int h = 0; h < vecLoci.size(); h++) {
                    if (g != h && !lociAlreadyCombined[h] && matrixLoci.get(g, h) < 100) {
                        Vector vecGenesToCombine = (Vector) vecLoci.get(h);
                        uniqueLocusName += ";" + (String) vecLociNames.get(h);
                        for (int v = 0; v < vecGenesToCombine.size(); v++) {
                            if (!vecGenes.contains((String) vecGenesToCombine.get(v))) {
                                vecGenes.add(vecGenesToCombine.get(v));
                            }
                        }
                        vecLoci.set(h, null);
                        vecLoci.set(g, vecGenes);
                        lociAlreadyCombined[h] = true;
                    }
                }
                vecUniqueLoci.add(vecGenes);
                vecUniqueLociNames.add(uniqueLocusName);
            }
        }
        int nrUniqueLoci = vecUniqueLoci.size();
        System.out.println("Number of unique, non-overlapping loci:\t" + nrUniqueLoci);

        int nrUniqueGenesInLoci = 0;
        HashMap<String, Integer> genesInLoci = new HashMap<String, Integer>(); 
        for (int l = 0; l < nrUniqueLoci; l++) {
            Vector vecGenes = (Vector) vecUniqueLoci.get(l);
            nrUniqueGenesInLoci += vecGenes.size();
            Iterator it = vecGenes.iterator();
            while (it.hasNext()) {
                genesInLoci.put((String) it.next(),1);
            }
        }
        System.out.println("Number of unique genes in the unique, non-overlapping loci:\t" + nrUniqueGenesInLoci);
        
        //We now have unique, non-overlapping loci, now sort these loci from large to small loci, to speed up the permutation of the loci.
        //Sort the loci from large to small loci, when permuting loci, we first randomly pick genomic regions for the largest genes, and then pick random, non-overlapping regions for the smaller genes.
        //If we would do it in another way, and e.g. first define permuted loci for the small loci, it is likely that when many loci are included for analysis, some of the larger loci can not be assigned to a genomic position, as these will overlap with permuted smaller loci:
        Vector vecUniqueLociNrGenesPerLocus = new Vector();
        for (int l = 0; l < nrUniqueLoci; l++) {
            Vector vecGenes = (Vector) vecUniqueLoci.get(l);
            StringIntegerObject object = new StringIntegerObject("" + l, vecGenes.size());
            vecUniqueLociNrGenesPerLocus.add(object);
        }
        sorter.sort(vecUniqueLociNrGenesPerLocus);
        Vector vecUniqueLociSorted = new Vector();
        Vector vecUniqueLociNamesSorted = new Vector();
        for (int l = vecUniqueLociNrGenesPerLocus.size() - 1; l >= 0; l--) {
            StringIntegerObject object = (StringIntegerObject) vecUniqueLociNrGenesPerLocus.get(l);
            int index = Integer.parseInt(object.stringValue);
            vecUniqueLociSorted.add(vecUniqueLoci.get(index));
            vecUniqueLociNamesSorted.add(vecUniqueLociNames.get(index));
        }
        vecUniqueLoci = vecUniqueLociSorted;
        vecUniqueLociNames = vecUniqueLociNamesSorted;

        int[] nrUniqueLociEmpiricPerms = null;
        int[][] nrGenesPerUniqueLociEmpiricPerms = null;
        int[][] geneStartIndexPermutedLociEmpiricPerms = null;

        if (nrUniqueLoci < 10) {
            System.out.println("There are too few independent loci (n = " + nrUniqueLoci + "). DEPICT needs at least 10 independent loci.");
            System.exit(-1);
        }
        if (directoryPermutedLociDefinedBySignificantSNPssAndLDInformation != null) {

            Vector vecPermutedLoci = new Vector();

            File dir = new File(directoryPermutedLociDefinedBySignificantSNPssAndLDInformation);
            if (!dir.isDirectory()) {
                System.out.println("DEPICT failed to construct backgound loci for " + nrUniqueLoci + " GWAS loci. Please contact Tune H. Pers (tunepers@broadinstitute.org).");
                System.exit(-1);
            }
            File[] files = dir.listFiles();
            Vector vecValidPermFiles = new Vector();
            for (int f = 0; f < files.length; f++) {
                File file = files[f];
                String fileName = file.getName();
                if (fileName.contains("permutation")) {
                    vecValidPermFiles.add(file);
                }
            }
            for (int f = 0; f < vecValidPermFiles.size(); f++) {
                File file = (File) vecValidPermFiles.get(f);
                try {
                    //java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(file));
                    InputStream fileStream = new FileInputStream(file);
                    InputStream gzipStream = new GZIPInputStream(fileStream);
                    Reader decoder = new InputStreamReader(gzipStream, gzip_encoding);
                    BufferedReader in = new BufferedReader(decoder);       
                    String str = in.readLine();
                    while ((str = in.readLine()) != null) {
                        String[] data = str.replace(" ", "").split("\t");
                        if ( str.startsWith("#") || data[locusFileChrCol].equals(chrToBeLeftOut)) {
                            //Skip this comment line in locus definition file
                        } else {
                            Vector vecGenes = new Vector();
                            int startIndex = Integer.MAX_VALUE;
                            int endIndex = -1;

                            if (data.length == locusFileGenesCol + 1) {
                                String[] dataGenes = data[locusFileGenesCol].split(";");    
                                for (int d = 0; d < dataGenes.length; d++) {
                                    if (hashEnsemblIndex.containsKey(dataGenes[d])) {
                                        int ensemblIndex = ((Integer) hashEnsemblIndex.get(dataGenes[d])).intValue();
                                        if (ensemblIndex < startIndex) {
                                            startIndex = ensemblIndex;
                                        }
                                        if (ensemblIndex > endIndex) {
                                            endIndex = ensemblIndex;
                                        }
                                    }
                                }
                            } else if (data.length == locusFileNearestCol + 1) { //If there are no genes in the locus, take the closest gene
                                String[] dataNearest = data[locusFileNearestCol].split(";");
                                for (int e = 0; e < dataNearest.length; e++) {
                                    if (hashEnsemblIndex.containsKey(dataNearest[e])) {                            
                                        int ensemblIndex = ((Integer) hashEnsemblIndex.get(dataNearest[e])).intValue();
                                        if (ensemblIndex < startIndex) {
                                            startIndex = ensemblIndex;
                                        }
                                        if (ensemblIndex > endIndex) {
                                            endIndex = ensemblIndex;
                                        }
                                    }
                                }
                            }

                            if (startIndex != Integer.MAX_VALUE) {
                                int nrGenesThisLocus = endIndex - startIndex + 1;
                                depict.math.DoubleDoubleObject object = new depict.math.DoubleDoubleObject(startIndex, nrGenesThisLocus);
                                vecPermutedLoci.add(object);
                            }
                        }
                    }
                } catch (Exception e) {
                    System.out.println("Error with locus file:\t" + e.getMessage());
                    e.printStackTrace();
                }
            }

            System.out.println("Total number of permuted loci in all available permutation files:\t" + vecPermutedLoci.size());
            depict.math.DoubleDoubleObjectSorter sorter2 = new depict.math.DoubleDoubleObjectSorter();
            sorter2.sort(vecPermutedLoci);

            nrUniqueLociEmpiricPerms = new int[nrPerms + nrReps];
            nrGenesPerUniqueLociEmpiricPerms = new int[nrPerms + nrReps][1000];
            geneStartIndexPermutedLociEmpiricPerms = new int[nrPerms + nrReps][1000];

            for (int perm = 0; perm < nrPerms + nrReps; perm++) {
                nrUniqueLociEmpiricPerms[perm] = nrUniqueLoci;
            }

            boolean[][] geneAlreadyInUseByLocus = new boolean[nrPerms + nrReps][nrGenes];
            for (int l = 0; l < nrUniqueLoci; l++) {
                Vector vecGenesL = (Vector) vecUniqueLoci.get(l);
                Vector vecPermutedLociWithSimilarNumberOfGenes = new Vector();

                int nrRequestedGenes = vecGenesL.size();
                double lower = (double) nrRequestedGenes * 0.9300000000000001D;
                double upper = (double) nrRequestedGenes * 1.1D;
                while (true) {
                    for (int q = 0; q < vecPermutedLoci.size(); q++) {
                        depict.math.DoubleDoubleObject object = (depict.math.DoubleDoubleObject) vecPermutedLoci.get(q);
                        if (object.doubleValueToSortOn >= lower && object.doubleValueToSortOn <= upper) {
                            vecPermutedLociWithSimilarNumberOfGenes.add(object);
                        }
                    }
                    if (vecPermutedLociWithSimilarNumberOfGenes.size() >= 10) break;
                    lower *= 0.9300000000000001D;
                    upper *= 1.1D;
                    vecPermutedLociWithSimilarNumberOfGenes.clear();
                }
                System.out.println("Locus:\t" + l + "\tcontains:\t" + nrRequestedGenes + "\tgenes. Number of permuted loci from pool that have a similar number of genes:\t" + vecPermutedLociWithSimilarNumberOfGenes.size());
                for (int perm = 0; perm < nrPerms + nrReps; perm++) {
                    boolean validLocus = false;
                    while (!validLocus) {
                        int random = (int) (Math.random() * (double) vecPermutedLociWithSimilarNumberOfGenes.size());
                        depict.math.DoubleDoubleObject object = (depict.math.DoubleDoubleObject) vecPermutedLociWithSimilarNumberOfGenes.get(random);
                        int nrGenesThisLocus = (int) object.doubleValueToSortOn;
                        int startIndexThisLocus = (int) object.doubleValue;
                        validLocus = true;
                        for (int g = 0; g < nrGenesThisLocus; g++) {
                            if (geneAlreadyInUseByLocus[perm][g + startIndexThisLocus]) {
                                validLocus = false;
                                break;
                            }
                        }
                        if (validLocus) {
                            nrGenesPerUniqueLociEmpiricPerms[perm][l] = nrGenesThisLocus;
                            geneStartIndexPermutedLociEmpiricPerms[perm][l] = startIndexThisLocus;
                            for (int g = 0; g < nrGenesThisLocus; g++) {
                                geneAlreadyInUseByLocus[perm][g + startIndexThisLocus] = true;
                            }
                        }
                    }
                }
            }
        }

        //Define HashMap with all the genes in the real loci:
        HashMap hashGenesInRealLoci = new HashMap();
        for (int l = 0; l < nrUniqueLoci; l++) {
            Vector vecGenesL = (Vector) vecUniqueLoci.get(l);
            for (int vl = 0; vl < vecGenesL.size(); vl++) {
                hashGenesInRealLoci.put(vecGenesL.get(vl), null);
            }
        }


        //Store the most significant prioritized gene per locus in the real analysis:
        Vector vecResultsRealAnalysis = new Vector();
        //Store the most significant prioritized gene per locus in the different repetitions of the algorithm:
	int nrRepsYDim = 5000;
        double[][] repPValues = new double[nrReps][nrRepsYDim];
	for(int i=0; i<nrReps; i++) {
		for(int j=0; j<nrRepsYDim; j++) {
			repPValues[i][j] = 1.01;
		}
	}
        // FDR for genes outside loci
	Vector vecResultsOutsideRealAnalysis = new Vector();
	int nrRepsOutsideYDim = 20000;
        double[][] repPValuesOutside = new double[nrReps][nrRepsOutsideYDim];
	for(int i=0; i<nrReps; i++) {
		for(int j=0; j<nrRepsOutsideYDim; j++) {
			repPValuesOutside[i][j] = 1.01;
		}
	}

        //Now start with the real analysis (rep = -1), and then run the repetitions (necessary for the eventual FDR calculation):
        Vector vecResultsRealPathwayAnalysis = new Vector();
        double[][] repPathwayPValues = new double[nrReps][dataset.nrSamples];

        // Anntations for pathway analysis
        HashMap hashPathwayDescription = new HashMap();
        if(conductTissueAnalysis || conductPathwayAnalysis) {
           
            // Read annotation
            try {
                java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(goMappingFile)));
                String str;
                while ((str = in.readLine()) != null) {
                    String[] data = str.split("\t");
                    hashPathwayDescription.put(data[0], data[2]);
                }
                in.close();
                in = new java.io.BufferedReader(new java.io.FileReader(new File(mgiMappingFile)));
                while ((str = in.readLine()) != null) {
                    String[] data = str.split("\t");
                    hashPathwayDescription.put(data[0], data[1]);
                }
                in.close();
                in = new java.io.BufferedReader(new java.io.FileReader(new File(inwebMappingFile)));
                while ((str = in.readLine()) != null) {
                    String[] data = str.split("\t");
                    hashPathwayDescription.put(data[0], data[1]);
                }
                in.close();
                in = new java.io.BufferedReader(new java.io.FileReader(new File(tissueMappingFile)));
                while ((str = in.readLine()) != null) {
                    String[] data = str.split(";");
                    //if (data.length > 4) {
                    //    hashPathwayDescription.put(data[0], data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4]);
                    //} else {
                    hashPathwayDescription.put(data[0], data[1]);
                    //}
                }
                in.close();
            
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }             
        }            
        
        java.util.concurrent.ExecutorService threadPool = Executors.newFixedThreadPool(nrCores);
        CompletionService<IntegerDoubleObject> completionService = new ExecutorCompletionService<IntegerDoubleObject>(threadPool);        
        for (int rep = -1; rep < nrReps; rep++) {

            //Store for each genes the prioritization z-scores, we will use this later on to assess whether individual GO-terms, MGI-phenotypes or pathways in the database strongly correlate with this z-score vector:
            double[] allGenesPrioritizationZScore = new double[nrGenes];

            if (rep == -1) {
                System.out.println("Running analysis of real data");
                //System.out.println("\n\n\n");
            } else {
                System.out.println("Running repetition:\t" + (rep + 1));
            }

            //If we start working on the repetitions (to eventually get accurate FDR estimates), we permute the loci:
            if (rep >= 0) {

                Vector vecUniqueLociPerm = new Vector();
                if (directoryPermutedLociDefinedBySignificantSNPssAndLDInformation == null) {

                    //Permute the loci for this repetition:
                    int[] nrGenesPerLocus = new int[nrUniqueLoci];
                    for (int l = 0; l < nrUniqueLoci; l++) {
                        Vector vecGenes = (Vector) vecUniqueLoci.get(l);
                        nrGenesPerLocus[l] = vecGenes.size();
                    }
                    Arrays.sort(nrGenesPerLocus);
                    int[] nrGenesPerLocusAsc = new int[nrUniqueLoci];
                    for (int n = 0; n < nrUniqueLoci; n++) {
                        nrGenesPerLocusAsc[n] = nrGenesPerLocus[nrUniqueLoci - n - 1];
                    }
                    nrGenesPerLocus = nrGenesPerLocusAsc;

                    if (hashClosestGeneToSNP != null) {
                        hashClosestGeneToSNP = new HashMap();
                    }

                    boolean[] genePartPermutedLocus = new boolean[nrGenes];
                    //Process each unique loci, start with the largest locus:
                    for (int l = 0; l < nrUniqueLoci; l++) {
                        int nrGenesThisLocus = nrGenesPerLocus[l];
                        //Assign a random position on the genome, check whether there is not yet a permuted locus overlapping this position:
                        boolean validPermutedLocus = false;
                        while (!validPermutedLocus) {

                            //Get a random start position for this locus. This procedure however does takes into account the length of genes. If we want to correct for gene size this method preferentially will include larger. This is proportional to the gene length of each gene:
                            int randomPosStart = randomGene[(int) (Math.random() * (double) totalLengthGenes)];
                            randomPosStart -= nrGenesThisLocus / 2;
                            if (randomPosStart < 0) {
                                randomPosStart = 0;
                            }
                            boolean locusFree = true;

                            //Check whether this position is available. I.e. ensure this requested position is not overlapping with an already permuted locus:
                            for (int i = randomPosStart; i < randomPosStart + nrGenesThisLocus; i++) {
                                if (i >= nrGenes || genePartPermutedLocus[i]) {
                                    locusFree = false;
                                    break;
                                }
                            }

                            //We are save, this locus is not overlapping, so defined this as the permuted position:
                            if (locusFree) {
                                validPermutedLocus = true;
                                Vector vecGenes = new Vector();
                                for (int i = randomPosStart; i < randomPosStart + nrGenesThisLocus; i++) {
                                    genePartPermutedLocus[i] = true;
                                    vecGenes.add((String) vecEnsemblIndex.get(i));
                                }
                                if (hashClosestGeneToSNP != null) {
                                    hashClosestGeneToSNP.put(vecEnsemblIndex.get(randomPosStart), null);
                                }
                                vecUniqueLociPerm.add(vecGenes);
                            }
                        }
                    }

                } else {

                    for (int l = 0; l < nrUniqueLociEmpiricPerms[nrPerms + rep]; l++) {
                        Vector vecGenes = new Vector();
                        for (int i = geneStartIndexPermutedLociEmpiricPerms[nrPerms + rep][l]; i < geneStartIndexPermutedLociEmpiricPerms[nrPerms + rep][l] + nrGenesPerUniqueLociEmpiricPerms[nrPerms + rep][l]; i++) {
                            vecGenes.add((String) vecEnsemblIndex.get(i));
                        }
                        vecUniqueLociPerm.add(vecGenes);
                    }
                    if (hashClosestGeneToSNP != null) {
                        hashClosestGeneToSNP = new HashMap();
                    }

                }

                //Done, the loci have been permuted:
                vecUniqueLoci = vecUniqueLociPerm;

            }

            int[] nrGenesPerUniqueLoci = new int[nrUniqueLoci]; //Score the number of genes per locus:
            String[] descriptionPerUniqueLoci = new String[nrUniqueLoci]; //Make a description of the locus
            int[] distanceGeneIndexToClosestDefinedLocus = new int[nrGenes]; //For genes outside the defined loci, score the distance to the closest defined locus
            int[] distanceGeneIndexToClosestDefinedLocusInBasePairs = new int[nrGenes]; //For genes outside the defined loci, score the distance to the closest defined locus
            int[] distanceGeneIndexToClosestDefinedLocusLocusIndex = new int[nrGenes]; //For genes outside the defined loci, score the locus ID for the closest defined locus
            for (int g = 0; g < nrGenes; g++) {
                distanceGeneIndexToClosestDefinedLocus[g] = Integer.MAX_VALUE;
                distanceGeneIndexToClosestDefinedLocusInBasePairs[g] = Integer.MAX_VALUE;
                distanceGeneIndexToClosestDefinedLocusLocusIndex[g] = -1;
            }

            //Now make a description for every locus: Include the chromosome, the start position and end position of each locus:
            //Also determine for every gene how far away that gene is mapping from a predefined locus:
            for (int l = 0; l < nrUniqueLoci; l++) {
                Vector vecGenesL = (Vector) vecUniqueLoci.get(l);
                nrGenesPerUniqueLoci[l] = vecGenesL.size();
                int minPos = Integer.MAX_VALUE;
                int maxPos = Integer.MIN_VALUE;
                int chr = -1;
                for (int vl = 0; vl < vecGenesL.size(); vl++) {
                    int geneIndex = ((Integer) hashEnsemblIndex.get((String) vecGenesL.get(vl))).intValue();
                    chr = ensemblChr[geneIndex];
                    if (ensemblChrPosStart[geneIndex] < minPos) {
                        minPos = ensemblChrPosStart[geneIndex];
                    }
                    if (ensemblChrPosStart[geneIndex] > maxPos) {
                        maxPos = ensemblChrPosStart[geneIndex];
                    }
                    if (ensemblChrPosEnd[geneIndex] < minPos) {
                        minPos = ensemblChrPosEnd[geneIndex];
                    }
                    if (ensemblChrPosEnd[geneIndex] > maxPos) {
                        maxPos = ensemblChrPosEnd[geneIndex];
                    }
                    for (int g = 0; g < nrGenes; g++) {
                        if (chr == ensemblChr[g]) {
                            int distance = Math.abs(geneIndex - g);
                            if (distance < distanceGeneIndexToClosestDefinedLocus[g]) {
                                distanceGeneIndexToClosestDefinedLocus[g] = distance;
                                distanceGeneIndexToClosestDefinedLocusLocusIndex[g] = l;
                            }
                            int distanceInBasePairs = Math.abs(ensemblChrPosStart[geneIndex] - ensemblChrPosStart[g]);
                            if (distanceInBasePairs < distanceGeneIndexToClosestDefinedLocusInBasePairs[g]) {
                                distanceGeneIndexToClosestDefinedLocusInBasePairs[g] = distanceInBasePairs;
                            }
                            distanceInBasePairs = Math.abs(ensemblChrPosStart[geneIndex] - ensemblChrPosEnd[g]);
                            if (distanceInBasePairs < distanceGeneIndexToClosestDefinedLocusInBasePairs[g]) {
                                distanceGeneIndexToClosestDefinedLocusInBasePairs[g] = distanceInBasePairs;
                            }
                            distanceInBasePairs = Math.abs(ensemblChrPosEnd[geneIndex] - ensemblChrPosStart[g]);
                            if (distanceInBasePairs < distanceGeneIndexToClosestDefinedLocusInBasePairs[g]) {
                                distanceGeneIndexToClosestDefinedLocusInBasePairs[g] = distanceInBasePairs;
                            }
                            distanceInBasePairs = Math.abs(ensemblChrPosEnd[geneIndex] - ensemblChrPosEnd[g]);
                            if (distanceInBasePairs < distanceGeneIndexToClosestDefinedLocusInBasePairs[g]) {
                                distanceGeneIndexToClosestDefinedLocusInBasePairs[g] = distanceInBasePairs;
                            }
                        }
                    }
                }
                descriptionPerUniqueLoci[l] = "chr" + chr + ":" + minPos + "-" + maxPos;
            }

            if ( ( conductTissueAnalysis || conductPathwayAnalysis ) && repPathwayPValues != null) {
                
                //Process each gene in this locus l:
                double[][] meanValPerLocus = new double[dataset.nrSamples][nrUniqueLoci];
                for (int l = 0; l < nrUniqueLoci; l++) {
                  
                    //Get all genes for locus l:
                    Vector vecGenesL = (Vector) vecUniqueLoci.get(l);
                    for (int vl = 0; vl < vecGenesL.size(); vl++) {
                    
                        //Get the gene index:
                        int geneL = ((Integer) dataset.hashProbes.get((String) vecGenesL.get(vl))).intValue();
                        for (int s = 0; s < dataset.nrSamples; s++) {
                            double zScore = dataset.rawData[geneL][s];
                            meanValPerLocus[s][l] += zScore;
                        }
                    }
                    for (int s = 0; s < dataset.nrSamples; s++) {
                        meanValPerLocus[s][l] /= (double) vecGenesL.size();
                    }
                }


                cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.DRand();
                for (int s = 0; s < dataset.nrSamples; s++) {

                    double realZScore = mean(meanValPerLocus[s]) / JSci.maths.ArrayMath.standardDeviation(meanValPerLocus[s]);

                    double[] zScoresPerms = new double[nrPerms];

                    for (int perm = 0; perm < nrPerms; perm++) {
                            PerformPathwayAnalysisSusceptibilityLociPermutationTask task = new PerformPathwayAnalysisSusceptibilityLociPermutationTask(perm, s, nrUniqueLociEmpiricPerms[perm], nrGenesPerUniqueLociEmpiricPerms[perm], geneStartIndexPermutedLociEmpiricPerms, geneIDToMatrixID, nrGenes, dataset);
                            completionService.submit(task);
                    }
                    for (int q = 0; q < nrPerms; q++) {
                        try {
                            IntegerDoubleObject result = completionService.take().get();
                            int perm = result.intValue;
                            zScoresPerms[perm] = result.doubleValue;
                        } catch (ExecutionException ex) {
                            Logger.getLogger(PerformPathwayAnalysisSusceptibilityLoci.class.getName()).log(Level.SEVERE, null, ex);
                        }

                    }

                    //Calculate the mean and stdev of the minLogPValueSum for the permutations:
                    double mean = JSci.maths.ArrayMath.mean(zScoresPerms);
                    double stdev = JSci.maths.ArrayMath.standardDeviation(zScoresPerms);

                    //Correct the real observed minLogPValueSum for the permutations mean and stdev.
                    //This procedure assumes a normal distribution of these values. I have not thoroughly tested this, but I think we are safe here, exploiting the central limit theorem:
                    double zScoreCorrectedForPermutedZScoreDist = (realZScore - mean) / stdev;
                    double pValue = cern.jet.stat.Probability.normal(-zScoreCorrectedForPermutedZScoreDist);
                    if (pValue < 2.0E-323) {
                        pValue = 2.0E-323;
                    }
                    
                    if (rep == -1) {
                        //Real analysis, both store the co-functionality P-Value and additional annotation information:                        
                        String description = dataset.sampleNames[s];
                        if (hashPathwayDescription.containsKey(description)) description = (String) hashPathwayDescription.get(description);
                        String pValueNice = pValue < 0.01 ? String.format("%6.2e",pValue) : String.format( "%.2f",pValue);
                        String output = dataset.sampleNames[s] + "\t" + description + "\t" + pValueNice;
                        if (conductTissueAnalysis) {
                            output = description + "\t" + pValueNice;
                        }
                        StringDoubleObject stringDoubleObject = new StringDoubleObject(output, pValue);
                        vecResultsRealPathwayAnalysis.add(stringDoubleObject);
                    } else {
                        //Results of the repetitions, only store the co-functionality P-Value, do not store additional annotation information:
                        repPathwayPValues[rep][s] = pValue;
                    }
                }
            }

            if (conductNetworkAnalysis) {

                // Prioritize genes in associated loci
		int geneItr = 0;
                for (int l = 0; l < nrUniqueLoci; l++) {

                    //Get all genes for locus l:
                    Vector vecGenesL = (Vector) vecUniqueLoci.get(l);
                    double maxZScoreThisLocus = -100000d;
                    String maxZScoreThisLocusGene = "-";

                    //Process each gene in this locus l:
                    for (int vl = 0; vl < vecGenesL.size(); vl++) {

                        //Get the gene index:
                        int geneL = ((Integer) dataset.hashProbes.get((String) vecGenesL.get(vl))).intValue();

                        //Calculate the co-functionality score for geneL with all other other loci:
                        double minLogPValueSum = 0;
                        double sumZScores = 0;

                        //Ascertain all other loci (m):
                        for (int m = 0; m < nrUniqueLoci; m++) {
                            if (l != m) {

                                //Get all genes for locus m:
                                Vector vecGenesM = (Vector) vecUniqueLoci.get(m);

                                //What is the P-Value of the gene in locus m that has the highest co-functionality with geneL:
                                double pValueMin = 1;

                                //Process each gene in this locus m:
                                for (int vm = 0; vm < nrGenesPerUniqueLoci[m]; vm++) {

                                    //Get the gene index:
                                    int geneM = ((Integer) dataset.hashProbes.get((String) vecGenesM.get(vm))).intValue();

                                    //Calculate the P-Value for the co-functionality between geneL and geneM:
                                    double pValue = 1;
                                    //If we have not enforced a uniform P-Value distribution on the co-functionality network, get the original P-Value, otherwise get the corrected P-Value:
                                    if (matrixUniform == null) {
                                        pValue = matrixValToPValue[matrix.get(geneL, geneM)];
                                    } else {
                                        pValue = matrixValToPValue[(int) matrixUniform[geneL][geneM]];
                                    }

                                    //If geneM has a higher co-functionality with geneL, store this minimal P-Value:
                                    if (pValue < pValueMin) {
                                        pValueMin = pValue;
                                    }

                                }

                                //We now have the most significant P-Value for geneL with genes in locus m. Correct for the number of ascertained genes in locus m:
                                pValueMin = 1d - Math.pow(1d - pValueMin, (double) nrGenesPerUniqueLoci[m]);
                                //There is a limit on the accuracy of this exponential function:
                                if (pValueMin < 1E-16) {
                                    pValueMin = 1E-16;
                                }

                                //Now increase the score of gene l for locus m:
                                minLogPValueSum += -Math.log(pValueMin);
                                if (nrPerms == 0) {
                                    sumZScores += cern.jet.stat.Probability.normalInverse(pValueMin);
                                }

                            }
                        }

                        double zScore = 0;
                        if (nrPerms == 0) {

                            //If we do not want to run permutations, use an unweighted Z-Score method:
                            zScore = -sumZScores / Math.sqrt((double) nrUniqueLoci);

                        } else {

                            //Now do exactly the same for the permuted loci.
                            //Ascertain for geneL what the minLogPValueSum would be, when permuting all loci, except for locus l:
                            //This is multithreaded code, so it looks a bit more complex, but this does exactly the same as the above procedure for the real loci:
                            double[] minLogPValueSumPerms = new double[nrPerms];
                            for (int perm = 0; perm < nrPerms; perm++) {
                                    PerformNetworkAnalysisSusceptibilityLociPermutationTask task = new PerformNetworkAnalysisSusceptibilityLociPermutationTask(perm, l, vl, geneL, nrUniqueLociEmpiricPerms[perm], nrGenesPerUniqueLociEmpiricPerms[perm], geneStartIndexPermutedLociEmpiricPerms, geneIDToMatrixID, nrGenes, matrix, matrixValToPValue, matrixUniform);
                                    completionService.submit(task);
                            }
                            for (int q = 0; q < nrPerms; q++) {
                                try {
                                    IntegerDoubleObject result = completionService.take().get();
                                    int perm = result.intValue;
                                    //Store for each of the permutations what the minLogPValueSum is:
                                    minLogPValueSumPerms[perm] = result.doubleValue;
                                } catch (ExecutionException ex) {
                                    Logger.getLogger(PerformPathwayAnalysisSusceptibilityLoci.class.getName()).log(Level.SEVERE, null, ex);
                                }

                            }

                            //We now have the real minLogPValueSum for geneL, and also have it for the permutations:

                            //Calculate the mean and stdev of the minLogPValueSum for the permutations:
                            double mean = JSci.maths.ArrayMath.mean(minLogPValueSumPerms);
                            double stdev = JSci.maths.ArrayMath.standardDeviation(minLogPValueSumPerms);

                            //Correct the real observed minLogPValueSum for the permutations mean and stdev.
                            //This procedure assumes a normal distribution of these values. I have not thoroughly tested this, but I think we are safe here, exploiting the central limit theorem:
                            zScore = (minLogPValueSum - mean) / stdev;
                        }

                        //Now calculate the single-sided P-Value: The higher the minLogPValue in the real analysis compared to the permutations, the more significant this gene is:
                        double pValue = cern.jet.stat.Probability.normal(-zScore);
                        if (pValue < 2.0E-323) {
                            pValue = 2.0E-323;
                        }

                        //Store the Z-Score in a vector, permitting subsequent individual pathway enrichment:
                        int geneIndex = ((Integer) hashEnsemblIndex.get((String) vecGenesL.get(vl))).intValue();
                        allGenesPrioritizationZScore[geneIndex] = zScore;

                        if (rep == -1) {

                            //In the real analysis, output for every included gene various statistics:
                            String ensemblGene = (String) vecGenesL.get(vl);
                            int geneSize = ensemblChrPosEnd[geneIndex] - ensemblChrPosStart[geneIndex];

                            // eQTL information
                            String topCisEQTLSNP = "-";
                            if (hashEnsemblGeneCisEQTLSNP.containsKey(ensemblGene)) {
                                topCisEQTLSNP = (String) hashEnsemblGeneCisEQTLSNP.get(ensemblGene);
                            }

                            // Nearest gene
                            String nearestGene = "-";
                            if (hashClosestGeneToSNP != null) {
                                    nearestGene = "" + hashClosestGeneToSNP.containsKey((String) vecGenesL.get(vl)); 
                            }

                            // Write output string
                            String hgnc = (String) hashEnsemblToHGNC.get((String) vecGenesL.get(vl));
                            String pValueNice = pValue < 0.01 ? String.format("%6.2e",pValue) : String.format( "%.2f",pValue);
                            String output = (String) vecUniqueLociNames.get(l) + "\t" + nrGenesPerUniqueLoci[l] + "\t" + descriptionPerUniqueLoci[l] + "\t" + hashClosestGeneToSNP.get((String) vecUniqueLociNames.get(l)) + "\t" + (String) vecGenesL.get(vl) + "\t" + hgnc + "\t" + pValueNice + "\t" + nearestGene + "\t" + topCisEQTLSNP;

                            //Real analysis, both store the co-functionality P-Value and additional annotation information:
                            StringDoubleObject stringDoubleObject = new StringDoubleObject(output, pValue);
                            vecResultsRealAnalysis.add(stringDoubleObject);
			} else {
				//Results of the repetitions, only store the co-functionality P-Value, do not store additional annotation information:
				repPValues[rep][geneItr++] = pValue;
			}
		   }
                }
                
                // Prioritize genes outside associated loci
                if (calculateGenePrioritizationPValueForGenesOutsideLoci) {
                    //We have now calculate the gene co-functionality P-Value for the genes inside the loci, now also calculate this for all the genes outside the loci:
                    //The procedure is identical to the above procedure:

                    //Now determine which genes map outside the loci:
                    Vector vecGenesOutsideLoci = new Vector();
                    for (int p = 0; p < dataset.nrProbes; p++) {
                        if (!genesInLoci.containsKey(dataset.probeNames[p]) && hashEnsemblIndex.containsKey(dataset.probeNames[p])) {
                            vecGenesOutsideLoci.add(dataset.probeNames[p]);
                        }
                    }

                    //Process all genes outside the loci:
                    HashMap hashCofunctionalityPValueGenesOutsideLoci = new HashMap();
                    for (int vl = 0; vl < vecGenesOutsideLoci.size(); vl++) {
                        int geneL = ((Integer) dataset.hashProbes.get((String) vecGenesOutsideLoci.get(vl))).intValue();
                        double minLogPValueSum = 0;
                        double sumZScores = 0;
                        for (int m = 0; m < nrUniqueLoci; m++) {
                            Vector vecGenesM = (Vector) vecUniqueLoci.get(m);
                            double pValueMin = 1;
                            for (int vm = 0; vm < nrGenesPerUniqueLoci[m]; vm++) {
                                int geneM = ((Integer) dataset.hashProbes.get((String) vecGenesM.get(vm))).intValue();
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
                            pValueMin = pValueMin;
                            if (pValueMin < 1E-16) {
                                pValueMin = 1E-16;
                            }
                            minLogPValueSum += -Math.log(pValueMin);
                            if (nrPerms == 0) {
                                sumZScores += cern.jet.stat.Probability.normalInverse(pValueMin);
                            }
                        }

                        double zScore = 0;
                        if (nrPerms == 0) {

                            //If we do not want to run permutations, use an unweighted Z-Score method:
                            zScore = -sumZScores / Math.sqrt((double) nrUniqueLoci);

                        } else {
                            for (int perm = 0; perm < nrPerms; perm++) {
                                    PerformNetworkAnalysisSusceptibilityLociPermutationTask task = new PerformNetworkAnalysisSusceptibilityLociPermutationTask(perm, -1, vl, geneL, nrUniqueLociEmpiricPerms[perm], nrGenesPerUniqueLociEmpiricPerms[perm], geneStartIndexPermutedLociEmpiricPerms, geneIDToMatrixID, nrGenes, matrix, matrixValToPValue, matrixUniform);
                                    completionService.submit(task);
                            }

                            double[] minLogPValueSumPerms = new double[nrPerms];
                            for (int q = 0; q < nrPerms; q++) {
                                try {
                                    IntegerDoubleObject result = completionService.take().get();
                                    int perm = result.intValue;
                                    minLogPValueSumPerms[perm] = result.doubleValue;
                                } catch (ExecutionException ex) {
                                    Logger.getLogger(PerformPathwayAnalysisSusceptibilityLoci.class.getName()).log(Level.SEVERE, null, ex);
                                }

                            }

                            double mean = JSci.maths.ArrayMath.mean(minLogPValueSumPerms);
                            double stdev = JSci.maths.ArrayMath.standardDeviation(minLogPValueSumPerms);
                            zScore = (minLogPValueSum - mean) / stdev;
                        }

                        int geneIndex = ((Integer) hashEnsemblIndex.get((String) vecGenesOutsideLoci.get(vl))).intValue();
                        allGenesPrioritizationZScore[geneIndex] = zScore;

                        double pValue = cern.jet.stat.Probability.normal(-zScore);
                        if (pValue < 2.0E-323) {
                       		 pValue = 2.0E-323;
                        }

                        if (rep == -1) {
                            String ensemblGene = (String) vecGenesOutsideLoci.get(vl);
                            String topCisEQTLSNP = "-";
                            if (hashEnsemblGeneCisEQTLSNP.containsKey(ensemblGene)) {
                                topCisEQTLSNP = (String) hashEnsemblGeneCisEQTLSNP.get(ensemblGene);
                            }
                            String output = vl + "\t" + ensemblGene + "\t" + ensemblChr[geneIndex] + "\t" + ensemblChrPosStart[geneIndex] + "\t" + ensemblChrPosEnd[geneIndex] + "\t" + pValue + "\t" + distanceGeneIndexToClosestDefinedLocusLocusIndex[geneIndex] + "\t" + distanceGeneIndexToClosestDefinedLocus[geneIndex] + "\t" + distanceGeneIndexToClosestDefinedLocusInBasePairs[geneIndex] + "\t" + topCisEQTLSNP;

                            if (distanceGeneIndexToClosestDefinedLocusInBasePairs[geneIndex] > 500000) {
                                hashCofunctionalityPValueGenesOutsideLoci.put(ensemblGene, pValue);
                            }

       		   	    //Real analysis, both store the co-functionality P-Value and additional annotation information:
                       	    StringDoubleObject stringDoubleObject = new StringDoubleObject(output, pValue);
	                    vecResultsOutsideRealAnalysis.add(stringDoubleObject);

                        } else {
				//Results of the repetitions, only store the co-functionality P-Value, do not store additional annotation information:
				repPValuesOutside[rep][vl] = pValue;
			}
                    }

                    /*
                    if (rep == -1 && vecGenePValuesOutsideLociNotUsedForPrioritization != null) {

                        //Take all the prioritized genes outside the loci, and make QQ plots in bins of 1000 genes, sorted on the prioritization P-Value, do we see enrichment of signals outside the significant GWAS genes?
                        Vector vectorGWASGenePValuePrioritizationPValue = new Vector();
                        for (int g = 0; g < vecGenePValuesOutsideLociNotUsedForPrioritization.size(); g++) {
                            String ensembl = (String) vecGenePValuesOutsideLociNotUsedForPrioritization.get(g);
                            if (hashCofunctionalityPValueGenesOutsideLoci.containsKey(ensembl)) {
                                double gwasGenePValue = ((Double) hashGenePValuesOutsideLociNotUsedForPrioritization.get(ensembl)).doubleValue();
                                double prioritizationPValue = ((Double) hashCofunctionalityPValueGenesOutsideLoci.get(ensembl)).doubleValue();
                                depict.math.DoubleDoubleObject object = new depict.math.DoubleDoubleObject(gwasGenePValue, prioritizationPValue);
                                vectorGWASGenePValuePrioritizationPValue.add(object);
                            }
                        }
                        depict.math.DoubleDoubleObjectSorter doubleDoubleObjectsorter = new depict.math.DoubleDoubleObjectSorter();
                        doubleDoubleObjectsorter.sort(vectorGWASGenePValuePrioritizationPValue);


                        int binIncrement = 500;
                        System.out.println("\n\nQQPLot of GWAS Gene PValues of genes with decreasing prioritization signficance:\n");
                        System.out.println("Bin\tNrGenesPerBin\tMedianChiSquareObserved\tMedianChiSquareExpected\tLambdaInflation");
                        for (int bin = 0; bin < vectorGWASGenePValuePrioritizationPValue.size(); bin += binIncrement) {
                            Vector vecObs = new Vector();
                            for (int a = 0; a < binIncrement; a++) {
                                if (bin + a < vectorGWASGenePValuePrioritizationPValue.size()) {
                                    depict.math.DoubleDoubleObject doubleDoubleObject = (depict.math.DoubleDoubleObject) vectorGWASGenePValuePrioritizationPValue.get(bin + a);
                                    vecObs.add(doubleDoubleObject.doubleValue);
                                }
                            }
                            if (vecObs.size() > 1) {
                                double[] valsObs = new double[vecObs.size()];
                                double[] valsExp = new double[vecObs.size()];
                                double startExpP = 0.5d / (double) vecObs.size();
                                JSci.maths.statistics.ChiSqrDistribution chiSqrDistribution = new JSci.maths.statistics.ChiSqrDistribution(1);
                                for (int v = 0; v < valsObs.length; v++) {
                                    valsObs[v] = ((Double) vecObs.get(v)).doubleValue();
                                    valsExp[v] = startExpP;
                                    valsObs[v] = chiSqrDistribution.inverse(1 - valsObs[v]);
                                    valsExp[v] = chiSqrDistribution.inverse(1 - valsExp[v]);
                                    startExpP += 1.0d / (double) vecObs.size();
                                }
                                Arrays.sort(valsObs);
                                Arrays.sort(valsExp);
                                double medianObs = cern.jet.stat.Descriptive.median(new cern.colt.list.DoubleArrayList(valsObs));
                                double medianExp = cern.jet.stat.Descriptive.median(new cern.colt.list.DoubleArrayList(valsExp));
                                double lambda = medianObs / medianExp;
                                System.out.println(bin + "\t" + vecObs.size() + "\t" + medianObs + "\t" + medianExp + "\t" + lambda);
                            }
                        }        
                    }

                    if (rep == -1) {
                        System.out.println("\n\n\n");
                    }
                    */    
                }
            }

        }
        threadPool.shutdown();

        if ( (conductTissueAnalysis || conductPathwayAnalysis ) && repPathwayPValues != null) {
            StringDoubleObjectSorter stringDoubleObjectSorter = new StringDoubleObjectSorter();
            stringDoubleObjectSorter.sort(vecResultsRealPathwayAnalysis);

            double[] pValues = new double[dataset.nrSamples];
            for (int l = 0; l < dataset.nrSamples; l++) {
                StringDoubleObject stringDoubleObject = (StringDoubleObject) vecResultsRealPathwayAnalysis.get(l);
                pValues[l] = stringDoubleObject.doubleValue;
            }

            for (int rep = 0; rep < nrReps; rep++) {
                Arrays.sort(repPathwayPValues[rep]);
            }

	    // Compute FDR and meta gene sets
            double maxNrTPs = 0;
            try {
                String outFile = outputGenesetFileName;
                if (conductTissueAnalysis) {
                    outFile = outputTissueFileName;
                }
                java.io.BufferedWriter out = new java.io.BufferedWriter(new java.io.FileWriter(new File(outFile)));
                
                
                // Header for gene set enrichment and tissue enrichment result files
                String genesInGeneSetHeader = "";
                for (int i = 0; i < Math.min(maxTopGenesPerGeneSet,limitMaxTopGenesPerGeneSet); i++) { 
                    if (!conductTissueAnalysis) {
                        genesInGeneSetHeader += "\tReconstituted gene set Z score gene " + String.valueOf(i+1); 
                    } else {
                        genesInGeneSetHeader += "\tTissue-specific expression Z score gene " + String.valueOf(i+1); 
                    }
                }
                if (conductTissueAnalysis) {
                    if (cofuncMatrixFile.toLowerCase().contains(gtex_identifier)) {
                        out.write("GTEx tissue\tNominal P value\tFalse discovery rate" + new String(genesInGeneSetHeader) + "\n");
                    } else {
                        out.write("MeSH term\tName\tMeSH first level term\tMeSH second level term\tNominal P value\tFalse discovery rate" + new String(genesInGeneSetHeader) + "\n");
                    }
                }  else {
                    out.write("Original gene set ID\tOriginal gene set description\tNominal P value\tFalse discovery rate" + new String(genesInGeneSetHeader) + "\n");
                }             
                
                // Loop over gene sets
                boolean fdrBinary_01 = true;
                boolean fdrBinary_05 = true;
                boolean fdrBinary_20 = true;                
                for (int i = 0; i < dataset.nrSamples; i++) {

                    StringDoubleObject stringDoubleObject = (StringDoubleObject) vecResultsRealPathwayAnalysis.get(i);
                    String output = stringDoubleObject.stringValue;
                    
                    // Top genes for gene set / tissue
                    int geneSetIndex = geneSetIDToMatrixID.get(output.split("\t")[0]);                    
                    Vector geneSetGenes = new Vector();
                    for (int j = 0; j < dataset.nrProbes; j++) { // Iterate over genes
                        geneSetGenes.add(new StringDoubleObject(dataset.probeNames[j], dataset.rawData[j][geneSetIndex])); // Save Z scores of all genes
                    }
                    StringDoubleObjectSorter mySorter = new StringDoubleObjectSorter(); // Sort genes according to Z score (low to high)
                    mySorter.sort(geneSetGenes);
                    int counter = 0;
                    ArrayList<String> geneSetGenesInLocus = new ArrayList<String>();
                    for (int j = (dataset.nrProbes-1); j >= 0; j--) { // Iterate over genes for gene set / tissue (high to low Z score)
                        StringDoubleObject geneAndZ = (StringDoubleObject) geneSetGenes.get(j);
                        if (counter < Math.min(maxTopGenesPerGeneSet,limitMaxTopGenesPerGeneSet)) {
                            if (genesInLoci.containsKey(geneAndZ.stringValue)) { // Limit to genes in associated loci
                                String geneID = (String) hashEnsemblToHGNC.get(geneAndZ.stringValue); // Get gene symbol
                                if (geneID.equals("-")) {
                                    geneID = geneAndZ.stringValue; // Use Ensembl identifier
                                }
                                geneSetGenesInLocus.add(geneID + " (" + String.format("%.1f", geneAndZ.doubleValue) + ")");
                                counter++;
                            }
                        } else {
                            break; // Maxium number of top genes found for gene set
                        }
                    } 
                    
                    // FDR
                    double fdr = 0;
                    for (int rep = 0; rep < nrReps; rep++) {
                        for (int m = 0; m < dataset.nrSamples; m++) {
                            if (repPathwayPValues[rep][m] <= pValues[i]) {
                                fdr++;
                            }
                        }
                    }
                    fdr /= (double) nrReps;
                    fdr /= (double) (i + 1);
                    double nrTPs = (1 - fdr) * (double) (i + 1);
                    if (nrTPs > maxNrTPs) {
                        maxNrTPs = nrTPs;
                    }
                    if (fdrBinary_01 && fdr >= 0.01) {
                        fdrBinary_01 = false;
                    }
                    if (fdrBinary_05 && fdr >= 0.05) {
                        fdrBinary_05 = false;
                    }
                    if (fdrBinary_20 && fdr >= 0.20) {
                        fdrBinary_20 = false;
                    }                    

                    // Format top genes from gene set overlapping with associated loci into string
                    String topGenes = "";
                    for (int j = 0; j < Math.min(maxTopGenesPerGeneSet,limitMaxTopGenesPerGeneSet); j++) {
                        if (j < geneSetGenesInLocus.size()) {
                            topGenes += "\t" + geneSetGenesInLocus.get(j);
                        } else { // If there are less than maxTopGenesPerGeneSet genes in the associated loci
                            topGenes += "\t-";
                        }
                    }
                    
                    if (fdrBinary_01) {
                        out.write(output + "\t<0.01" + topGenes + "\n");
                    } else if (fdrBinary_05) {
                        out.write(output + "\t<0.05" + topGenes + "\n");
                    } else if (fdrBinary_20) {
                        out.write(output + "\t<0.20" + topGenes + "\n");
                    } else {
                        out.write(output + "\t>=0.20" + topGenes + "\n");
                    }                   
                }
                out.flush();
                out.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }
        }

        HashMap hashGenesToInclude = hashGenesInRealLoci;
        if (conductNetworkAnalysis) {

            //We have now conducted a real analysis and a number of repetitions, now calculate the FDR:
            System.out.println("\n\n\n");
            StringDoubleObjectSorter stringDoubleObjectSorter = new StringDoubleObjectSorter();
            stringDoubleObjectSorter.sort(vecResultsRealAnalysis);
            double[] pValues = new double[nrUniqueGenesInLoci]; 
            for (int l = 0; l < nrUniqueGenesInLoci; l++) {
                StringDoubleObject stringDoubleObject = (StringDoubleObject) vecResultsRealAnalysis.get(l);
                pValues[l] = stringDoubleObject.doubleValue;
            }
            for (int rep = 0; rep < nrReps; rep++) {
                Arrays.sort(repPValues[rep]);
            }
            hashGenesToInclude = new HashMap();
            try {
                java.io.BufferedWriter out = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputGeneFileName)));
                out.write("Locus\tNr of genes in locus\tChromosome and position\tGWAS P value\tEnsembl gene ID\tGene symbol\tNominal P value\tGene closest to lead SNP\tTop cis eQTL SNP (Westra et al. Nature Genetics 2014)\tFalse discovery rate\n");
                double maxNrTPs = 0;
                boolean fdrBinary_01 = true;
                boolean fdrBinary_05 = true;
                boolean fdrBinary_20 = true;
                for (int l = 0; l < nrUniqueGenesInLoci; l++) {
                    double fdr = 0;
                    for (int rep = 0; rep < nrReps; rep++) {
                        for (int m = 0; m < nrUniqueGenesInLoci; m++) {
                            if (repPValues[rep][m] <= pValues[l]) {
                                fdr++;
                            }
                        }
                    }
                    fdr /= (double) nrReps;
                    fdr /= (double) (l + 1);
                    StringDoubleObject stringDoubleObject = (StringDoubleObject) vecResultsRealAnalysis.get(l);
                    double nrTPs = (1 - fdr) * (double) (l + 1);
                    if (nrTPs > maxNrTPs) {
                        maxNrTPs = nrTPs;
                    }
                    String output = stringDoubleObject.stringValue;
                    if (fdrBinary_01 && fdr >= 0.01) {
                        fdrBinary_01 = false;
                    }
                    if (fdrBinary_05 && fdr >= 0.05) {
                        fdrBinary_05 = false;
                    }
                    if (fdrBinary_20 && fdr >= 0.20) {
                        fdrBinary_20 = false;
                    } 
                    if (fdrBinary_01) {
                        out.write(output + "\t<=0.01\n");
                    } else if (fdrBinary_05) {
                        out.write(output + "\t<0.05\n");
                    } else if (fdrBinary_20) {
                        out.write(output + "\t<0.20\n");
                    } else {
                        out.write(output + "\t>0.20\n");
                    }
                }
                out.flush();
                out.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }
            
            // Write out genes outside loci
            if (calculateGenePrioritizationPValueForGenesOutsideLoci) {

                int nrUniqueGenesOutsideLoci = vecResultsOutsideRealAnalysis.size();
                System.out.println("\n\n\n");
                StringDoubleObjectSorter stringDoubleObjectSorterOutside = new StringDoubleObjectSorter();
                stringDoubleObjectSorterOutside.sort(vecResultsOutsideRealAnalysis); 
                double[] pValuesOutside = new double[nrUniqueGenesOutsideLoci]; 
                for (int l = 0; l < nrUniqueGenesOutsideLoci; l++) {
                        StringDoubleObject stringDoubleObject = (StringDoubleObject) vecResultsOutsideRealAnalysis.get(l);
                        pValuesOutside[l] = stringDoubleObject.doubleValue;
                }
                for (int rep = 0; rep < nrReps; rep++) {
                        Arrays.sort(repPValuesOutside[rep]);
                }

                try {
                    boolean fdrBinary_01 = true;
                    boolean fdrBinary_05 = true;
                    boolean fdrBinary_20 = true;   
                    java.io.BufferedWriter out = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputGeneGenomeFileName)));
                    out.write("Ensembl gene ID\tGene symbol\tChromosome\tTranscript start (hg19)\tTranscript end (hg19)\tP value\tFalse discovery rate\n");
                    for (int l = 0; l < nrUniqueGenesOutsideLoci; l++) { 
                            double fdr = 0;
                            for (int rep = 0; rep < nrReps; rep++) {
                                    for (int m = 0; m < nrUniqueGenesOutsideLoci; m++) {
                                            if (repPValuesOutside[rep][m] <= pValuesOutside[l]) {
                                                    fdr++;
                                            }
                                    }
                            }
                            fdr /= (double) nrReps;
                            fdr /= (double) (l + 1);
                            if (fdrBinary_01 && fdr >= 0.01) {
                                fdrBinary_01 = false;
                            }
                            if (fdrBinary_05 && fdr >= 0.05) {
                                fdrBinary_05 = false;
                            }
                            if (fdrBinary_20 && fdr >= 0.20) {
                                fdrBinary_20 = false;
                            }                    
                            StringDoubleObject stringDoubleObject = (StringDoubleObject) vecResultsOutsideRealAnalysis.get(l);
                            String[] tmp = stringDoubleObject.stringValue.split("\t");
                            String geneID = (String) hashEnsemblToHGNC.get(tmp[1]); // Get gene symbol
                            if (geneID.equals("-")) {
                                geneID = tmp[1];
                            }
                            String fdrString = ">=0.20";
                            if (fdrBinary_01) {
                                fdrString = "\t<0.01";
                            } else if (fdrBinary_05) {
                                fdrString = "\t<0.05";
                            } else if (fdrBinary_20) {
                                fdrString = "\t<0.20";
                            } 
                            out.write(tmp[1] + "\t" + geneID + "\t" + tmp[2] + "\t" + tmp[3] + "\t" + tmp[4] + "\t" + tmp[5] + fdrString + "\n");
                    }
                    out.flush();
                    out.close();
                } catch (Exception e) {
                    System.out.println("Error:\t" + e.getMessage());
                    e.printStackTrace();
                }
            }

        }
        System.exit(0);
    }

    public double[] getLinearRegressionCoefficients(double[] xVal, double[] yVal) {
        double n = (double) xVal.length;
        double sumX = 0;
        double sumXX = 0;
        double sumY = 0;
        double sumYY = 0;
        double sumXY = 0;
        for (int x = 0; x < xVal.length; x++) {
            sumX += xVal[x];
            sumXX += xVal[x] * xVal[x];
            sumY += yVal[x];
            sumYY += yVal[x] * yVal[x];
            sumXY += xVal[x] * yVal[x];
        }
        double sXX = sumXX - sumX * sumX / n;
        double sXY = sumXY - sumX * sumY / n;
        double a = sXY / sXX;
        double b = (sumY - a * sumX) / n;
        double[] regressionCoefficients = new double[2];
        regressionCoefficients[0] = a;
        regressionCoefficients[1] = b;
        return regressionCoefficients;
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
