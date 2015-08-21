package depict;

import depict.matrix.ExpressionDataset;
import java.io.File;
import java.io.InputStream;
import java.io.FileInputStream;
import java.util.zip.GZIPInputStream;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 *
 * @author DEPICTdevelopers
 */
public class Testing {
    
    public Testing() throws InterruptedException {
        //String filenameDatabaseToUse = "/Users/tp/Documents/Work/xDEPICT/data/GPL570EnsemblGeneExpressionPerTissue_DEPICT20130820_z.txt";
        //Load the predicted Z-Scores for a certain database (e.g. GO_BP, MGI, Reactome, KEGG etc.):
        //ExpressionDataset dataset = new ExpressionDataset(filenameDatabaseToUse, "\t", null, null);
        //dataset.adjustData(0);
        
        try {
            String filename = "/Users/tp/Downloads/test.tab.gz";
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(filename)));
            System.err.println(in.readLine());

            String encoding = "US-ASCII";
            InputStream fileStream = new FileInputStream(filename);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, encoding);
            BufferedReader in2 = new BufferedReader(decoder);
            System.err.println(in2.readLine());
            
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
            System.exit(-1);
        }
/*        
        cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.DRand();
        
        // T-test
        double mean1 = 3.3;
        double mean2 = 2.2;
        double var1 = 0.3;
        double var2 = 0.1;
        int n1 = 30;
        int n2 = 40;
        double t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);
        
        // Correlation
        double correlation = 0.8;
        int vals1_length = 10000;
        cern.jet.random.StudentT tDistColt = new cern.jet.random.StudentT(vals1_length/2 - 2, randomEngine);
        //double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (vals1_length/2 - 2)));        

        // t>0 scenario
        if (t<0) {
            System.out.println("t<0");
            double p = tDistColt.cdf(t);        
            System.out.println("t = " + t);
            System.out.println("p = " + p);
            System.out.println("z = " + cern.jet.stat.Probability.normalInverse(p));            
        } else {
            System.out.println("t>0");
            double p = tDistColt.cdf(-t);        
            System.out.println("t = " + t);
            System.out.println("p = " + p);
            System.out.println("z = " + -cern.jet.stat.Probability.normalInverse(p));            
        }
  
*/
        
    }
}
