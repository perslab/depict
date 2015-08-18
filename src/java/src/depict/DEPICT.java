
package depict;

/**
 *
 * @author DEPICTdevelopers
 */
public class DEPICT {

    public static void main(String[] args) {
     	try {
            new PerformPathwayAnalysisSusceptibilityLoci(args);
            //new PredictGenesetMemberBasedOnTCs(args);
            //new Testing();
        } catch (java.lang.InterruptedException e) {
		System.out.println(e.getMessage());
		e.printStackTrace();
	}
    }
}