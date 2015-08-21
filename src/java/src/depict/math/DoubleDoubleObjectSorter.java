package depict.math;

/**
 *
 * @author Lude Franke and DEPICTdevelopers
 */
public class DoubleDoubleObjectSorter extends VectorSorter {
    
    /** Creates a new instance of GeneLocationObjectSorter */
    public DoubleDoubleObjectSorter() {
        super();
    }

    /** Override object comparer
     * @param a the first GeneLocationObject to be compared
     * @param b the second GeneLocationObject to be compared
     * @return true if the first GeneLocationObject.getChrStart() is lower than the second one
     */
    protected boolean lt (Object a, Object b) {
        return (((DoubleDoubleObject)a).doubleValueToSortOn < ((DoubleDoubleObject)b).doubleValueToSortOn);
    }
    
}
