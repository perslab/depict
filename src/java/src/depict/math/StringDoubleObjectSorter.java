package depict.math;

/**
 *
 * @author Lude Franke and DEPICTdevelopers
 */
public class StringDoubleObjectSorter extends VectorSorter {
    
    /** Creates a new instance of GeneLocationObjectSorter */
    public StringDoubleObjectSorter() {
        super();
    }

    /** Override object comparer
     * @param a the first GeneLocationObject to be compared
     * @param b the second GeneLocationObject to be compared
     * @return true if the first GeneLocationObject.getChrStart() is lower than the second one
     */
    protected boolean lt (Object a, Object b) {
        return (((StringDoubleObject)a).doubleValue < ((StringDoubleObject)b).doubleValue);
    }
    
}
