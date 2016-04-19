package uninovo.util;

import java.util.Comparator;

/**
 * This class represents a pair of two objects.
 * 
 * @param <A> the first object
 * @param <B> the second object
 * 
 * @author sangtaekim
 */
public class Pair<A, B> {
    
    public static class PairComparator<A extends Comparable<? super A> ,B extends Comparable<? super B>> implements Comparator<Pair<A,B>> {
    	boolean useSecondForComprison;
    	public PairComparator()
    	{
    		this(false);
    	}
    	public PairComparator(boolean useSecondForComprison)
    	{
    		this.useSecondForComprison = useSecondForComprison;
    	}
        /**
         * Determines the order of Pair objects. If useSecondForComparison is set, use B for comparison, otherwise A is used.
         * @param o1 the first element.
         * @param o2 the second element.
         * @return 1 if p1 > p2, -1 if p2 > p1 and 0 otherwise.
         */
        @Override
		public int compare(Pair<A,B> p1, Pair<A,B> p2) {
        	if(!useSecondForComprison)
        		return p1.getFirst().compareTo(p2.getFirst());
        	else
        		return p1.getSecond().compareTo(p2.getSecond());
        }
    }
    
    /** The first. */
    private A first;

    /** The second. */
    private B second;

    /**
     * Instantiates a new pair.
     * 
     * @param first the first
     * @param second the second
     */
    public Pair(A first, B second) {
        super();
        this.first = first;
        this.second = second;
    }

    /* (non-Javadoc)
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
	public boolean equals(Object other) {
        if (other instanceof Pair<?,?>) {
                Pair<?,?> otherPair = (Pair<?,?>) other;
                return 
                ((  this.first == otherPair.first ||
                        ( this.first != null && otherPair.first != null &&
                          this.first.equals(otherPair.first))) &&
                 (      this.second == otherPair.second ||
                        ( this.second != null && otherPair.second != null &&
                          this.second.equals(otherPair.second))) );
        }

        return false;
    }

    /**
     * Gets the first.
     * 
     * @return the first
     */
    public A getFirst() {
        return first;
    }

    /**
     * Gets the second.
     * 
     * @return the second
     */
    public B getSecond() {
        return second;
    }

    /* (non-Javadoc)
     * @see java.lang.Object#hashCode()
     */
    @Override
	public int hashCode() {
        int hashFirst = first != null ? first.hashCode() : 0;
        int hashSecond = second != null ? second.hashCode() : 0;

        return (hashFirst + hashSecond) * hashSecond + hashFirst;
    }

    /**
     * Sets the first.
     * 
     * @param first the new first
     */
    public void setFirst(A first) {
        this.first = first;
    }

    /**
     * Sets the second.
     * 
     * @param second the new second
     */
    public void setSecond(B second) {
        this.second = second;
    }
    
    /* (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    @Override
	public String toString()
    { 
           return "(" + first + ", " + second + ")"; 
    }
}
