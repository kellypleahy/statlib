/*
 * ChiSquaredDist.java
 *
 * Created on September 30, 2002, 11:06 AM
 */

package statlib;

/**
 *
 * @author  KLeahy
 */
public class ChiSquaredDist extends GammaDist {
  
  /** Creates a new instance of ChiSquaredDist */
  public ChiSquaredDist(int r) {
    super(r / 2, 2);
  }
  
}
