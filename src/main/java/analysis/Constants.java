package analysis;

import org.jlab.clas.pdg.PhysicsConstants;

/**
 *
 * @author devita
 */
public class Constants {
    
    public static final int BMTSECTORS = 3;
    public static final int BMTREGIONS = 3;
    public static final int SVTLAYERS = 6;
    public static final int[] SVTSECTORS = {10, 10, 14, 14, 18, 18};
    public static final int SVTSTRIPS = 200;
    public static final int[] BMTCLAYERS = {1,4,6};
    public static final int[] BMTZLAYERS = {2,3,5};
    public static final int[] BMTCSTRIPS = {896, 1024, 1152};
    public static final int[] BMTZSTRIPS = {640, 640, 768};
    public static final int[] BMTZPITCH  = {487, 536, 529};
    public static final int[] BMTMEANPHI = {-150, 90, -30};
    public static final double[] BMTRADIUS = {14.7746, 16.2646, 17.7646, 19.2646, 20.7646, 22.2646}; // radius at half drift
    
    public static final double B = 5; //field in Tesla
    public static final double ALPHA = 1E6/PhysicsConstants.speedOfLight()/B;

    public static double TWINDOW = 250; // ns

    public static double getTimeWindow() {
        return TWINDOW;
    }

    public static void setTimeWindow(double window) {
        Constants.TWINDOW = window;
    }
    
}
