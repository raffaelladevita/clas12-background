package objects;

import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class Hit {
 
    private final int index;
    private final int sector;
    private final int layer;
    private final int strip;
    private final DetectorType type;
    
    private int adc;
    private int tdc;
    private double time;
    
    public Hit(int index, int sector, int layer, int component, DetectorType type) {
        this.index = index;
        this.sector = sector;
        this.layer = layer;
        this.strip = component;
        this.type = type;
    }

    public int getSector() {
        return sector;
    }

    public int getLayer() {
        return layer;
    }

    public int getComponent() {
        return strip;
    }

    public DetectorType getType() {
        return type;
    }

    public int getADC() {
        return adc;
    }

    public void setADC(int adc) {
        this.adc = adc;
    }

    public int getTDC() {
        return tdc;
    }

    public void setTDC(int tdc) {
        this.tdc = tdc;
    }

    public double getTime() {
        return time;
    }

    public void setTime(double time) {
        this.time = time;
    }

    public String getName() {
        if(this.type==DetectorType.BST)
            return "SVT";
        else {
            if(layer==1 || layer==4 || layer==6)
                return "BMTC";
            else
                return "BMTZ";
        }
    }
    
    public static Hit readADC(DataBank bank, int row, DetectorType type) {
        Hit hit = new Hit(row,
                          bank.getByte("sector", row),
                          bank.getByte("layer", row),
                          bank.getInt("component", row),
                          type);
        hit.setADC(bank.getInt("ADC", row));
        hit.setTime(bank.getFloat("time", row));
        return hit;
    }
    
    public static Hit readTDC(DataBank bank, int row, DetectorType type) {
        Hit hit = new Hit(row,
                          bank.getByte("sector", row),
                          bank.getByte("layer", row),
                          bank.getInt("component", row),
                          type);
        hit.setTDC(bank.getInt("TDC", row));
        return hit;
    }
    
}
