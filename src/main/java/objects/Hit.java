package objects;

import java.util.ArrayList;
import java.util.List;
import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

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

    private True tru;
    
    public Hit(int index, int sector, int layer, int component, DetectorType type) {
        this.index = index;
        this.sector = sector;
        this.layer = layer;
        this.strip = component;
        this.type = type;
    }

    public int getIndex() {
        return index;
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

    public True getTrue() {
        return tru;
    }

    public void setTrue(True t) {
        this.tru = t;
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
    
    public static List<Hit> readHits(DataEvent event, DetectorType type) {
        DataBank adcBank = event.getBank(type.getName()+"::adc");
        DataBank tdcBank = event.getBank(type.getName()+"::tdc");
        if(adcBank!=null && adcBank.rows()>0) {
            List<Hit> hits = new ArrayList<>();
            for(int row=0; row<adcBank.rows(); row++) {
                Hit hit = new Hit(row,
                                  adcBank.getByte("sector", row),
                                  adcBank.getByte("layer", row),
                                  adcBank.getInt("component", row),
                                  type);
                hit.setADC(adcBank.getInt("ADC", row));
                hit.setTime(adcBank.getFloat("time", row));
                hits.add(hit);
            }
            return hits;
        }
        else if(tdcBank!=null && tdcBank.rows()>0) {
            List<Hit> hits = new ArrayList<>();
            for(int row=0; row<tdcBank.rows(); row++) {
                Hit hit = new Hit(row,
                                  tdcBank.getByte("sector", row),
                                  tdcBank.getByte("layer", row),
                                  tdcBank.getInt("component", row),
                                  type);
                hit.setTDC(tdcBank.getInt("TDC", row));
                hits.add(hit);
            }
            return hits;                                        
        }
        return null;
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
