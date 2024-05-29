package objects;

import analysis.Constants;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class Event {
    private final boolean debug = false;
    
    private int run;
    private int event;
    private double startTime;
    private double[] raster = new double[2];
    private final Map<DetectorType, List<Hit>> hits   = new HashMap<>();
    private final Map<DetectorType, List<True>> trues = new HashMap<>();
    private DataEvent hipoEvent;

    public Event(DataEvent event, List<DetectorType> types) {
        this.hipoEvent = event;
        this.readEvent(event, types);
    }
    


    private DataBank getBank(DataEvent de, String bankName) {
        DataBank bank = null;
        if (de.hasBank(bankName)) {
            bank = de.getBank(bankName);
        }
        return bank;
    }

    private DataBank getBank(DataEvent de, String bankName1, String bankName2) {
        DataBank bank = null;
        if (de.hasBank(bankName1)) {
            bank = de.getBank(bankName1);
        }
        else if (de.hasBank(bankName2)) {
            bank = de.getBank(bankName2);
        }
        return bank;
    }

    private void readHeader(DataEvent event) {
        DataBank head = this.getBank(event, "RUN::config");
        if(head!=null) {
            this.run   = head.getInt("run", 0);
            this.event = head.getInt("event", 0);
        }
    }
    
    private void readStartTime(DataEvent event) {
        DataBank recEvent = this.getBank(event, "REC::Event");
        if(recEvent!=null) {
            startTime = recEvent.getFloat("startTime",0);
        }
    }
    
    
    private void readRaster(DataEvent event) {
        DataBank rasterPos = this.getBank(event, "RASTER::position");
        if(rasterPos!=null) {
            raster[0] = rasterPos.getFloat("x",0);
            raster[1] = rasterPos.getFloat("y",0);
        }
    }
    
    private void readHits(DataEvent event, DetectorType type) {
        DataBank bank = this.getBank(event, "DC::tdc");
        if(bank==null) return;
        for(int i=0; i<bank.rows(); i++) {
            Hit hit = Hit.readTDC(bank, i, DetectorType.DC);
            if(hits.get(type)==null) 
                hits.put(type, new ArrayList<>());
            hits.get(type).add(hit);
        }
    }

    private void readTrues(DataEvent event, DetectorType type) {
        DataBank mc  = this.getBank(event, "MC::True");
        if(mc==null) return;
        for(int i=0; i<mc.rows(); i++) {
            if(mc.getByte("detector", i)==type.getDetectorId()) {
                True t = True.readTruth(mc, i);
                if(!trues.containsKey(type)) 
                    trues.put(type, new ArrayList<>());
                trues.get(type).add(t);
            }
//        Collections.sort(trues);
        }
    }
    

    private void readEvent(DataEvent de, List<DetectorType> types) {
        this.readHeader(de);
        this.readStartTime(de);
        this.readRaster(de);
        for(DetectorType t : types) {
            this.readHits(de, t);
            this.readTrues(de, t);
        }
    }


    public List<Hit> getHits(DetectorType type) {
        return hits.get(type);
    }

    public List<True> getTrues(DetectorType type) {
        return trues.get(type);
    }

    public int getRun() {
        return run;
    }

    public int getEvent() {
        return event;
    }

    
    public double getStartTime() {
        return startTime;
    }

    public double[] getRaster() {
        return raster;
    }

    public DataEvent getHipoEvent() {
        return hipoEvent;
    }
    
    
}
