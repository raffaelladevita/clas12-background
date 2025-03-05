package modules;

import analysis.Constants;
import java.util.List;
import objects.Hit;
import analysis.Module;
import objects.Event;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;

/**
 *
 * @author devita
 */
public class CNDmodule extends Module {
    
    private static final int NPADDLE  = 48;
    private static final int NLAYER   = 3;
    
    private static final double THRESHOLD = 0.1; // MeV FIXME
    private static final double WEIGHT = 1.5*1.5*20*8.28/1000; //kg  FIXME
    
    public CNDmodule() {
        super(DetectorType.CND);
    }
    
    public DataGroup occupancies() {
        DataGroup dg = new DataGroup(3,2);
        H2F hi_occ_2D     = histo2D("hi_occ_2D", "Occupancy (%)", "X", "Y", NPADDLE+2 , 0, NPADDLE+1, NLAYER+2 , 0, NLAYER+1);           
        H2F hi_rate_2D    = histo2D("hi_rate_2D", "Rate (kHz)", "X", "Y", NPADDLE+2 , 0, NPADDLE+1, NLAYER+2 , 0, NLAYER+1);           
        H2F hi_dose_2D    = histo2D("hi_dose_2D", "Dose (rad/h)", "X", "Y", NPADDLE+2 , 0, NPADDLE+1, NLAYER+2 , 0, NLAYER+1);           
        H2F hi_edep_2D    = histo2D("hi_edep_2D", "Energy deposition rate (MeV/us)", "X", "Y", NPADDLE+2, 0, NPADDLE+1, NLAYER+2 , 0, NLAYER+1);           
        H1F hi_occ_1D     = histo1D("hi_occ_1D",  " ", "ID", "Occupancy (%)", NPADDLE*NLAYER, 0, NPADDLE*NLAYER, 1);           
        H1F hi_edep_1D    = histo1D("hi_edep_1D",  " ", "Energy (MeV)", "Rate kHz)", 100, 0, 500, 4); 
        H1F hi_time_1D    = histo1D("hi_time_1D",  " ", "Time (ns)", "Rate kHz)", 100, 0, Constants.getTimeWindow()*1.2, 4);           
        dg.addDataSet(hi_occ_2D,    0);
        dg.addDataSet(hi_rate_2D,   1);
        dg.addDataSet(hi_edep_2D, 2);
        dg.addDataSet(hi_edep_1D,   3);
        dg.addDataSet(hi_time_1D,   4);
        dg.addDataSet(hi_dose_2D,   5);
        return dg;
    }
  
    @Override
    public void createHistos() {
        this.getHistos().put("Occupancy", this.occupancies());
    }
    
    @Override
    public void fillHistos(Event event) {
        if (event.getHits(DetectorType.CND) != null) {
            this.fillOccupancies(this.getHistos().get("Occupancy"), event.getHits(DetectorType.CND));
        }
    }
    
    public void fillOccupancies(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            if(hit.getOrder()>1) continue;  // use ADCs only
            int idy = hit.getLayer();
            int idx = hit.getSector()*2+hit.getOrder();
            double edep = hit.getTrue().getEdep();
            if(edep>THRESHOLD) {
                group.getH2F("hi_occ_2D").fill(idx, idy);
                group.getH2F("hi_rate_2D").fill(idx, idy);
//                group.getH1F("hi_occ_1D").fill(hit.getComponent());
                group.getH1F("hi_time_1D").fill(hit.getTrue().getTime());
            }
            group.getH1F("hi_edep_1D").fill(edep);
            group.getH2F("hi_edep_2D").fill(idx, idy, edep);
            group.getH2F("hi_dose_2D").fill(idx, idy, edep/WEIGHT);
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.normalizeToEventsX100(this.getHistos().get("Occupancy").getH2F("hi_occ_2D"));
//        this.normalizeToEventsX100(this.getHistos().get("Occupancy").getH1F("hi_occ_1D"));
        this.normalizeToTime(this.getHistos().get("Occupancy").getH2F("hi_rate_2D"));
        this.normalizeToTime(this.getHistos().get("Occupancy").getH2F("hi_edep_2D"), 1E6);
        this.normalizeToTime(this.getHistos().get("Occupancy").getH1F("hi_edep_1D"));
        this.normalizeToTime(this.getHistos().get("Occupancy").getH1F("hi_time_1D"));
        this.toDose(this.getHistos().get("Occupancy").getH2F("hi_dose_2D"));
    }
    
    @Override
    public void setPlottingOptions(String key) {
        if(key.equals("Occupancy")) {
            this.getCanvas(key).getCanvasPads().get(1).getAxisZ().setLog(true);
            this.getCanvas(key).getCanvasPads().get(2).getAxisZ().setLog(true);
            this.getCanvas(key).getCanvasPads().get(3).getAxisY().setLog(true);
        }
    }
   
}
