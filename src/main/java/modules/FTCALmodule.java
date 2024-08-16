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
public class FTCALmodule extends Module {
    
    private static final int NCRYSTAL = 332;
    private static final int NCRYSTALX = 22;
    private static final int NCRYSTALY = 22;
    
    private static final double THRESHOLD = 10; // MeV
    private static final double WEIGHT = 1.5*1.5*20*8.28/1000; //kg
    
    public FTCALmodule(int residualScale) {
        super(DetectorType.FTCAL);
    }
    
    public DataGroup dcHitVertices(int col) {
        String[] names = new String[]{"SVT", "BMTC", "BMTZ"};
        double[] EMAX = { 1000, 5000, 5000};
        double[] TMAX = {  1200, 1200, 1200};
        double[] RMAX = { 1500, 3000, 3000};
        DataGroup dg = new DataGroup(3,3);
        for(int i=0; i<names.length; i++) {
            String name = names[i];
            H1F hi_energy  = histo1D("hi_energy_" + name, name, name + " Hit Energy",   "Counts", 100,   0, EMAX[i], col);
            H1F hi_time    = histo1D("hi_time_" + name,name,   name + " Hit Time",     "Counts", 10+(int) TMAX[i], -10, TMAX[i], col);
            H1F hi_resi    = histo1D("hi_resi_" + name, name,   name + " Hit Residual (um)", "Counts", 100, -RMAX[i], RMAX[i], col);

            dg.addDataSet(hi_energy, 0 + i*3);
            dg.addDataSet(hi_time,   1 + i*3);
            dg.addDataSet(hi_resi,   2 + i*3);
        }
        return dg;
    }

    

    public DataGroup occupancies() {
        DataGroup dg = new DataGroup(3,2);
        H2F hi_occ_2D     = histo2D("hi_occ_2D", "Occupancy (%)", "X", "Y", NCRYSTALX+3, 0, NCRYSTALX+2, NCRYSTALY+3, 0, NCRYSTALY+2);           
        H2F hi_rate_2D    = histo2D("hi_rate_2D", "Rate (MHz)", "X", "Y", NCRYSTALX+3, 0, NCRYSTALX+2, NCRYSTALY+3, 0, NCRYSTALY+2);           
        H2F hi_dose_2D    = histo2D("hi_dose_2D", "Dose (rad/h)", "X", "Y", NCRYSTALX+3, 0, NCRYSTALX+2, NCRYSTALY+3, 0, NCRYSTALY+2);           
        H2F hi_energy_2D  = histo2D("hi_energy_2D", "Energy deposition rate (MeV/us)", "X", "Y", NCRYSTALX+3, 0, NCRYSTALX+2, NCRYSTALY+3, 0, NCRYSTALY+2);           
        H1F hi_occ_1D     = histo1D("hi_occ_1D",  " ", "ID", "Occupancy (%)", NCRYSTALX*NCRYSTALY, 0, NCRYSTALX*NCRYSTALY, 1);           
        H1F hi_edep_1D    = histo1D("hi_edep_1D",  " ", "Energy (MeV)", "Rate MHz)", 100, 0, 500, 4);           
        H1F hi_time_1D    = histo1D("hi_time_1D",  " ", "Time (ns)", "Rate MHz)", 100, 0, Constants.getTimeWindow()*1.2, 4);           
        dg.addDataSet(hi_occ_2D,    0);
        dg.addDataSet(hi_rate_2D,   1);
        dg.addDataSet(hi_energy_2D, 2);
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
        if (event.getHits(DetectorType.FTCAL) != null) {
            this.fillOccupancies(this.getHistos().get("Occupancy"), event.getHits(DetectorType.FTCAL));
        }
    }
    
    public void fillOccupancies(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            int idy = hit.getComponent()/NCRYSTALY+1;
            int idx = hit.getComponent()%NCRYSTALY+1;
            double edep = hit.getTrue().getEdep();
            if(edep>THRESHOLD) {
                group.getH2F("hi_occ_2D").fill(idx, idy);
                group.getH2F("hi_rate_2D").fill(idx, idy);
//                group.getH1F("hi_occ_1D").fill(hit.getComponent());
                group.getH1F("hi_time_1D").fill(hit.getTrue().getTime());
            }
            group.getH1F("hi_edep_1D").fill(edep);
            group.getH2F("hi_energy_2D").fill(idx, idy, edep);
            group.getH2F("hi_dose_2D").fill(idx, idy, edep/WEIGHT);
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.normalizeToEventsX100(this.getHistos().get("Occupancy").getH2F("hi_occ_2D"));
        this.normalizeToEventsX100(this.getHistos().get("Occupancy").getH1F("hi_occ_1D"));
        this.normalizeToTime(this.getHistos().get("Occupancy").getH2F("hi_rate_2D"));
        this.normalizeToTime(this.getHistos().get("Occupancy").getH2F("hi_energy_2D"));
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
