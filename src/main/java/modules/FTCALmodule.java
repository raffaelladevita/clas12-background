package modules;

import analysis.Constants;
import java.util.List;
import objects.Hit;
import analysis.Module;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import objects.Event;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
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
    private static final int NCRYSTALX = 44;
    private static final int NCRYSTALY = 44;
    
    private static final double THRESHOLD = 20; // MeV
    private static final double WEIGHT = 1.5*1.5*20*8.28/1000; //kg
    
    private BufferedWriter lundFile = null;
    private PhysicsEvent physicsEvent = new PhysicsEvent();

    public FTCALmodule() {
        super(DetectorType.FTCAL);
        try {
            this.lundFile = new BufferedWriter(new FileWriter("pions.lund"));
            }
            catch(IOException e) {
                System.out.println(e.getMessage());
        }    
    }

    public DataGroup occupancies() {
        DataGroup dg = new DataGroup(3,2);
        H2F hi_occ_2D     = histo2D("hi_occ_2D", "Occupancy (%)", "X", "Y", NCRYSTALX , 0.5, NCRYSTALX+0.5, NCRYSTALY, 0.5, NCRYSTALY+0.5);           
        H2F hi_rate_2D    = histo2D("hi_rate_2D", "Rate (kHz)", "X", "Y", NCRYSTALX , 0.5, NCRYSTALX+0.5, NCRYSTALY, 0.5, NCRYSTALY+0.5);           
        H2F hi_dose_2D    = histo2D("hi_dose_2D", "Dose (rad/h)", "X", "Y", NCRYSTALX , 0.5, NCRYSTALX+0.5, NCRYSTALY, 0.5, NCRYSTALY+0.5);           
        H2F hi_edep_2D    = histo2D("hi_edep_2D", "Energy deposition rate (MeV/us)", "X", "Y", NCRYSTALX , 0.5, NCRYSTALX+0.5, NCRYSTALY, 0.5, NCRYSTALY+0.5);          
        H1F hi_occ_1D     = histo1D("hi_occ_1D",  " ", "ID", "Occupancy (%)", NCRYSTALX*NCRYSTALY, 0, NCRYSTALX*NCRYSTALY, 1);           
        H1F hi_edep_1D    = histo1D("hi_edep_1D",  " ", "Energy (MeV)", "Rate kHz)", 100, 0, 1000, 4); 
        H1F hi_time_1D    = histo1D("hi_time_1D",  " ", "Time (ns)", "Rate kHz)", 100, 0, Constants.getTimeWindow()*1.2, 4);           
        dg.addDataSet(hi_occ_2D,    0);
        dg.addDataSet(hi_rate_2D,   1);
        dg.addDataSet(hi_edep_2D, 2);
        dg.addDataSet(hi_edep_1D,   3);
        dg.addDataSet(hi_time_1D,   4);
        dg.addDataSet(hi_dose_2D,   5);
        return dg;
    }
  
    public DataGroup pions() {
        DataGroup dg = new DataGroup(2,2);
        H1F hi_mom_1D   = histo1D("hi_mom_1D",  " ", "p (GeV)", "Counts", 100, 0, 5, 4); 
        H1F hi_edep_1D  = histo1D("hi_edep_1D",  " ", "Energy (MeV)", "Counts", 100, 0, 5000, 4); 
        H1F hi_ecut_1D  = histo1D("hi_ecut_1D",  " ", "Energy (MeV)", "Counts", 100, 0, 5000, 2); 
        H2F hi_mom_2D   = histo2D("hi_mom_2D",  " ", "p (GeV)", "#theta (deg)", 100, 0, 5, 100, 5, 40); 
        H2F hi_cut_2D   = histo2D("hi_ecut_2D",  " ", "p (GeV)", "#theta (deg)", 100, 0, 5, 100, 5, 40); 
        dg.addDataSet(hi_mom_2D,   0);
        dg.addDataSet(hi_cut_2D,   1);
        dg.addDataSet(hi_mom_1D,   2);
        dg.addDataSet(hi_edep_1D,  3);
        dg.addDataSet(hi_ecut_1D,  3);
        return dg;
    }
    
    @Override
    public void createHistos() {
        this.getHistos().put("Occupancy", this.occupancies());
        this.getHistos().put("Pions", this.pions());
    }
    
    @Override
    public void fillHistos(Event event) {
        if (event.getHits(DetectorType.FTCAL) != null) {
            this.fillOccupancies(this.getHistos().get("Occupancy"), event.getHits(DetectorType.FTCAL));
            this.fillPions(this.getHistos().get("Pions"), event);
        }
    }
    
    public void fillOccupancies(DataGroup group, List<Hit> hits) {
//	int component = (IDY-1)*44+IDX-1;
        for (Hit hit : hits) {
            int idy = hit.getComponent()/NCRYSTALY+1;
            int idx = hit.getComponent()%NCRYSTALY+1;
            double edep = hit.getTrue()!=null ? hit.getTrue().getEdep() : THRESHOLD*1.01;
            double time = hit.getTrue()!=null ? hit.getTrue().getTime(): 0;
            if(edep==0)
                continue;
            else if(edep>THRESHOLD) {
                group.getH2F("hi_occ_2D").fill(idx, idy, 1/2.5);
                group.getH2F("hi_rate_2D").fill(idx, idy);
//                group.getH1F("hi_occ_1D").fill(hit.getComponent());
                group.getH1F("hi_time_1D").fill(time);
            }
            group.getH1F("hi_edep_1D").fill(edep);
            group.getH2F("hi_edep_2D").fill(idx, idy, edep);
            group.getH2F("hi_dose_2D").fill(idx, idy, edep/WEIGHT);
        }
    }
    
    public void fillPions(DataGroup group, Event event) {
        double etot = 0;
        Particle mc = event.getGeneratedParticle();
        for (Hit hit : event.getHits(DetectorType.FTCAL)) {
            double edep = hit.getTrue()!=null ? hit.getTrue().getEdep() : 0;
                etot +=edep;
        }
        if(true) {
            group.getH2F("hi_mom_2D").fill(mc.p(),Math.toDegrees(mc.theta()));
            group.getH1F("hi_mom_1D").fill(mc.p());
            group.getH1F("hi_edep_1D").fill(etot);
//            etot = 0;
            if(event.getHipoEvent().hasBank("FTCAL::clusters")) {
                int size = event.getHipoEvent().getBank("FTCAL::clusters").getShort("size", 0); 
                if(size>3)
                    etot = event.getHipoEvent().getBank("FTCAL::clusters").getFloat("recEnergy", 0)*1000;
            }
            if(etot>500) {
                group.getH2F("hi_ecut_2D").fill(mc.p(),Math.toDegrees(mc.theta()));
                group.getH1F("hi_ecut_1D").fill(etot);
            }
            mc.setTheta(Math.toRadians(3.5));
            this.physicsEvent.clear();
            this.physicsEvent.addParticle(mc);
            try {
                lundFile.write(this.physicsEvent.toLundString());
            } catch (IOException ex) {
                Logger.getLogger(FTCALmodule.class.getName()).log(Level.SEVERE, null, ex);
            }
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
        if(lundFile!=null)
            try {
                lundFile.close();
        } catch (IOException ex) {
            Logger.getLogger(FTCALmodule.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    @Override
    public void setPlottingOptions(String key) {
        if(key.equals("Occupancy")) {
//            this.getCanvas(key).getCanvasPads().get(1).getAxisZ().setLog(true);
            this.getCanvas(key).getCanvasPads().get(2).getAxisZ().setLog(true);
            this.getCanvas(key).getCanvasPads().get(3).getAxisY().setLog(true);
        }
    }
   
}
