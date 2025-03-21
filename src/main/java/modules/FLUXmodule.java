package modules;

import analysis.Constants;
import java.util.List;
import objects.Hit;
import analysis.Module;
import objects.Event;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;

/**
 *
 * @author devita
 */
public class FLUXmodule extends Module {
    
    private static final double DTHETA = 0.5;
    private static final double[] R = {50, 700}; // cm
    
    public FLUXmodule() {
        super(DetectorType.TARGET);
    }
    
    public DataGroup rates() {
        DataGroup dg = new DataGroup(2,2);
        
        for(int il=0; il<2; il++) {
            double range = 100-55*il;
            String title = il==0 ? "Central" : "Forward";
            for (int ip=0; ip<PNAMES.length; ip++) {
                H1F hi_all = histo1D("hi_all_1D_"+(il+1) + "_" + PNAMES[ip], title, "#theta (deg)", "Flux [Hz/d#Omega] ", (int) (range/DTHETA), 0, range, 0); 
                this.setHistoAttr(hi_all, ip<5 ? ip+1 : ip+3);
                H1F hi_bwd = histo1D("hi_bwd_1D_"+(il+1) + "_" + PNAMES[ip], title, "#theta (deg)", "Flux [Hz/d#Omega] ", (int) (range/DTHETA), 0, range, 0); 
                this.setHistoAttr(hi_bwd, ip<5 ? ip+1 : ip+3);
                dg.addDataSet(hi_all, 0 + il*2);
                dg.addDataSet(hi_bwd, 1 + il*2);
            }
        }
        return dg;
    }
  
    public DataGroup fluxes() {
        DataGroup dg = new DataGroup(3,2);
        
        for(int il=0; il<2; il++) {
            H2F hi_all = histo2D("hi_all_2D_"+(il+1), "Total Flux (Hz/cm2)", "x (cm)", "y (cm) ", (int) R[il]*2, -R[il], R[il], (int) R[il]*2, -R[il], R[il]); 
            H2F hi_crg = histo2D("hi_crg_2D_"+(il+1), "Charged Flux (Hz/cm2)", "x (cm)", "y (cm) ", (int) R[il]*2, -R[il], R[il], (int) R[il]*2, -R[il], R[il]); 
            H2F hi_neu = histo2D("hi_neu_2D_"+(il+1), "Neutron Flux (Hz/cm2)", "x (cm)", "y (cm) ", (int) R[il]*2, -R[il], R[il], (int) R[il]*2, -R[il], R[il]); 
            dg.addDataSet(hi_all, 0 + il*3);
            dg.addDataSet(hi_crg, 1 + il*3);
            dg.addDataSet(hi_neu, 2 + il*3);
        }
        return dg;
    }
  
    @Override
    public void createHistos() {
        this.getHistos().put("Rates", this.rates());
        this.getHistos().put("Fluxes", this.fluxes());
    }
    
    @Override
    public void fillHistos(Event event) {
        if (event.getHits(DetectorType.TARGET) != null) {
            this.fillRates(this.getHistos().get("Rates"), event.getHits(DetectorType.TARGET));
            this.fillFluxes(this.getHistos().get("Fluxes"), event.getHits(DetectorType.TARGET));
        }
    }
    
    public void fillRates(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            double theta = hit.getTrue().getPosition().toVector3D().theta();
            double domega = 2*Math.PI*Math.sin(theta)*Math.toRadians(DTHETA);
            int il = hit.getComponent()/10;
            boolean bwd = hit.getTrue().getPosition().toVector3D().dot(hit.getTrue().getMomentum())<0;
            group.getH1F("hi_all_1D_" + (il+1) + "_all").fill(Math.toDegrees(theta), 1/domega);
            if(bwd)
                group.getH1F("hi_bwd_1D_" + (il+1) + "_all").fill(Math.toDegrees(theta), 1/domega);                
            if(this.pidToName(Math.abs(Math.abs(hit.getTrue().getPid())))!=null) {
                group.getH1F("hi_all_1D_" + (il+1) + "_" + this.pidToName(Math.abs(hit.getTrue().getPid()))).fill(Math.toDegrees(theta), 1/domega);
                if(bwd)
                    group.getH1F("hi_bwd_1D_" + (il+1) + "_" + this.pidToName(Math.abs(hit.getTrue().getPid()))).fill(Math.toDegrees(theta), 1/domega);
            }
        }
    }
    
    
    public void fillFluxes(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            Vector3D proj = hit.getTrue().getPosition().toVector3D().asUnit();
            double theta = hit.getTrue().getPosition().toVector3D().theta();
            if(theta>Math.toRadians(70)) continue;
            int pid = Math.abs(hit.getTrue().getPid());
            int il = hit.getComponent()/10;
            group.getH2F("hi_all_2D_" + (il+1)).fill(proj.x()*R[il],proj.y()*R[il]);
            if(pid==11 || pid==2212 || pid==211) {
                group.getH2F("hi_crg_2D_" + (il+1)).fill(proj.x()*R[il],proj.y()*R[il]);
            }
            else if(pid==2112) {
                group.getH2F("hi_neu_2D_" + (il+1)).fill(proj.x()*R[il],proj.y()*R[il]);
            }
        }
    }
    
    @Override
    public void analyzeHistos() {
        for(int il=0; il<2; il++) {
            this.normalizeToTime(this.getHistos().get("Fluxes").getH2F("hi_all_2D_" + (il+1)), 1);
            this.normalizeToTime(this.getHistos().get("Fluxes").getH2F("hi_crg_2D_" + (il+1)), 1);
            this.normalizeToTime(this.getHistos().get("Fluxes").getH2F("hi_neu_2D_" + (il+1)), 1);
            for (String pn : PNAMES) {
                this.normalizeToTime(this.getHistos().get("Rates").getH1F("hi_all_1D_" + (il+1) + "_" + pn), 1);
                this.normalizeToTime(this.getHistos().get("Rates").getH1F("hi_bwd_1D_" + (il+1) + "_" + pn), 1);
            }
        }
    }
    
    @Override
    public void setPlottingOptions(String key) {
        for(EmbeddedPad pad : this.getCanvas(key).getCanvasPads()) {
            if(pad.getDatasetPlotters().get(0).getDataSet() instanceof H1F) {
                pad.getAxisY().setLog(true);
                pad.setOptStat("100001");
            }
            else {
                pad.getAxisZ().setLog(true);
            }
        }
    }
   
}
