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
    private static final double[] R = {30, 60, 600}; // cm
    private static final double[] tmin = {30,  0,  0}; // cm
    private static final double[] tmax = {90, 30, 40}; // cm
    private static final double THRESHOLD = 0;
    
    public FLUXmodule() {
        super(DetectorType.TARGET);
    }
    
    public DataGroup rates() {
        DataGroup dg = new DataGroup(3,2);
        
        String[] title = {"Central", "Calorimeter", "Forward carriage"};
        for(int il=0; il<3; il++) {
            for (int ip=0; ip<PNAMES.length; ip++) {
                H1F hi_all = histo1D("hi_all_1D_"+(il+1) + "_" + PNAMES[ip], title[il], "#theta (deg)", "Flux [Hz/d#Omega] ", (int) ((tmax[il]-tmin[il])/DTHETA), tmin[il], tmax[il], 0); 
                this.setHistoAttr(hi_all, ip<5 ? ip+1 : ip+3);
                H1F hi_bwd = histo1D("hi_bwd_1D_"+(il+1) + "_" + PNAMES[ip], title[il], "#theta (deg)", "Flux [Hz/d#Omega] ", (int) ((tmax[il]-tmin[il])/DTHETA), tmin[il], tmax[il], 0); 
                this.setHistoAttr(hi_bwd, ip<5 ? ip+1 : ip+3);
                dg.addDataSet(hi_all, il + 0);
                dg.addDataSet(hi_bwd, il + 3);
            }
        }
        return dg;
    }
  
    public DataGroup fluxes() {
        DataGroup dg = new DataGroup(3,3);
        
        for(int il=0; il<3; il++) {
            String xt = "x (cm)";
            String yt = "y (cm)";
            double xmin = -R[il]*(0.4+0.2*il);
            double xmax =  R[il]*(0.4+0.2*il);
            double ymin = -R[il]*(0.4+0.2*il);
            double ymax =  R[il]*(0.4+0.2*il);
            if(il==0) {
                yt = "R#phi (cm)";
                xt = "z (cm)";
                xmin = 0;
                xmax = R[il]*2;
                ymin = -R[il]*3.5;
                ymax =  R[il]*3.5;
            }
            H2F hi_all = histo2D("hi_all_2D_"+(il+1), "Total Flux (Hz/cm2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax); 
            H2F hi_crg = histo2D("hi_crg_2D_"+(il+1), "Charged Flux (Hz/cm2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax); 
            H2F hi_neu = histo2D("hi_neu_2D_"+(il+1), "Neutron Flux (Hz/cm2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax);  
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
            if(hit.getTrue().getKinEnergy()<THRESHOLD) continue;
            double theta = hit.getTrue().getPosition().toVector3D().theta();
            double domega = 2*Math.PI*Math.sin(theta)*Math.toRadians(DTHETA);
            int il = hit.getTrue().getPosition().z()<500 ? 0 : 1+hit.getComponent()/10;
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
            if(hit.getTrue().getKinEnergy()<THRESHOLD) continue;
            Vector3D proj = hit.getTrue().getPosition().toVector3D().asUnit();
            double theta = hit.getTrue().getPosition().toVector3D().theta();
            if(theta>Math.toRadians(70)) continue;
            int pid = Math.abs(hit.getTrue().getPid());
            int il = hit.getTrue().getPosition().z()<500 ? 0 : 1+hit.getComponent()/10;
            double x = il==0 ? hit.getTrue().getPosition().z() : proj.x()*R[il];
            double y = il==0 ? proj.phi()*R[il] : proj.y()*R[il];
            group.getH2F("hi_all_2D_" + (il+1)).fill(x, y);
            if(pid==11 || pid==2212 || pid==211) {
                group.getH2F("hi_crg_2D_" + (il+1)).fill(x, y);
            }
            else if(pid==2112) {
                group.getH2F("hi_neu_2D_" + (il+1)).fill(x, y);
            }
        }
    }
    
    @Override
    public void analyzeHistos() {
        System.out.println("aaa");
        for(int il=0; il<3; il++) {
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
                H1F h = (H1F) pad.getDatasetPlotters().get(0).getDataSet();
                pad.getAxisY().setLog(true);
                pad.getAxisY().setRange(pad.getAxisY().getMin(), 10*h.getMax());
                pad.setOptStat("100001");
            }
//           



        }
    }
   
}
