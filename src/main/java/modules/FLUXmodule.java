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
    
    public DataGroup fluxes() {
        DataGroup dg = new DataGroup(3,2);
        
        for(int il=0; il<2; il++) {
            for (int ip=0; ip<PNAMES.length; ip++) {
                H1F hi_all = histo1D("hi_allflux_1D_"+(il+1) + "_" + PNAMES[ip], "", "#theta (deg)", "Flux at R=" + R[il] + "cm [Hz/cm2] ", (int) (100/DTHETA), 0, 100, 0); 
                this.setHistoAttr(hi_all, ip<5 ? ip+1 : ip+3);
                H1F hi_fwd = histo1D("hi_fwdflux_1D_"+(il+1) + "_" + PNAMES[ip], "", "#theta (deg)", "Flux at R=" + R[il] + "cm [Hz/cm2] ", (int) (100/DTHETA), 0, 100, 0); 
                this.setHistoAttr(hi_fwd, ip<5 ? ip+1 : ip+3);
                dg.addDataSet(hi_all, 0 + il*3);
                dg.addDataSet(hi_fwd, 1 + il*3);
            }
            H2F hi_flux = histo2D("hi_allflux_2D_"+(il+1), "Flux (kHz/cm2)", "x (cm)", "y (cm) ", (int) R[il]*2, -R[il], R[il], (int) R[il]*2, -R[il], R[il]); 
            dg.addDataSet(hi_flux, 2 + il*3);
        }
        return dg;
    }
  
    @Override
    public void createHistos() {
        this.getHistos().put("Fluxes", this.fluxes());
    }
    
    @Override
    public void fillHistos(Event event) {
        if (event.getHits(DetectorType.TARGET) != null) {
            this.fillFluxes(this.getHistos().get("Fluxes"), event.getHits(DetectorType.TARGET));
        }
    }
    
    public void fillFluxes(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            Vector3D proj = hit.getTrue().getPosition().toVector3D().asUnit();
            double theta = hit.getTrue().getPosition().toVector3D().theta();
            if(theta>Math.toRadians(70)) continue;
            double domega = 2*Math.PI*Math.sin(theta)*Math.toRadians(DTHETA);
            int il = hit.getComponent()/10;
            boolean fwd = hit.getTrue().getPosition().toVector3D().dot(hit.getTrue().getMomentum())>0;
            group.getH2F("hi_allflux_2D_" + (il+1)).fill(proj.x()*R[il],proj.y()*R[il]);
            group.getH1F("hi_allflux_1D_" + (il+1) + "_all").fill(Math.toDegrees(theta), 1/domega/R[il]/R[il]);
            if(fwd)
                group.getH1F("hi_fwdflux_1D_" + (il+1) + "_all").fill(Math.toDegrees(theta), 1/domega/R[il]/R[il]);                
            if(this.pidToName(Math.abs(Math.abs(hit.getTrue().getPid())))!=null) {
                group.getH1F("hi_allflux_1D_" + (il+1) + "_" + this.pidToName(Math.abs(hit.getTrue().getPid()))).fill(Math.toDegrees(theta), 1/domega/R[il]/R[il]);
                if(fwd)
                    group.getH1F("hi_fwdflux_1D_" + (il+1) + "_" + this.pidToName(Math.abs(hit.getTrue().getPid()))).fill(Math.toDegrees(theta), 1/domega/R[il]/R[il]);
            }
        }
    }
    
    @Override
    public void analyzeHistos() {
        for(int il=0; il<2; il++) {
            this.normalizeToTime(this.getHistos().get("Fluxes").getH2F("hi_allflux_2D_" + (il+1)), 1);
            for (String pn : PNAMES) {
                this.normalizeToTime(this.getHistos().get("Fluxes").getH1F("hi_allflux_1D_" + (il+1) + "_" + pn), 1);
                this.normalizeToTime(this.getHistos().get("Fluxes").getH1F("hi_fwdflux_1D_" + (il+1) + "_" + pn), 1);
                this.getHistos().get("Fluxes").getH1F("hi_allflux_1D_" + (il+1) + "_" + pn).setOptStat("1000001");
                this.getHistos().get("Fluxes").getH1F("hi_fwdflux_1D_" + (il+1) + "_" + pn).setOptStat("1000001");
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
