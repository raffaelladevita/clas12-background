package modules;

import java.util.List;
import objects.Hit;
import analysis.Module;
import java.util.Arrays;
import objects.Event;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.TCanvas;

/**
 *
 * @author devita
 */
public class FLUXmodule extends Module {
    
    private static final double DTHETA = 0.5;
    private static final double DR = 0.1;
    private static final double RREC = 7.5;
    private static final double RHOD = 25;
    private static final double RTRK = 40;
    private static final double RCAL = 60;
    private static final double RFC  = 600;
    private static final double[] tmin = {35,  0,  0}; // cm
    private static final double[] tmax = {75, 40, 40}; // cm
    private static final double[] THRESHOLD = {0.01, 20};
    
    public FLUXmodule() {
        super(DetectorType.TARGET);
    }
    
    private int[] getLayerIndex(double[] radiuses) {
        int[] layers = new int[radiuses.length];
        for(int i=0; i<layers.length; i++) {
            layers[i] = Math.max(0, i-(layers.length-tmax.length));
        }
        return layers;
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
  
    public DataGroup energies() {
        DataGroup dg = new DataGroup(3,2);
        
        String[] title = {"Central", "Calorimeter", "Forward carriage"};
        for(int il=0; il<3; il++) {
            for (int ip=0; ip<PNAMES.length; ip++) {
                H1F hi_all = histo1D("hi_all_1D_"+(il+1) + "_" + PNAMES[ip], title[il], "E (MeV)", "Flux (Hz/cm^2) ", 100, 0, 1000, 0); 
                this.setHistoAttr(hi_all, ip<5 ? ip+1 : ip+3);
                H1F hi_bwd = histo1D("hi_bwd_1D_"+(il+1) + "_" + PNAMES[ip], title[il], "E (MeV)", "Flux (Hz/cm^2) ", 100, 0, 1, 0); 
                this.setHistoAttr(hi_bwd, ip<5 ? ip+1 : ip+3);
                dg.addDataSet(hi_all, il + 0);
                dg.addDataSet(hi_bwd, il + 3);
            }
        }
        return dg;
    }
  
    public DataGroup origins() {
        DataGroup dg = new DataGroup(3,2);
        
        String[] title = {"Central", "Calorimeter", "Forward carriage"};
        for(int il=0; il<3; il++) {
            for (int ip=0; ip<PNAMES.length; ip++) {
                H1F hi_all = histo1D("hi_all_1D_"+(il+1) + "_" + PNAMES[ip], title[il], "z (cm)", "Flux (Hz/cm^2) ", 500, -100, 900, 0); 
                this.setHistoAttr(hi_all, ip<5 ? ip+1 : ip+3);
                H1F hi_bwd = histo1D("hi_bwd_1D_"+(il+1) + "_" + PNAMES[ip], title[il], "z (cm)", "Flux (Hz/cm^2) ", 500, -100, 900, 0); 
                this.setHistoAttr(hi_bwd, ip<5 ? ip+1 : ip+3);
                dg.addDataSet(hi_all, il + 0);
                dg.addDataSet(hi_bwd, il + 3);
            }
        }
        return dg;
    }
  
    public DataGroup fluxes(double... R) {
        DataGroup dg = new DataGroup(5,R.length);
        
        int[] layers = this.getLayerIndex(R);
        
        for(int il=0; il<R.length; il++) {
            String xt = "x (cm)";
            String yt = "y (cm)";
            String rt = "r (cm)";
            int li = layers[il];
            double xmin = -Math.floor(R[il]*Math.sin(Math.toRadians(tmax[li])));
            double xmax =  Math.floor(R[il]*Math.sin(Math.toRadians(tmax[li])));
            double ymin = -Math.floor(R[il]*Math.sin(Math.toRadians(tmax[li])));
            double ymax =  Math.floor(R[il]*Math.sin(Math.toRadians(tmax[li])));
            double rmin =  Math.floor(R[il]*Math.sin(Math.toRadians(tmin[li])));
            double rmax =  Math.ceil(R[il]*Math.sin(Math.toRadians(tmax[li])));
            if(li==0) {
                yt = "R#phi (cm)";
                xt = "z (cm)";
                rt = "z (cm)";
                xmin = Math.ceil(R[il]/Math.tan(Math.toRadians(tmax[0]-5)));
                xmax = Math.floor(R[il]/Math.tan(Math.toRadians(tmin[0]+5)));
                ymin = -Math.floor(R[il]*2*Math.PI)/2;
                ymax =  Math.floor(R[il]*2*Math.PI)/2;
                rmin = Math.ceil(R[il]/Math.tan(Math.toRadians(tmax[0]-5)));
                rmax = Math.floor(R[il]/Math.tan(Math.toRadians(tmin[0]+5)));
            }
            System.out.println((int) (ymax-ymin) + " " + ymin + " " + ymax);
            H2F hi_all_2D = histo2D("hi_all_2D_"+(il+1), "Total Flux (Hz/cm^2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax); 
            H2F hi_pho_2D = histo2D("hi_pho_2D_"+(il+1), "Photon Flux (Hz/cm^2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax); 
            H2F hi_crg_2D = histo2D("hi_crg_2D_"+(il+1), "Charged Flux (Hz/cm^2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax); 
            H2F hi_neu_2D = histo2D("hi_neu_2D_"+(il+1), "Neutron Flux (Hz/cm^2)", xt, yt, (int) (xmax-xmin), xmin, xmax, (int) (ymax-ymin), ymin, ymax);  
            H1F hi_all_1D = histo1D("hi_all_1D_"+(il+1), "Total Flux", rt, "Flux (Hz/cm^2)", (int) ((rmax-rmin)/DR), rmin, rmax, -1);
            H1F hi_pho_1D = histo1D("hi_pho_1D_"+(il+1), "Photon Flux", rt, "Flux (Hz/cm^2)", (int) ((rmax-rmin)/DR), rmin, rmax, -3); 
            H1F hi_crg_1D = histo1D("hi_crg_1D_"+(il+1), "Charged Flux", rt, "Flux (Hz/cm^2)", (int) ((rmax-rmin)/DR), rmin, rmax, -2); 
            dg.addDataSet(hi_all_2D, 0 + il*5);
            dg.addDataSet(hi_pho_2D, 1 + il*5);
            dg.addDataSet(hi_crg_2D, 2 + il*5);
            dg.addDataSet(hi_neu_2D, 3 + il*5);
            dg.addDataSet(hi_all_1D, 4 + il*5);
            dg.addDataSet(hi_pho_1D, 4 + il*5);
            dg.addDataSet(hi_crg_1D, 4 + il*5);
        }
        return dg;
    }
  
    @Override
    public void createHistos() {
        this.getHistos().put("Rates", this.rates());
        this.getHistos().put("Energies", this.energies());
        this.getHistos().put("Origins", this.origins());
        this.getHistos().put("Fluxes", this.fluxes(RREC, RHOD, RTRK, RFC));
        this.getHistos().put("Fluxes E>"+THRESHOLD[1]+" MeV", this.fluxes(RHOD, RCAL, RFC));
    }
    
    @Override
    public void fillHistos(Event event) {
        if (event.getHits(DetectorType.TARGET) != null) {
            this.fillRates(this.getHistos().get("Rates"), event.getHits(DetectorType.TARGET));
            this.fillEnergies(this.getHistos().get("Energies"), event.getHits(DetectorType.TARGET), RHOD, RTRK, RFC);
            this.fillOrigins(this.getHistos().get("Origins"), event.getHits(DetectorType.TARGET), RHOD, RTRK, RFC);
            this.fillFluxes(this.getHistos().get("Fluxes"), event.getHits(DetectorType.TARGET), THRESHOLD[0], RREC, RHOD, RTRK, RFC);
            this.fillFluxes(this.getHistos().get("Fluxes E>"+THRESHOLD[1]+" MeV"), event.getHits(DetectorType.TARGET), THRESHOLD[1], RHOD, RCAL, RFC);
        }
    }
    
    private int getFluxIndex(Hit hit) {
        int il = -1;
        double theta = hit.getTrue().getPosition().toVector3D().theta();
        if(hit.getComponent()==1 && theta<Math.toRadians(70) && theta>Math.toRadians(40)) {
            il=0;
        }            
        else if(hit.getComponent()==1 && theta<Math.toRadians(35) && theta>Math.toRadians(7.5)) {            
            il = 1;
        }
        else if(hit.getComponent()==10 && theta<Math.toRadians(35) && theta>Math.atan(0.05)) {            
            il = 2;
        }
        return il;
    }
    
    public void fillRates(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            
            if(hit.getTrue().getKinEnergy()<THRESHOLD[0]) continue;
    
            int il = this.getFluxIndex(hit);
            if(il<0) continue;
            
            double theta = hit.getTrue().getPosition().toVector3D().theta();
            double domega = 2*Math.PI*Math.sin(theta)*Math.toRadians(DTHETA);
            
            group.getH1F("hi_all_1D_" + (il+1) + "_all").fill(Math.toDegrees(theta), 1/domega);
            if(hit.getTrue().getKinEnergy()>THRESHOLD[1])
                group.getH1F("hi_bwd_1D_" + (il+1) + "_all").fill(Math.toDegrees(theta), 1/domega);                
            String pname = this.pidToName(Math.abs(Math.abs(hit.getTrue().getPid()))); 
            if(pname!=null) {
                group.getH1F("hi_all_1D_" + (il+1) + "_" + pname).fill(Math.toDegrees(theta), 1/domega);
                if(hit.getTrue().getKinEnergy()>THRESHOLD[1])
                    group.getH1F("hi_bwd_1D_" + (il+1) + "_" + pname).fill(Math.toDegrees(theta), 1/domega);
            }
        }
    }

    public void fillEnergies(DataGroup group, List<Hit> hits, double... R) {
        for (Hit hit : hits) {
                        
            int il = this.getFluxIndex(hit);
            if(il<0) continue;

            double theta = hit.getTrue().getPosition().toVector3D().theta();
            double domega = 2*Math.PI*Math.sin(theta)*Math.toRadians(DTHETA);
            double ds = domega*R[il]*R[il];
            group.getH1F("hi_all_1D_" + (il+1) + "_all").fill(hit.getTrue().getKinEnergy(), 1/ds);
            group.getH1F("hi_bwd_1D_" + (il+1) + "_all").fill(hit.getTrue().getKinEnergy(), 1/ds);
            String pname = this.pidToName(Math.abs(hit.getTrue().getPid())); 
            if(pname!=null) {
                group.getH1F("hi_all_1D_" + (il+1) + "_" + pname).fill(hit.getTrue().getKinEnergy(), 1/ds);
                group.getH1F("hi_bwd_1D_" + (il+1) + "_" + pname).fill(hit.getTrue().getKinEnergy(), 1/ds);
            }
        }
    }
    
    
    public void fillOrigins(DataGroup group, List<Hit> hits, double... R) {
        for (Hit hit : hits) {
                        
            int il = this.getFluxIndex(hit);
            if(il<0) continue;

            double theta = hit.getTrue().getPosition().toVector3D().theta();
            double domega = 2*Math.PI*Math.sin(theta)*Math.toRadians(DTHETA);
            double ds = domega*R[il]*R[il];
            group.getH1F("hi_all_1D_" + (il+1) + "_all").fill(hit.getTrue().getVertex().z(), 1/ds);
            group.getH1F("hi_bwd_1D_" + (il+1) + "_all").fill(hit.getTrue().getVertex().z(), 1/ds);                
            String pname = this.pidToName(Math.abs(hit.getTrue().getPid())); 
            if(pname!=null) {
//                if(pname=="other")
//                    System.out.println(hit.getTrue().getPid());
                group.getH1F("hi_all_1D_" + (il+1) + "_" + pname).fill(hit.getTrue().getVertex().z(), 1/ds);
                group.getH1F("hi_bwd_1D_" + (il+1) + "_" + pname).fill(hit.getTrue().getVertex().z(), 1/ds);
            }
        }
    }
    
    
    public void fillFluxes(DataGroup group, List<Hit> hits, double threshold, double... R) {
        
        int[] layers = this.getLayerIndex(R);
        
        double dr=0.1;
        for (Hit hit : hits) {

            if(hit.getTrue().getKinEnergy()<threshold) continue;

            int il = this.getFluxIndex(hit);
            if(il<0) continue;

            Vector3D proj = hit.getTrue().getPosition().toVector3D().asUnit();

            int pid = Math.abs(hit.getTrue().getPid());
            
            for(int i=0; i<layers.length; i++) {
                if(il!=layers[i]) continue;
                double radius = R[i];
                double x = il==0 ? proj.z()*radius/Math.sqrt(proj.x()*proj.x()+proj.y()*proj.y()) : proj.x()*radius;
                double y = il==0 ? Math.rint(proj.phi()*radius) : proj.y()*radius;
                double r = il==0 ? proj.z()*radius/Math.sqrt(proj.x()*proj.x()+proj.y()*proj.y()) : Math.sqrt(proj.x()*proj.x()+proj.y()*proj.y())*radius;
                double w = il==0 ? 1/(2*Math.PI*radius) : 1/(2*Math.PI*r);
                w /= dr;
                group.getH2F("hi_all_2D_" + (i+1)).fill(x, y);
                group.getH1F("hi_all_1D_" + (i+1)).fill(r, w);
                if(pid==2112) {
                    group.getH2F("hi_neu_2D_" + (i+1)).fill(x, y);
                }
                else if(pid==22) {
                    group.getH2F("hi_pho_2D_" + (i+1)).fill(x, y);
                    group.getH1F("hi_pho_1D_" + (i+1)).fill(r, w);
                }
                else if(pid==11 || pid==2212 || pid==211 || pid==321 || pid==13) {
                        group.getH2F("hi_crg_2D_" + (i+1)).fill(x, y);
                        group.getH1F("hi_crg_1D_" + (i+1)).fill(r, w);
                }
            }
        }
    }
    
    @Override
    public void analyzeHistos() {
        this.normalizeToTime(this.getHistos().get("Fluxes"), 1);
        this.normalizeToTime(this.getHistos().get("Fluxes E>"+THRESHOLD[1]+" MeV"), 1);
        this.normalizeToTime(this.getHistos().get("Rates"), 1);
        this.normalizeToTime(this.getHistos().get("Energies"), 1);
    }
    
    @Override
    public void setPlottingOptions(String key) {
        for(EmbeddedPad pad : this.getCanvas(key).getCanvasPads()) {
            if(pad.getDatasetPlotters().get(0).getDataSet() instanceof H1F) {
                double max = 0;
                double min = Double.MAX_VALUE;
                for(int i=0; i<pad.getDatasetPlotters().size(); i++) {
                    H1F h = (H1F) pad.getDatasetPlotters().get(i).getDataSet();
                    if(h.getMax()>max) max = h.getMax();
                    if(h.getIntegral()>0 &&
                       h.getIntegral()/h.getDataSize(0)<min) min = h.getIntegral()/h.getDataSize(0);  
                }
                if(max<1) max = 1;
                if(min==Double.MAX_VALUE) min = 0;
                System.out.println(min + " " + max);
                pad.getAxisY().setLog(true);
                pad.getAxisY().setRange(min/10, 10*max);
//                pad.setOptStat("100001");
            }
            else {
                H2F h = (H2F) pad.getDatasetPlotters().get(0).getDataSet();
                double integral = 0;
                for(int b=0; b<h.getDataBufferSize(); b++)
                    integral += h.getDataBufferBinAsDouble(b);
                System.out.println(h.getName() + " " + integral);

            }
        }
        if(key.contains("MeV")) {            

            TCanvas recoil = new TCanvas("Recoil", 1150, 600);
            recoil.divide(2,1);
            recoil.getCanvas().setGridX(false);
            recoil.getCanvas().setGridY(false);
            recoil.getCanvas().setCanvasPads(Arrays.asList(this.getCanvas().getCanvas("Fluxes").getPad(4),
                                                           this.getCanvas().getCanvas("Fluxes").getPad(0)));
            for(EmbeddedPad pad : recoil.getCanvas().getCanvasPads()) 
                pad.setTitle("   ");
       
            TCanvas hodoscope = new TCanvas("Hodoscope", 1500, 600);
            hodoscope.divide(2,1);
            hodoscope.getCanvas().setGridX(false);
            hodoscope.getCanvas().setGridY(false);
            hodoscope.getCanvas().setCanvasPads(Arrays.asList(this.getCanvas().getCanvas("Rates").getPad(0),
                                                              this.getCanvas().getCanvas("Fluxes").getPad(6)));
            hodoscope.getCanvas().getPad(0).setTitle("Recoil");
            
            TCanvas vtracker = new TCanvas("Trk", 1150, 600);
            vtracker.divide(2,1);
            vtracker.getCanvas().setGridX(false);
            vtracker.getCanvas().setGridY(false);
            vtracker.getCanvas().setCanvasPads(Arrays.asList(this.getCanvas().getCanvas("Fluxes").getPad(14),
                                                             this.getCanvas().getCanvas("Fluxes").getPad(10)));
            for(EmbeddedPad pad : vtracker.getCanvas().getCanvasPads()) 
                pad.setTitle("   ");

            TCanvas cal = new TCanvas("Cal", 1150, 600);
            cal.divide(2,1);
            cal.getCanvas().setGridX(false);
            cal.getCanvas().setGridY(false);
            cal.getCanvas().setCanvasPads(Arrays.asList(this.getCanvas().getCanvas("Fluxes E>"+THRESHOLD[1]+" MeV").getPad(9),
                                                        this.getCanvas().getCanvas("Fluxes E>"+THRESHOLD[1]+" MeV").getPad(5)));
            for(EmbeddedPad pad : cal.getCanvas().getCanvasPads()) 
                pad.setTitle("   ");
        }      
    }
   
    private void setZMax(EmbeddedCanvas canvas) {
        
        double zmin = Double.MAX_VALUE;
        double zmax = 0;
        for(EmbeddedPad pad : canvas.getCanvasPads()) {
            zmin = Math.min(zmin,pad.getDatasetPlotters().getFirst().getDataRegion().getDimension(2).getMin());
            zmax = Math.max(zmax,pad.getDatasetPlotters().getFirst().getDataRegion().getDimension(2).getMax());
        }
        
        for(EmbeddedPad pad : canvas.getCanvasPads()) {
            pad.getAxisZ().setRange(zmin, zmax);
        }
    }
}
