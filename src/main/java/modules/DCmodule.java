package modules;

import analysis.Constants;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import objects.Hit;
import java.text.DecimalFormat;
import analysis.Module;
import objects.Event;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.LatexText;

/**
 *
 * @author devita
 */
public class DCmodule extends Module {
    
    private static final int NREGIONS = 3;
    private static final int NSECTORS = 6;
    private static final int NLAYERS = 36;
    private static final int NWIRES = 112;

    public DCmodule(int residualScale) {
        super(DetectorType.DC);
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

    

    /* occupacy layer vs wire for each sector - need to normalize wrt number of region (3) */
    public DataGroup occupancies() {
        DataGroup dg = new DataGroup(3,2);
        for(int is=0; is<NSECTORS; is++) {
            int sector = is+1;
            String name = "sector"+sector;
            H2F hi_occ  = histo2D("hi_occ_" + name, "Wire",   "layer", NWIRES, 1, NWIRES+1, NLAYERS, 1, NLAYERS+1);           
            dg.addDataSet(hi_occ, 0 + is);
        }
        return dg;
    }

    /*  occupancy vs sector for each region - need to normalize wrt number of layer (12) and wires (112)    */
    public DataGroup occupancy_region() {
        DataGroup dg = new DataGroup(1, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            String name = "region" + region;
            H1F hi_occ = histo1D("hi-occ-" + name, name, "Sector", "Occupancy[%] ", NSECTORS, 0.5, NSECTORS + 0.5, 0);
            hi_occ.setLineColor(region + 1);
            hi_occ.setLineWidth(4);

            dg.addDataSet(hi_occ, 0);

        }
        return dg;
    }
    
        /*  histos to understand origin of BG  */
    public DataGroup[] origin_bg() {
        DataGroup[] dg = new DataGroup[2];

        dg[0] = new DataGroup(3, 3);
        dg[1] = new DataGroup(3, 3);
  
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            String name = "region" + region;
            H2F hi_bg_origin_rz = histo2D("hi_bg_origin_rz_" + name, "Vz(mm)", "r(mm) ", 200, -500., 5500., 200, 0., 5500.);
            H2F hi_bg_origin_rzE = histo2D("hi_bg_origin_rzE_" + name, "z(mm)", "r(mm) ", 200, -500., 5500., 200, 0., 5500.);
            
            dg[0].addDataSet(hi_bg_origin_rz, 0 + ir);
            dg[1].addDataSet(hi_bg_origin_rzE, 0 + ir);

            H2F hi_bg_origin_xy = histo2D("hi_bg_origin_xy_" + name, "Vx(mm)", "Vy(mm) ", 200, -1000.,1000., 200, -1000.,1000.);
            H2F hi_bg_origin_xyE = histo2D("hi_bg_origin_xyE_" + name, "Vx(mm)", "Vy(mm) ", 200, -1000.,1000., 200, -1000.,1000.);
            
            dg[0].addDataSet(hi_bg_origin_xy, 3 + ir);
            dg[1].addDataSet(hi_bg_origin_xyE, 3 + ir);
            String name1 =null;
            
            for (int i = 0; i < 7; i++) {
                name = "region" + region;
                if (i == 0) {
                    name1 = "all";
                    name = name + "-"+name1;
                }
                if (i == 1) {
                    name1 = "electron";
                    name = name + "-" + name1;
                }
                if (i == 2) {
                    name1 = "gamma";
                    name = name + "-" + name1;
                }
                if (i == 3) {
                    name1 = "neutron";
                    name = name + "-" + name1;
                }
                if (i == 4) {
                    name1 = "proton";
                    name = name + "-" + name1;
                }
                if (i == 5) {
                    name1 = "pion";
                    name = name + "-" + name1;
                }
                if (i == 6) {
                    name1 = "others";
                    name = name + "-" + name1;
                }
                H1F hi_bg_z = histo1D("hi-bg-z-" + name, name1, "Vz(mm)", "Rate [MHz] ", 200, -500., 5500., 0);
                       
                H1F hi_bg_E = histo1D("hi-bg-E-" + name, name1, "p(MeV)", "Rate [MHz] ", 200, 0., 200., 0);
                
                if (i < 5) {
                    this.setHistoAttr(hi_bg_z, i + 1);
                    this.setHistoAttr(hi_bg_E, i + 1);
                } else {
                    this.setHistoAttr(hi_bg_z, i + 3);
                    this.setHistoAttr(hi_bg_E, i + 3);
                }
                   
                dg[0].addDataSet(hi_bg_z, 6 + ir);
                dg[1].addDataSet(hi_bg_E, 6 + ir);
            }

        }
            
        return dg;
    }


        
    public String PID_name(int pid) {
        String name = null;
        switch (pid) {
            case 11:
            case -11:
                name = "electron";
                break;
            case 22:
                name = "gamma";
                break;
            case 2112:
                name = "neutron";
                break;
            case 2212:
                name = "proton";
                break;
            case 211:
            case -211:
                name = "pion";
                break;
            default:
                name = "others";
                break;
        }
        return name;
    }

    
    public DataGroup pos_bg() {
        DataGroup dg = new DataGroup(3, 2);

        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            double dx = 2000+ir*1000;
            String name = "region" + region;
            H2F hi_posZ_posYR = histo2D("hi-posZ-posR-" + name, "posZ [mm]", "posR [mm]",300, 0., 6000., 300, 0, Math.sqrt(2)*dx);
            H2F hi_posX_posY = histo2D("hi-posX-posY-" + name, "poX [mm]", "posY [deg]",200, -dx , dx, 200, -dx, dx);
            dg.addDataSet(hi_posZ_posYR, ir);
            dg.addDataSet(hi_posX_posY, ir+3);
        }

        return dg;
    }
 
    
  
    @Override
    public void createHistos() {

        this.getHistos().put("Sector Occupancy", this.occupancies());
        this.getHistos().put("Region Occupancy", this.occupancy_region());
        this.getHistos().put("Origin of Bg", this.origin_bg()[0]);
        this.getHistos().put("Origin of Bg - Energy", this.origin_bg()[1]);
        this.getHistos().put("Position of Bg", this.pos_bg());
        
    }
    
    
    
    @Override
    public void fillHistos(Event event) {
        if (event.getHits(DetectorType.DC) != null) {
            this.fillOccupancies(this.getHistos().get("Sector Occupancy"), event.getHits(DetectorType.DC));
            this.fillOccupancy_region(this.getHistos().get("Region Occupancy"), event.getHits(DetectorType.DC));
            this.fillOrigin(this.getHistos().get("Origin of Bg"), event);
            this.fillOriginE(this.getHistos().get("Origin of Bg - Energy"), event);
            this.fillPosition(this.getHistos().get("Position of Bg"), event);
        }

    }
    
    public void fillOccupancies(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            group.getH2F("hi_occ_sector" + hit.getSector()).fill(hit.getComponent(), hit.getLayer());
        }
    }

    public void fillOccupancy_region(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            int region = (hit.getLayer() - 1) / 12 + 1;
            group.getH1F("hi-occ-region" + region).fill(hit.getSector());
        }
    }
    
        
    
    public void fillOrigin(DataGroup group, Event event){
         
        List<Hit> hits = event.getHits(DetectorType.DC);
        for(True t : event.getTrues(DetectorType.DC)){
         
          Hit hit = hits.get(t.getHitn()-1);
          int region = (hit.getLayer()-1)/12+1;
          String name = "region"+region;
          double r = Math.sqrt(t.getVertex().toVector3D().x()*t.getVertex().toVector3D().x() + t.getVertex().toVector3D().y()*t.getVertex().toVector3D().y() );
          group.getH2F("hi_bg_origin_rz_" + name).fill(t.getVertex().z(), r);
          
          if(t.getVertex().z()>1000*(region)) group.getH2F("hi_bg_origin_xy_" + name).fill(t.getVertex().x(), t.getVertex().y());
         
          group.getH1F("hi-bg-z-" + name+ "-all").fill(t.getVertex().z());
          
          
          String name_pid = PID_name(t.getPid());
     
          group.getH1F("hi-bg-z-" + name+ "-" +name_pid).fill(t.getVertex().z());
          
          }
   
    }
    
        public void fillOriginE(DataGroup group, Event event){
        List<Hit> hits = event.getHits(DetectorType.DC);
        for(True t : event.getTrues(DetectorType.DC)){
          Hit hit = hits.get(t.getHitn()-1);
          int region = (hit.getLayer()-1)/12+1;
          String name = "region"+region;
          double r = Math.sqrt(t.getVertex().toVector3D().x()*t.getVertex().toVector3D().x() + t.getVertex().toVector3D().y()*t.getVertex().toVector3D().y() );
          group.getH2F("hi_bg_origin_rzE_" + name).fill(t.getVertex().z(), r, t.getMomentum().mag()*1000);
          
           if(t.getVertex().z()>1000*(region)) group.getH2F("hi_bg_origin_xyE_" + name).fill(t.getVertex().x(), t.getVertex().y(),t.getMomentum().mag() );
          
          group.getH1F("hi-bg-E-" + name+ "-all").fill(t.getMomentum().mag()*1000);
          
          String name_pid = PID_name(t.getPid());
     
          group.getH1F("hi-bg-E-" + name+ "-" +name_pid).fill(t.getMomentum().mag()*1000);
          
        }           
    
    }
        
        
    public void fillPosition(DataGroup group, Event event) {
        List<Hit> hits = event.getHits(DetectorType.DC);

        for (True t : event.getTrues(DetectorType.DC)) {
            Hit hit = hits.get(t.getHitn() - 1);
            int region = (hit.getLayer() - 1) / 12 + 1;

            String name = "region" + region;
          
            double r = Math.sqrt(t.getPosition().y() * t.getPosition().y() + t.getPosition().x() * t.getPosition().x());
            if (t.getVertex().z() > -150 && t.getVertex().z() < 100 && (t.getPid() == 11 || t.getPid() == -11)) {
                group.getH2F("hi-posZ-posR-" + name).fill(t.getPosition().z(), r);

                group.getH2F("hi-posX-posY-" + name).fill(t.getPosition().x(), t.getPosition().y());

            }
        }

    }   
    
    
    @Override
    public void analyzeHistos() {
        this.normalizeToEvents(this.getHistos().get("Sector Occupancy"));
        this.normalize(this.getHistos().get("Sector Occupancy"), 1./100);
        double norm = 112 * 12 / 100;
        this.normalizeToEvents(this.getHistos().get("Region Occupancy"));
        this.normalize(this.getHistos().get("Region Occupancy"), norm);
        
        this.fitDataGroup(this.getHistos().get("Region Occupancy"));
   
       
        this.normalizeHisto(this.getHistos().get("Origin of Bg - Energy"), this.getHistos().get("Origin of Bg"));
        this.normalizeToTime(this.getHistos().get("Origin of Bg"));
        this.normalizeToTime(this.getHistos().get("Origin of Bg - Energy"));
        
    }
    
    public void normalizeHisto(DataGroup dg1, DataGroup dg2) {

        int nrow = dg1.getRows();
        int ncol = dg1.getColumns();
        
        for (int i = 0; i < nrow * ncol; i++) {
            int idx =0;
            for (IDataSet ds1 : dg1.getData(i)) {
                if (ds1 instanceof H2F) {
                    H2F h1 = (H2F) ds1;
                    IDataSet ds2 = dg2.getData(i).get(idx);
                    if (ds2 instanceof H2F) {
                         H2F h2 = (H2F) ds2;
                        h1.divide(h2);      
                    }
                }
           idx = idx+1;
            }
        }
    }
    
    
    @Override
    public void fitDataGroup(DataGroup dg) {
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            F1D f = fitPol0(dg.getH1F("hi-occ-region" + region));
            dg.addDataSet(f, 0);
        }
    }
    
    @Override
    public void drawHistos() {
        this.addCanvas("Sector Occupancy", "Region Occupancy", "Origin of Bg", "Origin of Bg - Energy", "Position of Bg");

        this.getCanvas().getCanvas("Sector Occupancy").draw(this.getHistos().get("Sector Occupancy"));
        super.setPlottingOptions("Sector Occupancy");

        this.getCanvas().getCanvas("Region Occupancy").draw(this.getHistos().get("Region Occupancy"));
        super.setPlottingOptions("Region Occupancy");
        this.setTextCanvas("Region Occupancy");

        this.getCanvas().getCanvas("Origin of Bg").draw(this.getHistos().get("Origin of Bg"));
        super.setPlottingOptions("Origin of Bg");
        this.setLogZ("Origin of Bg");
        this.setLegend("Origin of Bg");  

        this.getCanvas().getCanvas("Origin of Bg - Energy").draw(this.getHistos().get("Origin of Bg - Energy"));
        this.setPlottingOptions("Origin of Bg - Energy");
        this.setLogZ("Origin of Bg - Energy");
        this.setLegend("Origin of Bg - Energy"); 
        
        this.getCanvas().getCanvas("Position of Bg").draw(this.getHistos().get("Position of Bg"));
        super.setPlottingOptions("Position of Bg");
        this.setLogZ("Position of Bg");
        

      //  super.setTextCanvas("Origin of Bg");
    }
    
    public void setLegend(String text) {
        DataGroup dg = this.getHistos().get(text);
        int nx = dg.getColumns();
        int ny = dg.getRows();
        for (int i = 0; i < nx * ny; i++) {
            if(i>2*nx-1){
            this.getCanvas().getCanvas(text).getPad(i).setLegend(true);
            this.getCanvas().getCanvas(text).getPad(i).setLegendPosition(65, 100);
            }
        }
    }
   
    @Override
    public void setTextCanvas(String Canvas) {

        DataGroup dg = this.getHistos().get(Canvas);
        int nx = dg.getColumns();
        int ny = dg.getRows();
        LatexText LT = null;
        for (int i = 0; i < nx * ny; i++) {
            List<IDataSet> ds = dg.getData(i);
            for (IDataSet d : ds) {
                if (d instanceof F1D) {
                    
                    double par = ((F1D)d).getParameter(0);
                    String reg = this.extractNumber(((F1D)d).getName());
                    DecimalFormat df = new DecimalFormat("#.000");
                    String text = "Region" + reg + " " + df.format(par)+"%";
                    int col = Integer.parseInt(reg)+1;
                    LT = setLatex(text, col);
                    this.getCanvas().getCanvas(Canvas).draw(LT);
                }
            }
        }

    }
    
    
        
    @Override
    public void setPlottingOptions(String name) {
        this.getCanvas().getCanvas(name).setGridX(false);
        this.getCanvas().getCanvas(name).setGridY(false);
        for (int i = 0; i < 3; i++) {
            //this.getCanvas().getCanvas(name).getPad(i*3+0).getAxisY().setLog(true);
            //this.getCanvas().getCanvas(name).getPad(i*3+1).getAxisY().setLog(true);
            this.getCanvas().getCanvas(name).getPad(i * 1 + 6).getAxisY().setLog(true);
        }
    }

    @Override
    public void normalizeToTime(DataGroup dg) {
        double factor = this.getNevents() * (Constants.getTimeWindow() * 1E-9) / 1E-6; //Mz

        int nrow = dg.getRows();
        int ncol = dg.getColumns();
        for (int i = 0; i < nrow * ncol; i++) {
            for (IDataSet ds : dg.getData(i)) {
                if (ds instanceof H1F) {
                    H1F h = (H1F) ds;
                    h.divide(factor);
                }
            }
        }
    
    
    
    }



}
