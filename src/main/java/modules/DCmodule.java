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
public class DCmodule extends Module {
    
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
            H1F hi_energy  = histo1D("hi_energy_" + name, name + " Hit Energy",   "Counts", 100,   0, EMAX[i], col);
            H1F hi_time    = histo1D("hi_time_" + name,   name + " Hit Time",     "Counts", 10+(int) TMAX[i], -10, TMAX[i], col);
            H1F hi_resi    = histo1D("hi_resi_" + name,   name + " Hit Residual (um)", "Counts", 100, -RMAX[i], RMAX[i], col);

            dg.addDataSet(hi_energy, 0 + i*3);
            dg.addDataSet(hi_time,   1 + i*3);
            dg.addDataSet(hi_resi,   2 + i*3);
        }
        return dg;
    }

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


    @Override
    public void createHistos() {
//        this.getHistos().put("Hits",           this.dcHitVertices(44));
//        this.getHistos().put("HitsOnTrack",    this.dcHitVertices(3));
//        this.getHistos().put("HitsNotOnTrack", this.dcHitVertices(44));
        this.getHistos().put("Sector Occupancy",  this.occupancies());
        this.getHistos().put("SVTOnTrack",     this.occupancies());
    }
    
    @Override
    public void fillHistos(Event event) {
        if(event.getHits(DetectorType.DC)!=null) {
            this.fillOccupancies(this.getHistos().get("Sector Occupancy"), event.getHits(DetectorType.DC));
        }
    }
    
    public void fillOccupancies(DataGroup group, List<Hit> hits) {
        for(Hit hit : hits) {
            group.getH2F("hi_occ_sector" + hit.getSector()).fill(hit.getComponent(), hit.getLayer());
        }
    }

    public void fillSVTGroup(DataGroup group, List<Hit> hits) {
        for(Hit hit : hits) {
            if(hit.getType()==DetectorType.BST) 
                group.getH1F("hi_occ_layer"+hit.getLayer()).fill(hit.getComponent()+Constants.SVTSTRIPS*(hit.getSector()-1));  
        }
    }
        
    @Override
    public void analyzeHistos() {
        this.normalizeToEvents(this.getHistos().get("Sector Occupancy"));
    }
    
    @Override
    public void drawHistos() {
        this.addCanvas("Sector Occupancy", "SVT");
        this.getCanvas().getCanvas("Sector Occupancy").draw(this.getHistos().get("Sector Occupancy"));
        super.setPlottingOptions("Sector Occupancy");
    }
       
    @Override
    public void setPlottingOptions(String name) {
        this.getCanvas().getCanvas(name).setGridX(false);
        this.getCanvas().getCanvas(name).setGridY(false);
        for(int i=0; i<3; i++) {
            this.getCanvas().getCanvas(name).getPad(i*3+0).getAxisY().setLog(true);
            this.getCanvas().getCanvas(name).getPad(i*3+1).getAxisY().setLog(true);
        }
    }
}
