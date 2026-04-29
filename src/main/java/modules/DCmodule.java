package modules;

import analysis.Constants;
import java.util.List;
import objects.Hit;
import analysis.Module;
import java.util.ArrayList;
import objects.Event;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.LatexText;
import org.jlab.groot.ui.TCanvas;

/**
 *
 * @author devita
 */
public class DCmodule extends Module {

    private static final int NREGIONS = 3;
    private static final int NSECTORS = 6;
    private static final int NLAYERS = 36;
    private static final int NWIRES = 112;

    private static final double[] RWINDOWS = {500, 1400, 1200};
    private static final double[] DR = {2500, 4000, 5500};
    private static final double[] DZ = {3500, 5000, 6500};    
    
    public DCmodule() {
        super(DetectorType.DC);
    }

    /* occupacy layer vs wire for each sector - need to normalize wrt number of region (3) */
    public DataGroup occupancies() {
        DataGroup dg = new DataGroup(3, 2);
        for (int is = 0; is < NSECTORS; is++) {
            int sector = is + 1;
            String name = "sector" + sector;
            H2F hi_occ = histo2D("hi_occ_" + name, "Wire", "Layer", NWIRES, 1, NWIRES + 1, NLAYERS, 1, NLAYERS + 1);
            dg.addDataSet(hi_occ, 0 + is);
        }
        return dg;
    }

    /*  occupancy vs sector for each region - need to normalize wrt number of layer (12) and wires (112)    */
    public DataGroup occupancy_region() {
        DataGroup dg = new DataGroup(2, 2);
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            String name = "region" + region;
            H1F hi_occ = histo1D("hi_occ_" + name, " ", "Sector", "Occupancy[%] ", NSECTORS, 0.5, NSECTORS + 0.5, 0);
            hi_occ.setLineColor(region + 1);
            hi_occ.setLineWidth(2);
            H1F hi_wire = histo1D("hi_wire_" + name, " ", "Wire", "Occupancy[%] ", NWIRES, 0.5, NWIRES + 0.5, 0);
            hi_wire.setLineColor(region + 1);
            hi_wire.setLineWidth(2);
            dg.addDataSet(hi_occ, 0);
            dg.addDataSet(hi_wire, 1);
        }
        H1F hi_lay = histo1D("hi_layer", " ", "Layer", "Occupancy[%] ", NLAYERS, 0.5, NLAYERS + 0.5, 0);
        hi_lay.setLineWidth(2);
        H1F hi_sl = histo1D("hi_sl", " ", "Superlayer", "Occupancy[%] ", NREGIONS*2, 0.5, NREGIONS*2 + 0.5, 0);
        hi_sl.setLineWidth(2);
        dg.addDataSet(hi_lay, 2);
        dg.addDataSet(hi_sl, 3);
        return dg;
    }

    /*  histos to understand origin of BG  */
    public DataGroup[] origin_bg() {
        DataGroup[] dg = new DataGroup[2];

        for(int i=0; i<dg.length; i++) {
                    
            dg[i] = new DataGroup(3, 3);
            
            
            for (int ir = 0; ir < NREGIONS; ir++) {
                int region = ir + 1;
                H2F hi_bg_origin_rz = histo2D("hi_bg_origin_rz_region" + region, "Vz(mm)", "r(mm) ", 200, -500., DZ[ir], 200, 0., DR[ir]);
                H2F hi_bg_origin_xy = histo2D("hi_bg_origin_xy_region" + region, "Vx(mm)", "Vy(mm) ", 200, -DR[ir], DR[ir], 200, -DR[ir], DR[ir]);
                dg[i].addDataSet(hi_bg_origin_xy, 0 + ir);
                dg[i].addDataSet(hi_bg_origin_rz, 3 + ir);

                double min = -500;
                double max = DZ[ir];
                String name = "Vz(mm)";
                if(i==1) {
                    min = 0;
                    max = 200;
                    name = "E(MeV)";
                }
                for (int ip=0; ip<PNAMES.length; ip++) {
                    H1F hi_bg = histo1D("hi_bg_region" + region + "_" + PNAMES[ip], PNAMES[ip], name, "Rate [MHz] ", 200, min, max, 0);   
                    this.setHistoAttr(hi_bg, ip<5 ? ip+1 : ip+3);
                    dg[i].addDataSet(hi_bg, 6 + ir);
                }
            }

        }

        return dg;
    }

    public DataGroup sectorBG() {
        DataGroup dg = new DataGroup(2,3);
        
        for (int is = 0; is < NSECTORS; is++) {
            int sector = is + 1;
            for (int ip=0; ip<PNAMES.length; ip++) {
                H1F hi_bg = histo1D("hi_bg_r1_s" + sector + "_" + PNAMES[ip], "R1S" +sector + "-" + PNAMES[ip], "Vz(mm)", "Rate [MHz] ", 200, -500, DZ[0], 0);   
                this.setHistoAttr(hi_bg, ip<5 ? ip+1 : ip+3);
                dg.addDataSet(hi_bg, is);
            }
        }  
        return dg;
    }

    public DataGroup pos_bg() {
        DataGroup dg = new DataGroup(3, 2);

        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            String name = "region" + region;
            H2F hi_posZ_posYR = histo2D("hi-posZ-posR-" + name, "posZ [mm]", "posR [mm]", 300, 0., DZ[ir], 300, 0, DR[ir]);
            H2F hi_posX_posY = histo2D("hi-posX-posY-" + name, "poX [mm]", "posY [mm]", 200, -DR[ir], DR[ir], 200, -DR[ir], DR[ir]);
            dg.addDataSet(hi_posZ_posYR, ir);
            dg.addDataSet(hi_posX_posY, ir + 3);
        }

        return dg;
    }

    @Override
    public void createHistos() {

        this.getHistos().put("Sector Occupancy", this.occupancies());
        this.getHistos().put("Region Occupancy", this.occupancy_region());
        this.getHistos().put("Origin of Bg", this.origin_bg()[0]);
        this.getHistos().put("Origin of Bg - Energy", this.origin_bg()[1]);
        this.getHistos().put("Origin of Bg - Sector", this.sectorBG());
        this.getHistos().put("Position of Bg", this.pos_bg());

    }

    @Override
    public void fillHistos(Event event) {
        List<Hit> allhits = event.getHits(DetectorType.DC);
        if (allhits!=null) {
            List<Hit> hits = new ArrayList<>();
            for(Hit h : allhits) {
                if(h.getTrue()==null || h.getTrue().getEdep()<50E-6)
                    continue;
                hits.add(h);
            }
            this.fillOccupancies(this.getHistos().get("Sector Occupancy"), hits);
            this.fillOccupancy_region(this.getHistos().get("Region Occupancy"), hits);
            this.fillOrigin(this.getHistos().get("Origin of Bg"), hits, false);
            this.fillOrigin(this.getHistos().get("Origin of Bg - Energy"), hits, true);
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector"), hits);
            this.fillPosition(this.getHistos().get("Position of Bg"), hits);
        }

    }

    private double readoutWeight(int layer) {
        if(Constants.getTimeWindow()<0)
            return 1;
        else
            return RWINDOWS[(layer-1)/12]/Constants.getTimeWindow();
    }
    public void fillOccupancies(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
//            System.out.println(hit.getTrue().getEdep() + " " + hit.getTrue().getTime()+ " " + hit.getTDC());
            group.getH2F("hi_occ_sector" + hit.getSector()).fill(hit.getComponent(), hit.getLayer(),this.readoutWeight(hit.getLayer()));
        }
    }

    public void fillOccupancy_region(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            int region = (hit.getLayer() - 1) / 12 + 1;
            group.getH1F("hi_occ_region" + region).fill(hit.getSector(),this.readoutWeight(hit.getLayer()));
            group.getH1F("hi_wire_region" + region).fill(hit.getComponent(),this.readoutWeight(hit.getLayer())*NWIRES/NSECTORS);
            group.getH1F("hi_layer").fill(hit.getLayer(),this.readoutWeight(hit.getLayer())*NLAYERS/NREGIONS/NSECTORS);
            group.getH1F("hi_sl").fill((hit.getLayer()-1)/6+1,this.readoutWeight(hit.getLayer())/NSECTORS*2);
        }
    }

    public void fillOrigin(DataGroup group, List<Hit> hits, boolean energyWeight) {

        for (Hit hit : hits) {
            
            True t = hit.getTrue();
            if(t==null) continue;
            
            int region = (hit.getLayer() - 1) / 12 + 1;
            
            double r = Math.sqrt(t.getVertex().x() * t.getVertex().x() + t.getVertex().y() * t.getVertex().y());
            double weight = energyWeight ? t.getKinEnergy() : 1;
//            System.out.println(weight);
            group.getH2F("hi_bg_origin_rz_region" + region).fill(t.getVertex().z(), r, weight);
            if (t.getVertex().z() > 1000 * (region)) 
                group.getH2F("hi_bg_origin_xy_region" + region).fill(t.getVertex().x(), t.getVertex().y(), weight);
            
            double value = energyWeight ? t.getKinEnergy() : t.getVertex().z();
            group.getH1F("hi_bg_region" + region + "_all").fill(value);
            if(this.pidToName(Math.abs(Math.abs(t.getPid())))!=null) 
                group.getH1F("hi_bg_region" + region + "_" + this.pidToName(Math.abs(t.getPid()))).fill(value);
            else 
                group.getH1F("hi_bg_region" + region + "_other").fill(value);
        }

    }

    public void fillSectorBG(DataGroup group, List<Hit> hits) {

        for (Hit hit : hits) {
            
            True t = hit.getTrue();
            if(t==null) continue;

            int region = (hit.getLayer() - 1) / 12 + 1;
            int sector = hit.getSector();
            
            if(region!=1) continue;
            
             group.getH1F("hi_bg_r1_s" + sector + "_all").fill(t.getVertex().z());
            if(this.pidToName(Math.abs(Math.abs(t.getPid())))!=null) 
                group.getH1F("hi_bg_r1_s" + sector + "_" + this.pidToName(Math.abs(t.getPid()))).fill(t.getVertex().z());
            else 
                group.getH1F("hi_bg_r1_s" + sector + "_other").fill(t.getVertex().z());
        }

    }

    public void fillPosition(DataGroup group, List<Hit> hits) {

        for (Hit hit : hits) {

            True t = hit.getTrue();
            if(t==null) continue;

            int region = (hit.getLayer() - 1) / 12 + 1;
 
            double r = Math.sqrt(t.getPosition().y() * t.getPosition().y() + t.getPosition().x() * t.getPosition().x());
//            if (t.getVertex().z() > -150 && t.getVertex().z() < 100 && (t.getPid() == 11 || t.getPid() == -11)) {
                group.getH2F("hi-posZ-posR-region" + region).fill(t.getPosition().z(), r);
                group.getH2F("hi-posX-posY-region" + region).fill(t.getPosition().x(), t.getPosition().y());
//            }
        }

    }

    @Override
    public void analyzeHistos() {
        this.normalizeToEventsX100(this.getHistos().get("Sector Occupancy"));
        double norm = 112 * 12 / 100;
        this.normalizeToEvents(this.getHistos().get("Region Occupancy"));
        this.normalize(this.getHistos().get("Region Occupancy"), norm);

        this.fitDataGroup(this.getHistos().get("Region Occupancy"));

        this.divide(this.getHistos().get("Origin of Bg - Energy"), this.getHistos().get("Origin of Bg"));
        this.normalizeToTime(this.getHistos().get("Origin of Bg"));
        this.normalizeToTime(this.getHistos().get("Origin of Bg - Energy"));

    }


    @Override
    public void fitDataGroup(DataGroup dg) {
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            F1D f = fitPol0(dg.getH1F("hi_occ_region" + region));
            dg.addDataSet(f, 0);
        }
    }


    @Override
    public void setPlottingOptions(String name) {
        if(name.equals("Region Occupancy")) {
            if(this.getHistos().get("Sector Occupancy").getH2F("hi_occ_sector1c")!=null) {
                TCanvas canvas = new TCanvas("DC", 1000, 800);
                canvas.divide(2, 2);
                canvas.getCanvas().setGridX(false);
                canvas.getCanvas().setGridY(false);
                canvas.cd(1);
                canvas.draw(this.getHistos().get("Sector Occupancy").getH2F("hi_occ_sector1"));
                canvas.getCanvas().getPad().getAxisZ().setRange(0, 6);
                canvas.cd(2);
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region1"));
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region2"),"same");
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region3"),"same");
                canvas.getCanvas().getPad().getAxisZ().setRange(0, 7);
                canvas.cd(3);
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_layer"));
                canvas.cd(0);
                canvas.draw(this.getHistos().get("Sector Occupancy").getH2F("hi_occ_sector1c"));
                canvas.getCanvas().getPad().getAxisZ().setRange(0, 7);            
                canvas.cd(2);
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region1c"),"same");
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region2c"),"same");
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region3c"),"same");
                canvas.getCanvas().getPad().getAxisZ().setRange(0, 7);
                canvas.cd(3);
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_layerc"),"same");
                this.getHistos().get("Region Occupancy").getH1F("hi_layer").setLineWidth(1);
                this.getHistos().get("Region Occupancy").getH1F("hi_wire_region1").setLineWidth(1);
                this.getHistos().get("Region Occupancy").getH1F("hi_wire_region2").setLineWidth(1);
                this.getHistos().get("Region Occupancy").getH1F("hi_wire_region3").setLineWidth(1);
            }
            else {
                TCanvas canvas = new TCanvas("DC", 1300, 600);
                canvas.divide(2, 1);
                canvas.getCanvas().setGridX(false);
                canvas.getCanvas().setGridY(false);
                canvas.cd(0);
                canvas.draw(this.getHistos().get("Sector Occupancy").getH2F("hi_occ_sector1"));
                canvas.cd(1);
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region1"));
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region2"),"same");
                canvas.draw(this.getHistos().get("Region Occupancy").getH1F("hi_wire_region3"),"same");
            }
        } 
        for(EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
            pad.setTitle("");
        }
        if(name.equals("Region Occupancy")) {
            DataGroup dg = this.getHistos().get(name);
            for (int ir = 0; ir<NREGIONS; ir++) {
                List<IDataSet> ds = dg.getData(0);
                for (IDataSet d : ds) {
                    if (d instanceof F1D && d.getName().contains(""+(ir+1))) {
                        double par = ((F1D) d).getParameter(0);
                        String text = String.format("\tRegion %d: %.3f", ir+1, par)+"%";
                        LatexText latexText = new LatexText(text, 70, (ir+1)*30);
                        latexText.setColor(ir+2);
                        latexText.setFontSize(30);
                        latexText.setFont("Arial");
                        this.getCanvas().getCanvas(name).cd(0);
                        this.getCanvas().getCanvas(name).draw(latexText);
                    }
//                    else if(d instanceof H1F && d.getName().endsWith("c"))
//                        ((H1F) d).getAttributes().setLineStyle(2);
                }
            }
        }
        else if(name.contains("Origin of Bg"))
            if(name.contains("Sector"))
                this.setLegend(name, 450, 150);
            else
                this.setLegend(name, 250, 140);
        if(name.contains("Energy")) {
            for(EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
               if(pad.getDatasetPlotters().get(0).getDataSet() instanceof H1F)
                       pad.getAxisY().setLog(true);
        }
        if(!name.contains("Occupancy"))
            this.setLogZ(name);
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
