package modules;

import analysis.Constants;
import java.util.List;
import objects.Hit;
import analysis.Module;
import java.util.ArrayList;
import objects.Event;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.LatexText;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.graphics.GraphicsAxis;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author devita
 */
public class DCmodule extends Module {

    private static final int NREGIONS = 3;
    private static final int NSECTORS = 6;
    private static final int NLAYERS = 36;
    private static final int NWIRES = 112;

    // private static final double[] RWINDOWS = {500, 1400, 1200};
    private static final double[] RWINDOWS = { 250, 500, 500 };
    private static final double[] DR = { 2500, 4000, 5500 };
    private static final double[] DZ = { 3500, 5000, 6500 };

    /////////////// To fix max ///////
    private static final double OCC_Z_MAX = 20; // tune once
    private static final double BG_Y_MAX = 2000; // tune once
    // private static final double BG_SEC_Y_MAX = 50000; // tune once
    private static final double POS_Z_MAX = 2000; // for position plots
    private static final double Region_OCC_Max = 30;
    private static final double FLAME_VX_CUT = 0; // tune this
    private static final double FLAME_VY_CUT = 100; // tune this

    private static final double FLAME_VX_LEFT_CUT = -1200; // x left cut
    private static final double FLAME_VX_RIGHT_CUT = 300; // x right cut

    private static final double FLAME_POS_X_LEFT = -5000; // mm
    private static final double FLAME_POS_X_RIGHT = 5000; // mm
    private static final double FLAME_POS_Y_ABS = 60; // mm

    private boolean passVertexFlameCut(True t) {
        return (t.getVertex().x() > FLAME_VX_LEFT_CUT &&
                t.getVertex().x() < 0 &&
                Math.abs(t.getVertex().y()) < FLAME_VY_CUT);
    }

    private boolean passPositionFlameCut(True t) {
        return (t.getPosition().x() > FLAME_POS_X_LEFT &&
                t.getPosition().x() < FLAME_POS_X_RIGHT &&
                Math.abs(t.getPosition().y()) < FLAME_POS_Y_ABS);
    }

    private boolean passPositionFlameAntiCut(True t) {
        return !passPositionFlameCut(t);
    }

    private static final int CUT_MODE_NONE = 0;
    private static final int CUT_MODE_VERTEX = 1;
    private static final int CUT_MODE_POSITION = 2;
    private static final int CUT_MODE_POSITION_ANTI = 3;

    // private static final int CUT_MODE = CUT_MODE_POSITION;
    // private static final int CUT_MODE = CUT_MODE_POSITION_ANTI;
    private static final int CUT_MODE = CUT_MODE_NONE;

    // Raw hit counts per (region, sector) for Poisson error bars on Region
    // Occupancy
    private final long[][] rawDC = new long[NREGIONS][NSECTORS];
    // Region-averaged occupancy [%] saved before DataGroup swap; used in
    // setPlottingOptions
    private final double[] regionAvg = new double[NREGIONS];

    // Species-split Region Occupancy: electron vs photon/neutral (gamma+neutron)
    private static final String[] OCC_SPECIES = { "electron", "photonNeutral" };
    private static final String[] OCC_SPECIES_LABELS = { "Electron", "Photon/Neutral" };
    private static final int[] OCC_SPECIES_COLOR = { 2, 4 };
    // [species][region][sector] raw hit counts for Poisson error bars
    private final long[][][] rawDCspecies = new long[OCC_SPECIES.length][NREGIONS][NSECTORS];

    private String occSpeciesName(True t) {
        if (t == null)
            return null;
        int pid = Math.abs(t.getPid());
        if (pid == 11)
            return "electron";
        if (pid == 22 || pid == 2112 || pid == 111)
            return "photonNeutral";
        return null;
    }

    // [species][region][sector] raw hit counts for the mpid-aware split below
    private final long[][][] rawDCspeciesOrigin = new long[OCC_SPECIES.length][NREGIONS][NSECTORS];

    /**
     * Same electron vs photon/neutral split as occSpeciesName(), but traces a
     * pid==11 hit back to its mother (mpid): a "electron" hit whose mother is a
     * photon or neutron (e.g. a Compton/photoelectric/pair-production electron,
     * or a neutron-induced secondary) is reclassified as photonNeutral, since the
     * actual background source is the neutral, not the electron it kicked out.
     * Mirrors the mpid fallback HTCCmodule.particleName() uses for optical photons.
     */
    private String occSpeciesNameOrigin(True t) {
        if (t == null)
            return null;
        int pid = Math.abs(t.getPid());
        if (pid == 22 || pid == 2112 || pid == 111)
            return "photonNeutral";
        if (pid == 11) {
            int mpid = Math.abs(t.getMPID());
            return (mpid == 22 || mpid == 2112 || mpid == 111) ? "photonNeutral" : "electron";
        }
        return null;
    }

    private boolean passSelectedCut(True t) {
        if (CUT_MODE == CUT_MODE_NONE)
            return true;
        if (CUT_MODE == CUT_MODE_VERTEX)
            return passVertexFlameCut(t);
        if (CUT_MODE == CUT_MODE_POSITION)
            return passPositionFlameCut(t);
        if (CUT_MODE == CUT_MODE_POSITION_ANTI)
            return passPositionFlameAntiCut(t);
        return true;
    }

    public DCmodule() {
        super(DetectorType.DC);
    }

    /*
     * occupacy layer vs wire for each sector - need to normalize wrt number of
     * region (3)
     */
    public DataGroup occupancies() {
        DataGroup dg = new DataGroup(3, 2);
        for (int is = 0; is < NSECTORS; is++) {
            int sector = is + 1;
            String name = "sector" + sector;
            H2F hi_occ = histo2D("hi_occ_" + name, "Wire", "layer", NWIRES, 1, NWIRES + 1, NLAYERS, 1, NLAYERS + 1);
            dg.addDataSet(hi_occ, 0 + is);
        }
        return dg;
    }

    /*
     * occupancy vs sector for each region - need to normalize wrt number of layer
     * (12) and wires (112)
     */
    public DataGroup occupancy_region() {
        DataGroup dg = new DataGroup(1, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            String name = "region" + region;
            H1F hi_occ = histo1D("hi_occ_" + name, name, "Sector", "Occupancy[%] ", NSECTORS, 0.5, NSECTORS + 0.5, 0);
            hi_occ.setLineColor(region + 1);
            hi_occ.setLineWidth(4);

            dg.addDataSet(hi_occ, 0);

        }
        return dg;
    }

    /*
     * occupancy vs sector for each region, split by species (electron vs
     * photon/neutral) - one pad per region
     */
    public DataGroup occupancy_region_species() {
        DataGroup dg = new DataGroup(NREGIONS, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            for (int isp = 0; isp < OCC_SPECIES.length; isp++) {
                H1F hi_occ = histo1D("hi_occ_species_region" + region + "_" + OCC_SPECIES[isp],
                        OCC_SPECIES_LABELS[isp], "Sector", "Occupancy [%] ",
                        NSECTORS, 0.5, NSECTORS + 0.5, 0);
                hi_occ.setLineColor(OCC_SPECIES_COLOR[isp]);
                hi_occ.setLineWidth(4);
                dg.addDataSet(hi_occ, ir);
            }
        }
        return dg;
    }

    /*
     * same as occupancy_region_species(), but using the mpid-aware
     * occSpeciesNameOrigin() classification - one pad per region
     */
    public DataGroup occupancy_region_species_origin() {
        DataGroup dg = new DataGroup(NREGIONS, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            for (int isp = 0; isp < OCC_SPECIES.length; isp++) {
                H1F hi_occ = histo1D("hi_occ_species_origin_region" + region + "_" + OCC_SPECIES[isp],
                        OCC_SPECIES_LABELS[isp], "Sector", "Occupancy [%] ",
                        NSECTORS, 0.5, NSECTORS + 0.5, 0);
                hi_occ.setLineColor(OCC_SPECIES_COLOR[isp]);
                hi_occ.setLineWidth(4);
                dg.addDataSet(hi_occ, ir);
            }
        }
        return dg;
    }

    /* histos to understand origin of BG */
    public DataGroup[] origin_bg() {
        DataGroup[] dg = new DataGroup[2];

        for (int i = 0; i < dg.length; i++) {

            dg[i] = new DataGroup(3, 3);

            for (int ir = 0; ir < NREGIONS; ir++) {
                int region = ir + 1;
                H2F hi_bg_origin_rz = histo2D("hi_bg_origin_rz_region" + region, "Vz(m)", "r(m) ", 200, -0.5,
                        DZ[ir] / 1000,
                        200, -0.50, DR[ir] / 1000);
                H2F hi_bg_origin_xy = histo2D("hi_bg_origin_xy_region" + region, "Vx(m)", "Vy(m) ", 200, -DR[ir] / 1000,
                        DR[ir] / 1000, 200, -DR[ir] / 1000, DR[ir] / 1000);
                dg[i].addDataSet(hi_bg_origin_xy, 0 + ir);
                dg[i].addDataSet(hi_bg_origin_rz, 3 + ir);

                double min = -0.5;
                double max = DZ[ir] / 1000;
                String name = "Vz(m)";
                if (i == 1) {
                    min = 0;
                    max = 200;
                    name = "E(MeV)";
                }
                for (int ip = 0; ip < PNAMES.length; ip++) {
                    H1F hi_bg = histo1D("hi_bg_region" + region + "_" + PNAMES[ip], PNAMES[ip], name, "Rate [MHz] ",
                            200, -1.0, max, 0);
                    this.setHistoAttr(hi_bg, ip < 5 ? ip + 1 : ip + 3);
                    dg[i].addDataSet(hi_bg, 6 + ir);
                }
            }

        }

        return dg;
    }

    public DataGroup sectorBG(int region) {
        DataGroup dg = new DataGroup(2, 3);
        for (int is = 0; is < NSECTORS; is++) {
            int sector = is + 1;
            for (int ip = 0; ip < PNAMES.length; ip++) {
                H1F hi_bg = histo1D(
                        "hi_bg_r" + region + "_s" + sector + "_" + PNAMES[ip],
                        "R" + region + "S" + sector + "-" + PNAMES[ip],
                        "Vz(m)", "Rate [MHz]",
                        200, -1.0, DZ[region - 1] / 1000, 0);
                this.setHistoAttr(hi_bg, ip < 5 ? ip + 1 : ip + 3);
                dg.addDataSet(hi_bg, is);
            }
        }
        return dg;
    }
    // public DataGroup sectorBG() {
    // DataGroup dg = new DataGroup(2,3);

    // for (int is = 0; is < NSECTORS; is++) {
    // int sector = is + 1;
    // for (int ip=0; ip<PNAMES.length; ip++) {
    // H1F hi_bg = histo1D("hi_bg_r1_s" + sector + "_" + PNAMES[ip], "R1S" +sector +
    // "-" + PNAMES[ip], "Vz(mm)", "Rate [MHz] ", 200, -1000, DZ[0], 0);
    // this.setHistoAttr(hi_bg, ip<5 ? ip+1 : ip+3);
    // dg.addDataSet(hi_bg, is);
    // }
    // }
    // return dg;
    // }

    public DataGroup pos_bg() {
        DataGroup dg = new DataGroup(3, 2);

        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            String name = "region" + region;
            H2F hi_posZ_posYR = histo2D("hi-posZ-posR-" + name, "posZ [m]", "posR [m]", 300, 0., DZ[ir] / 1000, 300, 0,
                    DR[ir] / 1000);
            H2F hi_posX_posY = histo2D("hi-posX-posY-" + name, "posX [m]", "posY [m]", 200, -DR[ir] / 1000,
                    DR[ir] / 1000, 200,
                    -DR[ir] / 1000, DR[ir] / 1000);
            dg.addDataSet(hi_posZ_posYR, ir);
            dg.addDataSet(hi_posX_posY, ir + 3);
        }

        return dg;
    }

    public void fixOccupancyZAxis(DataGroup dg) {
        double globalMax = 0;

        // Find global max across all sectors
        for (int i = 0; i < 6; i++) {
            H2F h = dg.getH2F("hi_occ_sector" + (i + 1));
            if (h.getMaximum() > globalMax) {
                globalMax = h.getMaximum();
            }
        }

        // Apply same Z scale to all
        for (int i = 0; i < 6; i++) {
            H2F h = dg.getH2F("hi_occ_sector" + (i + 1));
            // h.setMaximum(globalMax);
            // h.setMinimum(0); // optional but recommended
        }
    }

    @Override
    public void createHistos() {
        this.getHistos().put("Sector Occupancy", this.occupancies());
        this.getHistos().put("Region Occupancy", this.occupancy_region());
        this.getHistos().put("Region Occupancy - Species", this.occupancy_region_species());
        this.getHistos().put("Region Occupancy - Species (Origin)", this.occupancy_region_species_origin());
        this.getHistos().put("Origin of Bg", this.origin_bg()[0]);
        this.getHistos().put("Origin of Bg - Energy", this.origin_bg()[1]);
        this.getHistos().put("Origin of Bg - Sector R1", this.sectorBG(1));
        this.getHistos().put("Origin of Bg - Sector R2", this.sectorBG(2));
        this.getHistos().put("Origin of Bg - Sector R3", this.sectorBG(3));
        this.getHistos().put("Position of Bg", this.pos_bg());
    }

    @Override
    public void fillHistos(Event event) {
        List<Hit> allhits = event.getHits(DetectorType.DC);
        if (allhits != null) {
            List<Hit> hits = new ArrayList<>();

            for (Hit h : allhits) {
                int region = (h.getLayer() - 1) / 12 + 1;

                // Exclude first 7 wires of Region 1, Sector 4
                boolean isR1S4HotWire = (region == 1
                        && h.getSector() == 4
                        && h.getComponent() <= 7); // change to <= 8 if needed

                boolean isR2S1Hot = (region == 2
                        && h.getSector() == 1
                        && h.getTrue() != null && h.getTrue().getPosition().x() > 2350);

                boolean isR2S1HotWire_upper = (region == 2
                        && h.getSector() == 1
                        && h.getLayer() >= 18
                        && h.getComponent() >= 95); // layers 18-24, last 12 wires

                boolean isR2S1HotWire_lower = (region == 2
                        && h.getSector() == 1
                        && h.getLayer() < 18
                        && h.getComponent() >= 101); // layers 13-17, adjust wire threshold

                // if (h.getTrue().getEdep() > 50E-6 && !isR1S4HotWire && !isR2S1HotWire_upper
                // && !isR2S1HotWire_lower)
                if (h.getTrue() == null)
                    continue;
                if (h.getTrue().getEdep() > 50E-6)
                    hits.add(h);
            }

            this.fillOccupancies(this.getHistos().get("Sector Occupancy"), hits);
            this.fillOccupancy_region(this.getHistos().get("Region Occupancy"), hits);
            this.fillOccupancyRegionSpecies(this.getHistos().get("Region Occupancy - Species"), hits);
            this.fillOccupancyRegionSpeciesOrigin(this.getHistos().get("Region Occupancy - Species (Origin)"), hits);
            this.fillOrigin(this.getHistos().get("Origin of Bg"), hits, false);
            this.fillOrigin(this.getHistos().get("Origin of Bg - Energy"), hits, true);
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector R1"), hits, 1);
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector R2"), hits, 2);
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector R3"), hits, 3);
            this.fillPosition(this.getHistos().get("Position of Bg"), hits);
        }
    }

    public void fillOccupancies(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {

            True t = hit.getTrue(); ///////// kr add
            if (!passSelectedCut(t))
                continue;
            // System.out.println(hit.getTrue().getEdep() + " " + hit.getTrue().getTime()+ "
            // " + hit.getTDC());
            group.getH2F("hi_occ_sector" + hit.getSector()).fill(hit.getComponent(), hit.getLayer(),
                    RWINDOWS[(hit.getLayer() - 1) / 12] / Constants.getTimeWindow());
        }
    }

    public void fillOccupancy_region(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            True t = hit.getTrue(); ////// kr add
            if (!passSelectedCut(t))
                continue;
            int region = (hit.getLayer() - 1) / 12 + 1;
            group.getH1F("hi_occ_region" + region).fill(hit.getSector(),
                    RWINDOWS[region - 1] / Constants.getTimeWindow());
            rawDC[region - 1][hit.getSector() - 1]++;
        }
    }

    public void fillOccupancyRegionSpecies(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            True t = hit.getTrue();
            if (!passSelectedCut(t))
                continue;
            String species = occSpeciesName(t);
            if (species == null)
                continue;
            int region = (hit.getLayer() - 1) / 12 + 1;
            int sector = hit.getSector();
            group.getH1F("hi_occ_species_region" + region + "_" + species).fill(sector,
                    RWINDOWS[region - 1] / Constants.getTimeWindow());
            int speciesIdx = species.equals("electron") ? 0 : 1;
            rawDCspecies[speciesIdx][region - 1][sector - 1]++;
        }
    }

    public void fillOccupancyRegionSpeciesOrigin(DataGroup group, List<Hit> hits) {
        for (Hit hit : hits) {
            True t = hit.getTrue();
            if (!passSelectedCut(t))
                continue;
            String species = occSpeciesNameOrigin(t);
            if (species == null)
                continue;
            int region = (hit.getLayer() - 1) / 12 + 1;
            int sector = hit.getSector();
            group.getH1F("hi_occ_species_origin_region" + region + "_" + species).fill(sector,
                    RWINDOWS[region - 1] / Constants.getTimeWindow());
            int speciesIdx = species.equals("electron") ? 0 : 1;
            rawDCspeciesOrigin[speciesIdx][region - 1][sector - 1]++;
        }
    }

    public void fillOrigin(DataGroup group, List<Hit> hits, boolean energyWeight) {

        for (Hit hit : hits) {

            True t = hit.getTrue();
            // ===== ADD THIS CUT for just sheet of flame =====
            if (!passSelectedCut(t))
                continue;
            int region = (hit.getLayer() - 1) / 12 + 1;

            double vx = t.getVertex().x() / 1000.0;
            double vy = t.getVertex().y() / 1000.0;
            double vz = t.getVertex().z() / 1000.0;
            double r = Math.sqrt(vx * vx + vy * vy);
            double weight = energyWeight ? t.getKinEnergy() : 1;
            // System.out.println(weight);
            group.getH2F("hi_bg_origin_rz_region" + region).fill(vz, r, weight);
            if (vz > region)
                group.getH2F("hi_bg_origin_xy_region" + region).fill(vx, vy, weight);

            double value = energyWeight ? t.getKinEnergy() : vz;
            group.getH1F("hi_bg_region" + region + "_all").fill(value);
            if (this.pidToName(Math.abs(Math.abs(t.getPid()))) != null)
                group.getH1F("hi_bg_region" + region + "_" + this.pidToName(Math.abs(t.getPid()))).fill(value);
            else
                group.getH1F("hi_bg_region" + region + "_other").fill(value);
        }

    }

    public void fillSectorBG(DataGroup group, List<Hit> hits, int targetRegion) {
        for (Hit hit : hits) {
            True t = hit.getTrue();
            if (!passSelectedCut(t))
                continue;

            int region = (hit.getLayer() - 1) / 12 + 1;
            int sector = hit.getSector();

            if (region != targetRegion)
                continue; // ← only fill matching region

            double vz = t.getVertex().z() / 1000.0;
            group.getH1F("hi_bg_r" + region + "_s" + sector + "_all")
                    .fill(vz);
            if (this.pidToName(Math.abs(t.getPid())) != null)
                group.getH1F("hi_bg_r" + region + "_s" + sector + "_"
                        + this.pidToName(Math.abs(t.getPid())))
                        .fill(vz);
            else
                group.getH1F("hi_bg_r" + region + "_s" + sector + "_other")
                        .fill(vz);
        }
    }

    public void fillPosition(DataGroup group, List<Hit> hits) {

        for (Hit hit : hits) {

            True t = hit.getTrue();
            if (!passSelectedCut(t))
                continue;
            int region = (hit.getLayer() - 1) / 12 + 1;

            double px = t.getPosition().x() / 1000.0;
            double py = t.getPosition().y() / 1000.0;
            double pz = t.getPosition().z() / 1000.0;
            double r = Math.sqrt(py * py + px * px);
            // if (t.getVertex().z() > -150 && t.getVertex().z() < 100 && (t.getPid() == 11
            // || t.getPid() == -11)) {
            group.getH2F("hi-posZ-posR-region" + region).fill(pz, r);
            group.getH2F("hi-posX-posY-region" + region).fill(px, py);
            // }
        }

    }

    @Override
    public void analyzeHistos() {
        //////// luminosity scale
        // Parse -lumi flag from command line args
        double x = Double.parseDouble(System.getProperty("lumi", "117000.0")); // default: no scaling
        double lumiScale = x / 117000.0;

        this.normalizeToEventsX100(this.getHistos().get("Sector Occupancy"));
        this.normalize(this.getHistos().get("Sector Occupancy"), lumiScale);

        fixOccupancyZAxis(this.getHistos().get("Sector Occupancy")); // to fix the Z axis scale for all sectors to be
                                                                     // the same and comparable
        // Rescale RGA to RGH effective luminosity
        ////// this.normalize(this.getHistos().get("Sector Occupancy"), 1.0 /
        // RGA_TO_RGH_LUMI_SCALE);

        double norm = 112 * 12 / 100;
        this.normalizeToEvents(this.getHistos().get("Region Occupancy"));
        this.normalize(this.getHistos().get("Region Occupancy"), norm);
        this.normalize(this.getHistos().get("Region Occupancy"), lumiScale);

        this.fitDataGroup(this.getHistos().get("Region Occupancy"));

        // Set Poisson errors on each sector bin and save per-region averages
        DataGroup regGroup = this.getHistos().get("Region Occupancy");
        for (int ir = 0; ir < NREGIONS; ir++) {
            H1F h = regGroup.getH1F("hi_occ_region" + (ir + 1));
            double avg = 0;
            for (int is = 0; is < NSECTORS; is++) {
                int bin = h.getAxis().getBin(is + 1.0);
                long raw = rawDC[ir][is];
                if (raw > 0)
                    h.setBinError(bin, h.getBinContent(bin) / Math.sqrt(raw));
                avg += h.getBinContent(bin);
            }
            regionAvg[ir] = avg / NSECTORS;
        }

        // Region Occupancy split by species (electron vs photon/neutral), same
        // normalization + Poisson-error treatment as the "all" Region Occupancy above
        this.normalizeToEvents(this.getHistos().get("Region Occupancy - Species"));
        this.normalize(this.getHistos().get("Region Occupancy - Species"), norm);
        this.normalize(this.getHistos().get("Region Occupancy - Species"), lumiScale);

        DataGroup speciesGroup = this.getHistos().get("Region Occupancy - Species");
        for (int isp = 0; isp < OCC_SPECIES.length; isp++) {
            for (int ir = 0; ir < NREGIONS; ir++) {
                H1F h = speciesGroup.getH1F("hi_occ_species_region" + (ir + 1) + "_" + OCC_SPECIES[isp]);
                for (int is = 0; is < NSECTORS; is++) {
                    int bin = h.getAxis().getBin(is + 1.0);
                    long raw = rawDCspecies[isp][ir][is];
                    if (raw > 0)
                        h.setBinError(bin, h.getBinContent(bin) / Math.sqrt(raw));
                }
            }
        }

        // Same split, but mpid-aware (traces pid==11 hits back to a photon/neutron
        // mother) - tests whether "electron" occupancy is actually neutral-induced
        this.normalizeToEvents(this.getHistos().get("Region Occupancy - Species (Origin)"));
        this.normalize(this.getHistos().get("Region Occupancy - Species (Origin)"), norm);
        this.normalize(this.getHistos().get("Region Occupancy - Species (Origin)"), lumiScale);

        DataGroup speciesGroupOrigin = this.getHistos().get("Region Occupancy - Species (Origin)");
        for (int isp = 0; isp < OCC_SPECIES.length; isp++) {
            for (int ir = 0; ir < NREGIONS; ir++) {
                H1F h = speciesGroupOrigin.getH1F("hi_occ_species_origin_region" + (ir + 1) + "_" + OCC_SPECIES[isp]);
                for (int is = 0; is < NSECTORS; is++) {
                    int bin = h.getAxis().getBin(is + 1.0);
                    long raw = rawDCspeciesOrigin[isp][ir][is];
                    if (raw > 0)
                        h.setBinError(bin, h.getBinContent(bin) / Math.sqrt(raw));
                }
            }
        }

        this.divide(this.getHistos().get("Origin of Bg - Energy"), this.getHistos().get("Origin of Bg"));
        this.normalizeToTime(this.getHistos().get("Origin of Bg"));
        this.normalize(this.getHistos().get("Origin of Bg"), lumiScale);

        this.normalizeToTime(this.getHistos().get("Origin of Bg - Energy"));

        // ← ADD HERE — normalize all three sector BG groups
        for (int ir = 1; ir <= 3; ir++) {
            this.normalizeToTime(this.getHistos().get("Origin of Bg - Sector R" + ir));
            this.normalize(this.getHistos().get("Origin of Bg - Sector R" + ir), lumiScale);
        }

        // // FIX: Origin of Bg - Sector (same treatment as Origin of Bg)
        // this.normalizeToTime(this.getHistos().get("Origin of Bg - Sector"));
        // this.normalize(this.getHistos().get("Origin of Bg - Sector"), lumiScale);

        // FIX: Position of Bg (if you want it comparable across runs)
        this.normalizeToTime(this.getHistos().get("Position of Bg"));
        this.normalize(this.getHistos().get("Position of Bg"), lumiScale);

        // Replace H1F DataGroup with GraphErrors + F1D average lines so that
        // drawHistos() renders PE-style markers instead of bar histograms
        DataGroup regGroupFinal = this.getHistos().get("Region Occupancy");
        int[] rColors = { 2, 3, 4 };
        DataGroup regGroupPE = new DataGroup(1, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            H1F h = regGroupFinal.getH1F("hi_occ_region" + (ir + 1));
            GraphErrors gr = new GraphErrors("Region " + (ir + 1));
            gr.setTitle("Region " + (ir + 1));
            gr.setTitleX("Sector");
            gr.setTitleY("Occupancy [%]");
            gr.setMarkerColor(rColors[ir]);
            gr.setMarkerStyle(20);
            gr.setMarkerSize(8);
            gr.setLineColor(rColors[ir]);
            for (int ib = 0; ib < h.getAxis().getNBins(); ib++) {
                double xc = h.getAxis().getBinCenter(ib);
                double ex = h.getAxis().getBinWidth(ib) / 2.0;
                double y = h.getBinContent(ib);
                if (y > 0)
                    gr.addPoint(xc, y, ex, h.getBinError(ib));
            }
            regGroupPE.addDataSet(gr, 0);
        }
        // preserve the F1D average lines from fitDataGroup
        for (IDataSet ds : regGroupFinal.getData(0)) {
            if (ds instanceof F1D)
                regGroupPE.addDataSet(ds, 0);
        }
        this.getHistos().put("Region Occupancy", regGroupPE);

        // Same H1F -> GraphErrors (PE-style) swap for the species-split plot,
        // one pad per region, electron vs photon/neutral markers
        DataGroup speciesGroupFinal = this.getHistos().get("Region Occupancy - Species");
        DataGroup speciesGroupPE = new DataGroup(NREGIONS, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            for (int isp = 0; isp < OCC_SPECIES.length; isp++) {
                H1F h = speciesGroupFinal.getH1F("hi_occ_species_region" + (ir + 1) + "_" + OCC_SPECIES[isp]);
                GraphErrors gr = new GraphErrors(OCC_SPECIES_LABELS[isp]);
                gr.setTitle(OCC_SPECIES_LABELS[isp]);
                gr.setTitleX("Sector");
                gr.setTitleY("Occupancy [%]");
                gr.setMarkerColor(OCC_SPECIES_COLOR[isp]);
                gr.setMarkerStyle(20);
                gr.setMarkerSize(8);
                gr.setLineColor(OCC_SPECIES_COLOR[isp]);
                for (int ib = 0; ib < h.getAxis().getNBins(); ib++) {
                    double xc = h.getAxis().getBinCenter(ib);
                    double ex = h.getAxis().getBinWidth(ib) / 2.0;
                    double y = h.getBinContent(ib);
                    if (y > 0)
                        gr.addPoint(xc, y, ex, h.getBinError(ib));
                }
                speciesGroupPE.addDataSet(gr, ir);
            }
        }
        this.getHistos().put("Region Occupancy - Species", speciesGroupPE);

        // Same swap for the mpid-aware (Origin) species split
        DataGroup speciesGroupOriginFinal = this.getHistos().get("Region Occupancy - Species (Origin)");
        DataGroup speciesGroupOriginPE = new DataGroup(NREGIONS, 1);
        for (int ir = 0; ir < NREGIONS; ir++) {
            for (int isp = 0; isp < OCC_SPECIES.length; isp++) {
                H1F h = speciesGroupOriginFinal
                        .getH1F("hi_occ_species_origin_region" + (ir + 1) + "_" + OCC_SPECIES[isp]);
                GraphErrors gr = new GraphErrors(OCC_SPECIES_LABELS[isp]);
                gr.setTitle(OCC_SPECIES_LABELS[isp]);
                gr.setTitleX("Sector");
                gr.setTitleY("Occupancy [%]");
                gr.setMarkerColor(OCC_SPECIES_COLOR[isp]);
                gr.setMarkerStyle(20);
                gr.setMarkerSize(8);
                gr.setLineColor(OCC_SPECIES_COLOR[isp]);
                for (int ib = 0; ib < h.getAxis().getNBins(); ib++) {
                    double xc = h.getAxis().getBinCenter(ib);
                    double ex = h.getAxis().getBinWidth(ib) / 2.0;
                    double y = h.getBinContent(ib);
                    if (y > 0)
                        gr.addPoint(xc, y, ex, h.getBinError(ib));
                }
                speciesGroupOriginPE.addDataSet(gr, ir);
            }
        }
        this.getHistos().put("Region Occupancy - Species (Origin)", speciesGroupOriginPE);
    }

    @Override
    public void fitDataGroup(DataGroup dg) {
        for (int ir = 0; ir < NREGIONS; ir++) {
            int region = ir + 1;
            // F1D f = fitPol0(dg.getH1F("hi_occ_region" + region));
            // dg.addDataSet(f, 0);
            H1F h = dg.getH1F("hi_occ_region" + region);

            double avg = 0.0;
            double sum = 0.0;

            for (int bin = 0; bin < h.getDataSize(0); bin++) {
                avg += h.getBinContent(bin);
                sum += h.getBinContent(bin);
            }

            avg /= h.getDataSize(0);

            // Create constant line at average value
            F1D f = new F1D("avg_region" + region, "[a]", 0.5, 6.5);
            // String text = String.format("Region%d Avg: %.1f%% Sum: %.1f%%", ir + 1, avg,
            // sum);
            String text = String.format("Region%d Avg: %.1f%%", ir + 1, avg);

            f.setParameter(0, avg);

            f.setLineColor(region + 1);
            f.setLineWidth(3);
            f.setLineStyle(2);
            dg.addDataSet(f, 0);
        }
    }

    @Override
    public void setPlottingOptions(String name) {

        for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
            pad.setTitle("");
        }

        if (name.equals("Sector Occupancy")) {
            GStyle.getAxisAttributesX().setTitleFontSize(20);
            GStyle.getAxisAttributesY().setTitleFontSize(20);
            GStyle.getAxisAttributesX().setLabelFontSize(16);
            GStyle.getAxisAttributesY().setLabelFontSize(16);
            // for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
            // pad.getAxisZ().setRange(1, OCC_Z_MAX);
            // }
            this.setupColorAxis(name, "Occupancy [%]");
        }

        else if (name.equals("Region Occupancy")) {
            GStyle.getAxisAttributesX().setTitleFontSize(36);
            GStyle.getAxisAttributesY().setTitleFontSize(36);
            GStyle.getAxisAttributesX().setLabelFontSize(32);
            GStyle.getAxisAttributesY().setLabelFontSize(32);
            for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
                // pad.getAxisY().setRange(0, Region_OCC_Max);
                pad.getAxisX().setRange(0.5, 6.5);
            }
            // regionAvg[] was computed in analyzeHistos() before the DataGroup swap.
            // LatexText y-coordinate is in canvas pixels (y=0 = top, y=900 = bottom).
            // Using 800 - ir*30 places labels near the bottom of the plot area (~y=1-3%).
            for (int ir = 0; ir < NREGIONS; ir++) {
                String text = String.format("Region%d Avg: %.1f%%", ir + 1, regionAvg[ir]);
                LatexText latexText = new LatexText(text, 100, 60 + (ir) * 30);
                latexText.setColor(ir + 2);
                latexText.setFontSize(32);
                latexText.setFont("Arial");
                this.getCanvas().getCanvas(name).draw(latexText);
            }
        }

        else if (name.equals("Region Occupancy - Species")) {
            GStyle.getAxisAttributesX().setTitleFontSize(28);
            GStyle.getAxisAttributesY().setTitleFontSize(28);
            GStyle.getAxisAttributesX().setLabelFontSize(24);
            GStyle.getAxisAttributesY().setLabelFontSize(24);
            for (int ir = 0; ir < NREGIONS; ir++) {
                EmbeddedPad pad = this.getCanvas(name).getPad(ir);
                // pad.getAxisY().setRange(0, Region_OCC_Max);
                pad.getAxisX().setRange(0.5, 6.5);
                pad.setTitle("Region " + (ir + 1));
                // GraphErrors datasets, so setLegend() helper (which only recognizes H1F)
                // won't enable this - set it directly instead
                pad.setLegend(true);
                pad.setLegendPosition(250, 70);
            }
        }

        else if (name.equals("Region Occupancy - Species (Origin)")) {
            GStyle.getAxisAttributesX().setTitleFontSize(28);
            GStyle.getAxisAttributesY().setTitleFontSize(28);
            GStyle.getAxisAttributesX().setLabelFontSize(24);
            GStyle.getAxisAttributesY().setLabelFontSize(24);
            for (int ir = 0; ir < NREGIONS; ir++) {
                EmbeddedPad pad = this.getCanvas(name).getPad(ir);
                // pad.getAxisY().setRange(0, Region_OCC_Max);
                pad.getAxisX().setRange(0.5, 6.5);
                pad.setTitle("Region " + (ir + 1));
                pad.setLegend(true);
                pad.setLegendPosition(250, 70);
            }
        }

        // ← Sector check BEFORE the generic "Origin of Bg" check
        else if (name.contains("Origin of Bg - Sector")) {
            GStyle.getAxisAttributesX().setTitleFontSize(20);
            GStyle.getAxisAttributesY().setTitleFontSize(20);
            GStyle.getAxisAttributesX().setLabelFontSize(16);
            GStyle.getAxisAttributesY().setLabelFontSize(16);
            // for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
            // pad.getAxisY().setRange(0, 500);
            // }
            this.setLegend(name, 250, 70);
        }

        else if (name.equals("Origin of Bg") || name.equals("Origin of Bg - Energy")) {
            GStyle.getAxisAttributesX().setTitleFontSize(20);
            GStyle.getAxisAttributesY().setTitleFontSize(20);
            GStyle.getAxisAttributesX().setLabelFontSize(16);
            GStyle.getAxisAttributesY().setLabelFontSize(16);
            this.setupColorAxis(name, "Rate [kHz]");
            for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
                // pad.getAxisZ().setRange(1, BG_Y_MAX);
                if (pad.getDatasetPlotters().size() == 0)
                    continue;
                IDataSet ds = pad.getDatasetPlotters().get(0).getDataSet();
                if (ds instanceof H1F) {
                    if (name.contains("Energy")) {
                        pad.getAxisY().setLog(true);
                        // pad.getAxisY().setRange(1, 300);
                        // } else {
                        // pad.getAxisY().setRange(0, 300);
                    }

                    this.setLegend(name, 250, 70);

                }
            }
        }

        else if (name.equals("Position of Bg")) {
            GStyle.getAxisAttributesX().setTitleFontSize(22);
            GStyle.getAxisAttributesY().setTitleFontSize(22);
            GStyle.getAxisAttributesX().setLabelFontSize(22);
            GStyle.getAxisAttributesY().setLabelFontSize(22);
            GStyle.getAxisAttributesZ().setTitleFontSize(22);
            GStyle.getAxisAttributesZ().setLabelFontSize(22);
            // for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
            // pad.getAxisZ().setRange(1, POS_Z_MAX);
            // }
            this.setupColorAxis(name, "Rate [kHz]");
            for (int ir = 0; ir < NREGIONS; ir++) {
                EmbeddedPad padBot = this.getCanvas(name).getPad(ir + 3);
                padBot.getAxisX().setRange(-DR[ir] / 1000, DR[ir] / 1000);
                padBot.getAxisY().setRange(-DR[ir] / 1000, DR[ir] / 1000);
                padBot.getAxisX().setAxisDivisions(5);
                padBot.getAxisY().setAxisDivisions(5);
            }
        }

        // if (!name.contains("Occupancy"))
        // this.setLogZ(name);
    }

    @Override
    public void normalizeToTime(DataGroup dg) {
        double factor = this.getNevents() * (Constants.getTimeWindow() * 1E-9) / 1E-6; // Mz

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
