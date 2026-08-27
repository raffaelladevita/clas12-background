package modules;

import analysis.Constants;
import analysis.Module;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import objects.Event;
import objects.Hit;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.LatexText;

/**
 * Background module for the CLAS12 Forward Time-of-Flight (FTOF) detector.
 *
 * FTOF geometry (MC::True detector ID = 12):
 * Layer 1 = Panel 1a (23 paddles per sector)
 * Layer 2 = Panel 1b (62 paddles per sector)
 * Layer 3 = Panel 2 ( 5 paddles per sector)
 * 6 sectors total
 *
 * Plots produced:
 * Sector-Layer Map — 2D occupancy [%]: sector vs layer (from FTOF::adc/tdc)
 * Rate by Sector — rate [kHz] vs sector, one histogram per panel (from
 * FTOF::adc/tdc)
 * Paddle Map — 2D occupancy [%]: paddle vs layer, one pad per sector (from
 * FTOF::adc/tdc)
 * Position Map — Y vs X and R vs Z from MC::True avgX/Y/Z (always filled when
 * MC::True present)
 *
 * The first three groups require a digitized FTOF::adc or FTOF::tdc bank.
 * The Position Map is filled directly from MC::True and works for pure-MC
 * studies.
 * All position coordinates are in mm.
 */
public class FTOFmodule extends Module {

    private static final int NSECTORS = 6;
    private static final int NLAYERS = 3; // panels: 1=1a, 2=1b, 3=P2
    private static final int NPADDLES = 62; // max paddles (panel 1b)
    private static final String[] PANELS = { "Panel1a", "Panel1b", "Panel2" };

    // Concatenated x-axis for the Rate Map: Panel-1B | Panel-1A | Panel-2
    // Layer 2 (1b, 62 paddles): bins 1-62 offset=0
    // Layer 1 (1a, 23 paddles): bins 71-93 offset=70
    // Layer 3 (P2, 5 paddles): bins 101-106 offset=100
    // PADDLE_XOFFSET index = layer-1
    private static final int[] PADDLE_XOFFSET = { 70, 0, 100 };
    private static final int RATEMAP_NBINS = 120;

    // Diagnostic counters — printed once in analyzeHistos()
    private long diagTotal = 0, diagOrder1 = 0, diagLowEdep = 0, diagFilled = 0, diagVxVyVeto = 0,
            diagVxVyVetoWouldBeHit = 0;
    // Per-panel breakdown of the two counters above: index 0=P1a, 1=P1b, 2=P2.
    private final long[] diagFilledByPanel = new long[3];
    private final long[] diagVxVyVetoWouldBeHitByPanel = new long[3];
    // Per-panel-per-sector breakdown: [panel-1][sector-1].
    private final long[][] diagFilledByPanelSector = new long[3][NSECTORS];
    private final long[][] diagVxVyVetoWouldBeHitByPanelSector = new long[3][NSECTORS];

    // Carman formula: current [µA] = Edep[MeV] × yield[ph/MeV] × 1.6e-7[µC/ph] /
    // (2h[cm])
    // MC::True::totEdep is in MeV. Yields and half-thicknesses from Carman (CLAS12
    // FTOF note 2014):
    // P1a: BC-408 + EMI/Philips PMTs, yield=373 ph/MeV, 2h=10 cm
    // P1b: BC-404/BC-408 + Hamamatsu R9779, yield=1158 ph/MeV, 2h=12 cm
    // P2: BC-408 + Photonis/EMI PMTs, yield=373 ph/MeV, 2h=10 cm ← same
    // scintillator/thickness
    // as P1a; update if Carman note gives a separate P2 yield.
    private static final double YIELD_P1A = 373.0;
    private static final double YIELD_P1B = 1158.0;
    private static final double YIELD_P2 = 373.0; // BC-408, Photonis PMTs — same as P1a pending measurement
    private static final double TWO_H_P1A = 10.0; // 2 × 5 cm
    private static final double TWO_H_P1B = 12.0; // 2 × 6 cm
    private static final double TWO_H_P2 = 10.0; // 2 × 5 cm

    // Physical extent of FTOF in metres (positions from MC::True are in mm → divide
    // by 1000 when filling)
    private static final double FTOF_Z_MIN = 3.5;
    private static final double FTOF_Z_MAX = 9.0;
    private static final double FTOF_R_MAX = 6.0;
    private static final double FTOF_XY_MAX = 5.5;

    // Vertex-origin ranges for origin_bg (m) — vertices from MC::True are in mm,
    // divided by 1000 before filling so all origin plots use metres.
    private static final double ORIG_VZ_MIN = -1.0; // m
    private static final double ORIG_VZ_MAX = 9.0; // m
    private static final double ORIG_DR = 6.0; // m
    private static final double ORIG_R_MIN = -1.0; // m (negative space for detector labels)
    // Fixed bin width keeps Rate [kHz] per bin independent of zoom level.
    private static final double ORIG_VZ_BIN_WIDTH = 0.05; // m
    private static final int ORIG_VZ_NBINS = (int) Math.round((ORIG_VZ_MAX - ORIG_VZ_MIN) / ORIG_VZ_BIN_WIDTH);

    public FTOFmodule() {
        super(DetectorType.FTOF);
    }

    // -------------------------------------------------------------------------
    // Histogram group definitions
    // -------------------------------------------------------------------------

    /** 2D occupancy [%] map: sector (x) vs layer/panel (y). */
    public DataGroup sectorLayerMap() {
        DataGroup dg = new DataGroup(1, 1);
        H2F hi = histo2D("hi_sec_lay", "Occupancy (%)",
                "Sector", "Layer  (1=1a  2=1b  3=P2)",
                NSECTORS, 0.5, NSECTORS + 0.5,
                NLAYERS, 0.5, NLAYERS + 0.5);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /** Hit rate [kHz] vs sector — one histogram per panel. */
    public DataGroup rateBySector() {
        DataGroup dg = new DataGroup(3, 1);
        int[] colors = { 4, 2, 3 };
        for (int i = 0; i < NLAYERS; i++) {
            H1F hi = histo1D("hi_rate_sec_" + PANELS[i], PANELS[i],
                    "Sector", "Rate [kHz]",
                    NSECTORS, 0.5, NSECTORS + 0.5, 0);
            hi.setLineColor(colors[i]);
            hi.setLineWidth(2);
            dg.addDataSet(hi, i);
        }
        return dg;
    }

    /**
     * Per-sector 2D occupancy map: paddle/component (x) vs layer/panel (y).
     * Six pads, one per sector.
     */
    public DataGroup paddleMap() {
        DataGroup dg = new DataGroup(3, 2);
        for (int s = 1; s <= NSECTORS; s++) {
            H2F hi = histo2D("hi_paddle_s" + s, "Occupancy (%)",
                    "Paddle", "Layer  (1=1a  2=1b  3=P2)",
                    NPADDLES, 0.5, NPADDLES + 0.5,
                    NLAYERS, 0.5, NLAYERS + 0.5);
            hi.setTitle("Sector " + s);
            dg.addDataSet(hi, s - 1);
        }
        return dg;
    }

    /**
     * 2D rate map [kHz]: all panels concatenated on x (1B | 1A | P2) vs sector on
     * y.
     * Mirrors Figure 7 of CLAS12 Note 2017-016.
     * X layout: Panel-1B paddles 1-62 (bins 1-62), gap, Panel-1A paddles 1-23
     * (bins 73-95), gap, Panel-2 paddles 1-5 (bins 107-111).
     */
    public DataGroup rateMap() {
        DataGroup dg = new DataGroup(1, 1);
        H2F hi = histo2D("hi_rate_map", "Rate [kHz]",
                "Paddle  (1B  |  1A  |  P2)", "Sector",
                RATEMAP_NBINS, 0.5, RATEMAP_NBINS + 0.5,
                NSECTORS, 0.5, NSECTORS + 0.5);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /**
     * 2D PMT current map [µA]: same paddle layout as Rate Map (x) vs sector (y).
     * Weight per hit: Edep_MeV × photon_yield × 1.6e-7 [nC/ph] / (2×h_cm).
     * Photon yields: P1a/P2 = 373 ph/MeV, P1b = 1158 ph/MeV.
     * 2×h: P1a=10 (2×5cm thickness), P1b=12 (2×6cm thickness), P2=10 (2×5cm
     * thickness).
     */
    public DataGroup pmtCurrentMap() {
        DataGroup dg = new DataGroup(1, 1);
        H2F hi = histo2D("hi_pmt_curr", "PMT Current [#muA]",
                "Paddle  (1B  |  1A  |  P2)", "Sector",
                RATEMAP_NBINS, 0.5, RATEMAP_NBINS + 0.5,
                NSECTORS, 0.5, NSECTORS + 0.5);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /**
     * 2D position maps from MC::True avgX/Y/Z.
     * Pad 0: Y vs X (transverse plane).
     * Pad 1: R vs Z (longitudinal view), R = sqrt(X²+Y²).
     */
    public DataGroup positionMap() {
        DataGroup dg = new DataGroup(2, 1);
        H2F hiXY = histo2D("hi_posXY", "Y vs X",
                "X [m]", "Y [m]",
                200, -FTOF_XY_MAX, FTOF_XY_MAX,
                200, -FTOF_XY_MAX, FTOF_XY_MAX);
        H2F hiRZ = histo2D("hi_posRZ", "R vs Z",
                "Z [m]", "R [m]",
                200, FTOF_Z_MIN, FTOF_Z_MAX,
                200, 0, FTOF_R_MAX);
        dg.addDataSet(hiXY, 0);
        dg.addDataSet(hiRZ, 1);
        return dg;
    }

    /**
     * Origin of background: 2D vertex maps (Vx/Vy and Vz/r) plus 1D
     * distributions per particle species, one set of pads per FTOF panel.
     * Layout (3 columns × 3 rows):
     * Row 0 (pads 0-2): Vx vs Vy per panel
     * Row 1 (pads 3-5): Vz vs r per panel
     * Row 2 (pads 6-8): 1D by species per panel — Vz for index 0, kinetic
     * energy for index 1.
     * Index 0 is filled unweighted (counts); index 1 is filled weighted by
     * kinetic energy so that, after dividing 1 by 0 in analyzeHistos(), the 2D
     * maps become mean-energy-per-bin (mirrors DCmodule's "Origin of Bg -
     * Energy").
     */
    public DataGroup[] origin_bg() {
        DataGroup[] dg = new DataGroup[2];
        for (int i = 0; i < dg.length; i++) {
            dg[i] = new DataGroup(3, 3);
            for (int il = 0; il < NLAYERS; il++) {
                int panel = il + 1;
                H2F hi_xy = histo2D("hi_bg_origin_xy_panel" + panel,
                        "Vx(m)", "Vy(m)",
                        200, -ORIG_DR, ORIG_DR,
                        200, -ORIG_DR, ORIG_DR);
                H2F hi_rz = histo2D("hi_bg_origin_rz_panel" + panel,
                        "Vz(m)", "r(m)",
                        ORIG_VZ_NBINS, ORIG_VZ_MIN, ORIG_VZ_MAX,
                        200, ORIG_R_MIN, ORIG_DR);
                dg[i].addDataSet(hi_xy, 0 + il);
                dg[i].addDataSet(hi_rz, 3 + il);

                double min = ORIG_VZ_MIN;
                double max = ORIG_VZ_MAX;
                String xname = "Vz(m)";
                int nbins = ORIG_VZ_NBINS;
                if (i == 1) {
                    min = -1.0;
                    max = 200;
                    xname = "E(MeV)";
                    nbins = 200;
                }
                for (int ip = 0; ip < PNAMES.length; ip++) {
                    H1F hi_bg = histo1D("hi_bg_panel" + panel + "_" + PNAMES[ip],
                            PNAMES[ip], xname, "Rate [kHz]",
                            nbins, min, max, 0);
                    this.setHistoAttr(hi_bg, ip < 5 ? ip + 1 : ip + 3);
                    dg[i].addDataSet(hi_bg, 6 + il);
                }
            }
        }
        return dg;
    }

    /**
     * Species-isolated origin maps (Vx/Vy and Vz/r), one set of pads per FTOF
     * panel — lets a single species (e.g. "electron") be inspected in
     * isolation instead of only as one line among several in "Origin of Bg".
     * Layout (3 columns × 2 rows): row 0 = Vx vs Vy, row 1 = Vz vs r.
     */
    public DataGroup speciesOriginMap(String species) {
        DataGroup dg = new DataGroup(3, 2);
        for (int il = 0; il < NLAYERS; il++) {
            int panel = il + 1;
            H2F hi_xy = histo2D("hi_bg_origin_xy_panel" + panel + "_" + species,
                    "Vx(m)", "Vy(m)",
                    200, -ORIG_DR, ORIG_DR,
                    200, -ORIG_DR, ORIG_DR);
            H2F hi_rz = histo2D("hi_bg_origin_rz_panel" + panel + "_" + species,
                    "Vz(m)", "r(m)",
                    ORIG_VZ_NBINS, ORIG_VZ_MIN, ORIG_VZ_MAX,
                    200, ORIG_R_MIN, ORIG_DR);
            dg.addDataSet(hi_xy, 0 + il);
            dg.addDataSet(hi_rz, 3 + il);
        }
        return dg;
    }

    /**
     * Per-sector 1D Vz background origin for one FTOF panel.
     * Layout: DataGroup(2,3) — 6 pads, one per sector.
     */
    public DataGroup sectorBG(int panel) {
        DataGroup dg = new DataGroup(2, 3);
        for (int is = 0; is < NSECTORS; is++) {
            int sector = is + 1;
            for (int ip = 0; ip < PNAMES.length; ip++) {
                H1F hi_bg = histo1D(
                        "hi_bg_p" + panel + "_s" + sector + "_" + PNAMES[ip],
                        PANELS[panel - 1] + " S" + sector + "-" + PNAMES[ip],
                        "Vz(m)", "Rate [kHz]",
                        ORIG_VZ_NBINS, ORIG_VZ_MIN, ORIG_VZ_MAX, 0);
                this.setHistoAttr(hi_bg, ip < 5 ? ip + 1 : ip + 3);
                dg.addDataSet(hi_bg, is);
            }
        }
        return dg;
    }

    public void fillSectorBG(DataGroup group, List<Hit> hits, int targetPanel) {
        for (Hit hit : hits) {
            True t = hit.getTrue();
            if (t == null)
                continue;
            int panel = hit.getLayer();
            int sector = hit.getSector();
            if (panel != targetPanel)
                continue;

            double vz = t.getVertex().z() / 1000.0; // mm → m
            group.getH1F("hi_bg_p" + panel + "_s" + sector + "_all").fill(vz);
            if (this.pidToName(Math.abs(t.getPid())) != null)
                group.getH1F("hi_bg_p" + panel + "_s" + sector + "_"
                        + this.pidToName(Math.abs(t.getPid()))).fill(vz);
            else
                group.getH1F("hi_bg_p" + panel + "_s" + sector + "_other").fill(vz);
        }
    }

    public void fillOrigin(DataGroup group, List<Hit> hits, boolean energyWeight) {
        for (Hit hit : hits) {
            True t = hit.getTrue();
            if (t == null)
                continue;
            int panel = hit.getLayer();

            double vx = t.getVertex().x() / 1000.0; // mm → m
            double vy = t.getVertex().y() / 1000.0;
            double vz = t.getVertex().z() / 1000.0;
            double r = Math.sqrt(vx * vx + vy * vy);
            double weight = energyWeight ? t.getKinEnergy() : 1;

            group.getH2F("hi_bg_origin_rz_panel" + panel).fill(vz, r, weight);
            group.getH2F("hi_bg_origin_xy_panel" + panel).fill(vx, vy, weight);

            double value = energyWeight ? t.getKinEnergy() : vz;
            group.getH1F("hi_bg_panel" + panel + "_all").fill(value);
            if (this.pidToName(Math.abs(t.getPid())) != null)
                group.getH1F("hi_bg_panel" + panel + "_" + this.pidToName(Math.abs(t.getPid()))).fill(value);
            else
                group.getH1F("hi_bg_panel" + panel + "_other").fill(value);
        }
    }

    public void fillSpeciesOrigin(DataGroup group, List<Hit> hits, String species) {
        for (Hit hit : hits) {
            True t = hit.getTrue();
            if (t == null)
                continue;
            if (!species.equals(this.pidToName(Math.abs(t.getPid()))))
                continue;
            int panel = hit.getLayer();

            double vx = t.getVertex().x() / 1000.0;
            double vy = t.getVertex().y() / 1000.0;
            double vz = t.getVertex().z() / 1000.0;
            double r = Math.sqrt(vx * vx + vy * vy);

            group.getH2F("hi_bg_origin_rz_panel" + panel + "_" + species).fill(vz, r);
            group.getH2F("hi_bg_origin_xy_panel" + panel + "_" + species).fill(vx, vy);
        }
    }

    // -------------------------------------------------------------------------
    // Module interface
    // -------------------------------------------------------------------------

    @Override
    public void createHistos() {
        this.getHistos().put("Sector-Layer Map", this.sectorLayerMap());
        this.getHistos().put("Rate by Sector", this.rateBySector());
        this.getHistos().put("Paddle Map", this.paddleMap());
        this.getHistos().put("Rate Map", this.rateMap());
        this.getHistos().put("PMT Current Map", this.pmtCurrentMap());
        this.getHistos().put("Position Map", this.positionMap());
        DataGroup[] originGroups = this.origin_bg();
        this.getHistos().put("Origin of Bg", originGroups[0]);
        this.getHistos().put("Origin of Bg - Energy", originGroups[1]);
        this.getHistos().put("Origin of Bg - electron", this.speciesOriginMap("electron"));
        this.getHistos().put("Origin of Bg - gamma", this.speciesOriginMap("gamma"));
        this.getHistos().put("Origin of Bg - Sector P1a", this.sectorBG(1));
        this.getHistos().put("Origin of Bg - Sector P1b", this.sectorBG(2));
        this.getHistos().put("Origin of Bg - Sector P2", this.sectorBG(3));
    }

    @Override
    public void fillHistos(Event event) {
        List<Hit> allHits = event.getHits(DetectorType.FTOF);
        if (allHits != null) {
            // hits: order=0, edep>1 MeV (hardware discriminator), used for
            // rate/occupancy/origin.
            // hitsForCurrent: order=0, no edep cut — all energy deposits contribute to PMT
            // charge.
            List<Hit> hits = new ArrayList<>();
            List<Hit> hitsForCurrent = new ArrayList<>();
            for (Hit h : allHits) {
                diagTotal++;
                if (h.getOrder() != 0) {
                    diagOrder1++;
                    continue;
                }
                hitsForCurrent.add(h);
                if (h.getTrue() == null || h.getTrue().getEdep() <= 1.0) {
                    diagLowEdep++;
                    continue;
                }
                hits.add(h);
                if (h.getLayer() >= 1 && h.getLayer() <= 3) {
                    diagFilledByPanel[h.getLayer() - 1]++;
                    if (h.getSector() >= 1 && h.getSector() <= NSECTORS)
                        diagFilledByPanelSector[h.getLayer() - 1][h.getSector() - 1]++;
                }
            }
            diagFilled += hits.size();

            DataGroup dgMap = this.getHistos().get("Sector-Layer Map");
            DataGroup dgRate = this.getHistos().get("Rate by Sector");
            DataGroup dgPaddle = this.getHistos().get("Paddle Map");
            DataGroup dgRateMap = this.getHistos().get("Rate Map");
            DataGroup dgCurr = this.getHistos().get("PMT Current Map");

            // Max paddles per panel: P1a=23, P1b=62, P2=5
            final int[] MAX_PADDLES = { 23, 62, 5 };

            for (Hit h : hits) {
                int sector = h.getSector();
                int layer = h.getLayer();
                int paddle = h.getComponent();
                if (sector < 1 || sector > NSECTORS)
                    continue;
                if (layer < 1 || layer > NLAYERS)
                    continue;
                if (paddle < 1 || paddle > MAX_PADDLES[layer - 1])
                    continue;

                dgMap.getH2F("hi_sec_lay").fill(sector, layer);
                dgRate.getH1F("hi_rate_sec_" + PANELS[layer - 1]).fill(sector);
                dgPaddle.getH2F("hi_paddle_s" + sector).fill(paddle, layer);
                dgRateMap.getH2F("hi_rate_map").fill(paddle + PADDLE_XOFFSET[layer - 1], sector);
            }

            // PMT current: no edep threshold — every deposited photon contributes to
            // charge.
            for (Hit h : hitsForCurrent) {
                if (h.getTrue() == null)
                    continue;
                int sector = h.getSector();
                int layer = h.getLayer();
                int paddle = h.getComponent();
                if (sector < 1 || sector > NSECTORS)
                    continue;
                if (layer < 1 || layer > NLAYERS)
                    continue;
                if (paddle < 1 || paddle > MAX_PADDLES[layer - 1])
                    continue;

                double edep_MeV = h.getTrue().getEdep();
                double weightCurrent;
                if (layer == 1) {
                    weightCurrent = edep_MeV * YIELD_P1A * 1.6e-7 / TWO_H_P1A;
                } else if (layer == 2) {
                    weightCurrent = edep_MeV * YIELD_P1B * 1.6e-7 / TWO_H_P1B;
                } else {
                    weightCurrent = edep_MeV * YIELD_P2 * 1.6e-7 / TWO_H_P2;
                }
                dgCurr.getH2F("hi_pmt_curr").fill(paddle + PADDLE_XOFFSET[layer - 1], sector, weightCurrent);
            }

            this.fillOrigin(this.getHistos().get("Origin of Bg"), hits, false);
            this.fillOrigin(this.getHistos().get("Origin of Bg - Energy"), hits, true);
            this.fillSpeciesOrigin(this.getHistos().get("Origin of Bg - electron"), hits, "electron");
            this.fillSpeciesOrigin(this.getHistos().get("Origin of Bg - gamma"), hits, "gamma");
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector P1a"), hits, 1);
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector P1b"), hits, 2);
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector P2"), hits, 3);

            // Position maps — edep > 1 MeV threshold applied via hits list
            DataGroup dgPos = this.getHistos().get("Position Map");
            for (Hit h : hits) {
                if (h.getTrue() == null)
                    continue;
                double x = h.getTrue().getPosition().x() / 1000.0;
                double y = h.getTrue().getPosition().y() / 1000.0;
                double z = h.getTrue().getPosition().z() / 1000.0;
                double r = Math.sqrt(x * x + y * y);
                dgPos.getH2F("hi_posXY").fill(x, y);
                dgPos.getH2F("hi_posRZ").fill(z, r);
            }
        }
    }

    @Override
    public void analyzeHistos() {
        double x = Double.parseDouble(System.getProperty("lumi", "117000.0"));
        double lumiScale = x / 117000.0;

        // Occupancy [%] = hits / (events * 0.01), then lumi-scaled
        this.normalizeToEventsX100(
                this.getHistos().get("Sector-Layer Map").getH2F("hi_sec_lay"));
        this.normalize(this.getHistos().get("Sector-Layer Map"), lumiScale);

        for (int s = 1; s <= NSECTORS; s++) {
            this.normalizeToEventsX100(
                    this.getHistos().get("Paddle Map").getH2F("hi_paddle_s" + s));
        }
        this.normalize(this.getHistos().get("Paddle Map"), lumiScale);

        // Rates [kHz], lumi-scaled
        this.normalizeToTime(this.getHistos().get("Rate by Sector"));
        this.normalize(this.getHistos().get("Rate by Sector"), lumiScale);

        // 2D rate map [kHz], lumi-scaled
        this.normalizeToTime(this.getHistos().get("Rate Map"));
        this.normalize(this.getHistos().get("Rate Map"), lumiScale);

        // PMT current map [µA]: weight is in µC (Carman formula), divide by time_s →
        // µC/s = µA
        this.normalizeToTime(this.getHistos().get("PMT Current Map").getH2F("hi_pmt_curr"), 1);
        this.normalize(this.getHistos().get("PMT Current Map"), lumiScale);

        // Position map [kHz/bin], lumi-scaled
        this.normalizeToTime(this.getHistos().get("Position Map"));
        this.normalize(this.getHistos().get("Position Map"), lumiScale);

        // Origin of Bg - Energy: convert the energy-weighted sum per bin into mean
        // energy per bin by dividing by the (not-yet-normalized) counts version —
        // normalizing either side first would cancel out in the ratio anyway.
        this.divide(this.getHistos().get("Origin of Bg - Energy"), this.getHistos().get("Origin of Bg"));

        // Origin of Bg [kHz/bin], lumi-scaled
        this.normalizeToTime(this.getHistos().get("Origin of Bg"));
        this.normalize(this.getHistos().get("Origin of Bg"), lumiScale);

        this.normalizeToTime(this.getHistos().get("Origin of Bg - Energy"));

        for (String species : new String[] { "electron", "gamma" }) {
            this.normalizeToTime(this.getHistos().get("Origin of Bg - " + species));
            this.normalize(this.getHistos().get("Origin of Bg - " + species), lumiScale);
        }

        for (String panel : new String[] { "P1a", "P1b", "P2" }) {
            this.normalizeToTime(this.getHistos().get("Origin of Bg - Sector " + panel));
            this.normalize(this.getHistos().get("Origin of Bg - Sector " + panel), lumiScale);
        }

    }

    /**
     * Draws colored vertical tick lines at each major CLAS12 detector z-position
     * and places the detector label in the negative-r band (r = 0 to -1),
     * staggered across three depths to prevent overlap between close detectors.
     */
    private static final int DET_COLOR_BASE = 100;
    private static final java.awt.Color[] DET_COLORS = {
            new java.awt.Color(0, 0, 0), // 100 Tgt — black
            new java.awt.Color(160, 0, 0), // 101 HTCC M — dark red
            new java.awt.Color(0, 0, 160), // 102 DC-R1 — dark blue
            new java.awt.Color(0, 110, 0), // 103 DC-R2 — dark green
            new java.awt.Color(90, 0, 140), // 104 DC-R3 — dark violet
            new java.awt.Color(140, 55, 0), // 105 FTOF 1b — dark brown
            new java.awt.Color(0, 95, 95), // 106 FTOF 1a — dark teal
            new java.awt.Color(120, 0, 80), // 107 PCAL — dark maroon
            new java.awt.Color(70, 90, 0), // 108 EC — dark olive
    };

    static {
        // Inject custom colors into ColorPalette's static TreeMap once per JVM.
        try {
            java.lang.reflect.Field f = org.jlab.groot.base.ColorPalette.class.getDeclaredField("colorPalette");
            f.setAccessible(true);
            @SuppressWarnings("unchecked")
            java.util.TreeMap<Integer, java.awt.Color> map = (java.util.TreeMap<Integer, java.awt.Color>) f.get(null);
            for (int i = 0; i < DET_COLORS.length; i++)
                map.put(DET_COLOR_BASE + i, DET_COLORS[i]);
        } catch (Exception e) {
            System.err.println("FTOFmodule: could not register custom palette colors: " + e.getMessage());
        }
    }

    private void addDetectorLabels(String canvasName) {
        // DSShieldFrontLead lab z from STL: front face 4960 mm, back face 5322 mm (fc
        // frame + 180-deg Y rotation)
        final double[] detZ = { 0.0, 1.6, 2.5, 3.9, 5.4, 7.2, 7.3, 7.7, 8.0, 4.960, 5.322 };
        final String[] detLabel = { "Tgt", "HTCC M", "DC-R1", "DC-R2", "DC-R3",
                "FTOF 1b", "FTOF 1a", "PCAL", "EC" };

        EmbeddedCanvas canvas = this.getCanvas(canvasName);

        // X_LEGEND / Y_START are in pad-pixel coords (include margins).
        // Frame (gray area) starts at roughly x≈70 px, y≈35 px.
        // Increase X_LEGEND to move labels right; increase Y_START to move down.
        final int FONT = 11;
        final int LINE_H = 13;
        final int X_LEGEND = 110; // px from pad left — puts labels just inside frame
        final int Y_START = 60; // px from pad top — puts labels just inside frame

        for (int padIdx = 3; padIdx <= 5; padIdx++) {
            canvas.cd(padIdx);

            // All ticks same full depth — no staggering
            for (int i = 0; i < detZ.length; i++) {
                double z = detZ[i];
                if (z < ORIG_VZ_MIN || z > ORIG_VZ_MAX)
                    continue;
                DataLine tick = new DataLine(z, 0.0, z, ORIG_R_MIN);
                tick.setLineColor(DET_COLOR_BASE + i);
                tick.setLineWidth(3);
                canvas.draw(tick);
            }

            for (int i = 0; i < detLabel.length; i++) {
                LatexText lbl = new LatexText(detLabel[i], X_LEGEND, Y_START + i * LINE_H);
                lbl.setFontSize(FONT);
                lbl.setColor(DET_COLOR_BASE + i);
                canvas.draw(lbl);
            }
        }
    }

    @Override
    public void setPlottingOptions(String name) {
        GStyle.getAxisAttributesX().setTitleFontSize(20);
        GStyle.getAxisAttributesY().setTitleFontSize(20);
        GStyle.getAxisAttributesX().setLabelFontSize(20);
        GStyle.getAxisAttributesY().setLabelFontSize(20);
        GStyle.getAxisAttributesZ().setTitleFontSize(20);
        GStyle.getAxisAttributesZ().setLabelFontSize(20);

        for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
            pad.setTitle("");

        if (name.equals("Sector-Layer Map") || name.equals("Paddle Map")) {
            this.setupColorAxis(name, "Occupancy [%]");
            this.setWideColorBar(name);
        }

        if (name.equals("Position Map")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
            // this.setLogZ(name);
        }

        if (name.equals("Rate Map")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
        }

        if (name.equals("PMT Current Map")) {
            this.setupColorAxis(name, "Current [#muA]");
            this.setWideColorBar(name);
        }

        if (name.equals("Rate by Sector"))
            this.setLegend(name, 250, 50);

        if (name.contains("Origin of Bg - Sector")) {
            GStyle.getAxisAttributesX().setTitleFontSize(20);
            GStyle.getAxisAttributesY().setTitleFontSize(20);
            GStyle.getAxisAttributesX().setLabelFontSize(16);
            GStyle.getAxisAttributesY().setLabelFontSize(16);
            this.setLegend(name, 350, 70);
        }

        if (name.equals("Origin of Bg")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
            this.setLogZ(name);
            for (EmbeddedPad p : this.getCanvas(name).getCanvasPads())
                p.getAxisZ().setRange(0.1, 500);
            // // Fix y-axis of 1D Vz rate plots (row 2: pads 6=P1a, 7=P1b, 8=P2)
            // List<EmbeddedPad> pads = this.getCanvas(name).getCanvasPads();
            // pads.get(6).getAxisY().setRange(0, 7000);
            // pads.get(7).getAxisY().setRange(0, 7000);
            // pads.get(8).getAxisY().setRange(0, 5000);
            this.setLegend(name, 200, 70);
            this.addDetectorLabels(name);
        }

        if (name.equals("Origin of Bg - Energy")) {
            this.setupColorAxis(name, "Mean E [MeV]");
            this.setWideColorBar(name);
            this.setLogZ(name);
            for (EmbeddedPad p : this.getCanvas(name).getCanvasPads())
                p.getAxisZ().setRange(0.1, 200);
            // Row 2 (pads 6=P1a, 7=P1b, 8=P2) are energy spectra — log-y, same range
            List<EmbeddedPad> pads = this.getCanvas(name).getCanvasPads();
            for (int idx = 6; idx <= 8; idx++) {
                pads.get(idx).getAxisY().setLog(true);
                pads.get(idx).getAxisY().setRange(1, 500);
            }
            this.setLegend(name, 200, 70);
            this.addDetectorLabels(name);
        }

        if (name.equals("Origin of Bg - electron") || name.equals("Origin of Bg - gamma")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
            this.setLogZ(name);
            for (EmbeddedPad p : this.getCanvas(name).getCanvasPads())
                p.getAxisZ().setRange(0.1, 500);
            this.addDetectorLabels(name);
        }

        this.saveCanvas(name);
    }
}
