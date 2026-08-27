package modules;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import analysis.Constants;
import analysis.Module;
import objects.Event;
import objects.Hit;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.LatexText;

/**
 * Background module for the CLAS12 electromagnetic calorimeter system.
 *
 * All three sub-detectors share the ECAL::adc bank, distinguished by layer:
 * layers 1-3 → PCAL (preshower, U/V/W views)
 * layers 4-6 → ECIN (inner EC, U/V/W views)
 * layers 7-9 → ECOUT (outer EC, U/V/W views)
 *
 * Strip counts per view:
 * PCAL U=68, V=62, W=62
 * ECIN U=36, V=36, W=36
 * ECOUT U=36, V=36, W=36
 *
 * Plots produced:
 * Sector-Layer Map — 2D occupancy [%]: sector (x) vs layer 1-9 (y)
 * Rate by Sector — 3 pads: kHz vs sector for PCAL / ECIN / ECOUT
 * PCAL Strip Map — 6 pads (one per sector): strip (x) vs view 1-3 (y)
 * ECIN Strip Map — same for ECIN
 * ECOUT Strip Map — same for ECOUT
 * Position Map — 6 pads: R vs Z (top) and Y vs X (bottom) per sub-detector
 * from MC::True avgX/avgY/avgZ; empty if MC::True not in file
 * Rate Map — 3 pads: strip (x) vs sector (y), rate [kHz], one per sub-detector
 * PMT Current Map — 3 pads: same layout, weighted PMT current [µA]
 * Origin of Bg — 3×3: Vx/Vy (row 0), Vz/r (row 1), 1D Vz species (row 2)
 * per sub-detector column (PCAL | ECIN | ECOUT)
 */
public class ECmodule extends Module {

    private static final int NSECTORS = 6;
    private static final int NLAYERS = 9; // 3 sub-dets × 3 views
    private static final int NSTRIPS = 70; // generous upper bound (PCAL U = 68)
    private static final String[] SUBDETS = { "PCAL", "ECIN", "ECOUT" };

    // Physical extent of ECAL in metres (MC::True positions are mm → divide by 1000
    // when filling)
    private static final double EC_Z_MIN = 5.0;
    private static final double EC_Z_MAX = 9.5;
    private static final double EC_R_MAX = 5.50;
    private static final double EC_XY_MAX = 5.50;

    // Vertex-origin ranges for origin_bg (m)
    private static final double ORIG_VZ_MIN = -1.0;
    private static final double ORIG_VZ_MAX = 10.0;
    private static final double ORIG_DR = 7.0;
    private static final double ORIG_R_MIN = -1.0; // negative band holds detector labels

    // PMT current: Carman-style formula — current[µA] = Edep[MeV] × yield[ph/MeV] ×
    // 1.6e-7 / (2h[cm])
    // Yields and half-thicknesses are placeholder estimates; update from EC
    // technical notes.
    // PCAL: BC-408-like scintillator strips + WLS fibers + SiPMs
    // ECIN/ECOUT: BC-412 scintillator + Photonis XP2262 PMTs
    private static final double[] YIELD = { 500.0, 800.0, 800.0 }; // ph/MeV — PCAL, ECIN, ECOUT
    private static final double[] TWO_H = { 2.0, 4.4, 4.4 }; // 2 × half-thickness [cm]

    // Minimum edep [MeV] — MC::True::totEdep is in MeV (same as FTOFmodule).
    // Set to 0.0 to skip only zero/negative entries; raise to 1.0 to match the
    // 1 MeV threshold used in RGA background notes.
    private static final double EDEP_THRESHOLD = 1.0; // MeV

    // Custom dark colors registered into GROOT's ColorPalette — same indices as
    // FTOFmodule
    // so both modules produce matching label/tick colors across canvases.
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
        try {
            java.lang.reflect.Field f = org.jlab.groot.base.ColorPalette.class.getDeclaredField("colorPalette");
            f.setAccessible(true);
            @SuppressWarnings("unchecked")
            java.util.TreeMap<Integer, java.awt.Color> map = (java.util.TreeMap<Integer, java.awt.Color>) f.get(null);
            for (int i = 0; i < DET_COLORS.length; i++)
                map.put(DET_COLOR_BASE + i, DET_COLORS[i]);
        } catch (Exception e) {
            System.err.println("ECmodule: could not register custom palette colors: " + e.getMessage());
        }
    }

    private void addDetectorLabels(String canvasName) {
        final double[] detZ = { 0.0, 1.6, 2.5, 3.9, 5.4, 7.2, 7.3, 7.7, 8.0 };
        final String[] detLabel = { "Tgt", "HTCC M", "DC-R1", "DC-R2", "DC-R3",
                "FTOF 1b", "FTOF 1a", "PCAL", "EC" };

        EmbeddedCanvas canvas = this.getCanvas(canvasName);

        final int FONT = 11;
        final int LINE_H = 13;
        final int X_LEGEND = 110;
        final int Y_START = 60;

        // Pads 3-5: Vz vs r — draw ticks into the negative-r band and add legend
        for (int padIdx = 3; padIdx <= 5; padIdx++) {
            canvas.cd(padIdx);
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

    public ECmodule() {
        super(DetectorType.ECAL);
    }

    // -----------------------------------------------------------------------
    // Sub-detector helpers
    // -----------------------------------------------------------------------

    /** Returns 0/1/2 for PCAL/ECIN/ECOUT, or -1 for invalid layer. */
    private static int subDetIdx(int layer) {
        if (layer >= 1 && layer <= 3)
            return 0;
        if (layer >= 4 && layer <= 6)
            return 1;
        if (layer >= 7 && layer <= 9)
            return 2;
        return -1;
    }

    /** Local view number 1-3 (U/V/W) within the sub-detector. */
    private static int localView(int layer) {
        return ((layer - 1) % 3) + 1;
    }

    // -----------------------------------------------------------------------
    // Histogram definitions
    // -----------------------------------------------------------------------

    /**
     * 2D occupancy map: sector (x) vs layer 1-9 (y).
     * Colour = fraction of events with at least one hit in that cell [%].
     */
    public DataGroup sectorLayerMap() {
        DataGroup dg = new DataGroup(1, 1);
        H2F hi = histo2D("hi_sec_lay", "Sector",
                "Layer  (1-3=PCAL  4-6=ECIN  7-9=ECOUT)",
                NSECTORS, 0.5, NSECTORS + 0.5,
                NLAYERS, 0.5, NLAYERS + 0.5);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /**
     * Hit rate [kHz] vs sector — one pad per sub-detector.
     */
    public DataGroup rateBySector() {
        DataGroup dg = new DataGroup(3, 1);
        int[] colors = { 4, 2, 3 };
        for (int i = 0; i < 3; i++) {
            H1F hi = histo1D("hi_rate_sec_" + SUBDETS[i], SUBDETS[i],
                    "Sector", "Rate [kHz]",
                    NSECTORS, 0.5, NSECTORS + 0.5, 0);
            hi.setLineColor(colors[i]);
            hi.setLineWidth(2);
            dg.addDataSet(hi, i);
        }
        return dg;
    }

    /**
     * Per-sector 2D strip map for one sub-detector (6 pads, one per sector).
     * X = strip / component number. Y = view 1/2/3 = U/V/W.
     */
    public DataGroup sectorStripMap(int sdIdx) {
        DataGroup dg = new DataGroup(3, 2);
        for (int s = 1; s <= NSECTORS; s++) {
            H2F hi = histo2D("hi_strip_" + SUBDETS[sdIdx] + "_s" + s,
                    "Strip", "View (1=U  2=V  3=W)",
                    NSTRIPS, 0.5, NSTRIPS + 0.5,
                    3, 0.5, 3.5);
            hi.setTitle("Sector " + s);
            dg.addDataSet(hi, s - 1);
        }
        return dg;
    }

    /**
     * 2D position maps from MC::True — one pad per sub-detector.
     * Top row: R vs Z. Bottom row: Y vs X.
     */
    public DataGroup positionMap() {
        DataGroup dg = new DataGroup(3, 2);
        for (int i = 0; i < 3; i++) {
            H2F hiRZ = histo2D("hi_posRZ_" + SUBDETS[i],
                    SUBDETS[i],
                    "Z [m]", "R [m]",
                    200, EC_Z_MIN, EC_Z_MAX,
                    200, 0, EC_R_MAX);
            H2F hiXY = histo2D("hi_posXY_" + SUBDETS[i],
                    SUBDETS[i],
                    "X [m]", "Y [m]",
                    200, -EC_XY_MAX, EC_XY_MAX,
                    200, -EC_XY_MAX, EC_XY_MAX);
            dg.addDataSet(hiRZ, i);
            dg.addDataSet(hiXY, i + 3);
        }
        return dg;
    }

    /**
     * 2D rate map [kHz]: strip (x) vs view/U/V/W (y), averaged over sectors.
     * Matches the RGA reference plot layout (Figure 10 in RGA background note).
     * Strip counts: PCAL up to 68 (U), ECIN/ECOUT up to 36.
     */
    private static final int[] NSTRIPS_VIEW = { 68, 36, 36 }; // PCAL, ECIN, ECOUT

    public DataGroup viewRateMap() {
        DataGroup dg = new DataGroup(1, 3);
        for (int i = 0; i < 3; i++) {
            H2F hi = histo2D("hi_view_rate_" + SUBDETS[i], "Rate [kHz]",
                    "Strip", "View (1=U  2=V  3=W)",
                    NSTRIPS_VIEW[i], 0.5, NSTRIPS_VIEW[i] + 0.5,
                    3, 0.5, 3.5);
            hi.setTitle(SUBDETS[i]);
            dg.addDataSet(hi, i);
        }
        return dg;
    }

    /**
     * 2D rate map [kHz]: strip (x) vs sector (y) — one pad per sub-detector.
     * Mirrors the FTOF Rate Map layout but per-sub-detector.
     */
    public DataGroup rateMap() {
        DataGroup dg = new DataGroup(3, 1);
        for (int i = 0; i < 3; i++) {
            H2F hi = histo2D("hi_rate_map_" + SUBDETS[i], "Rate [kHz]",
                    "Strip", "Sector",
                    NSTRIPS, 0.5, NSTRIPS + 0.5,
                    NSECTORS, 0.5, NSECTORS + 0.5);
            hi.setTitle(SUBDETS[i]);
            dg.addDataSet(hi, i);
        }
        return dg;
    }

    /**
     * 2D PMT current map [µA]: same layout as Rate Map, weighted by
     * Edep[MeV] × yield[ph/MeV] × 1.6e-7 / (2h[cm]).
     * Placeholder yields — update from EC technical notes.
     */
    public DataGroup pmtCurrentMap() {
        DataGroup dg = new DataGroup(3, 1);
        for (int i = 0; i < 3; i++) {
            H2F hi = histo2D("hi_curr_map_" + SUBDETS[i], "PMT Current [#muA]",
                    "Strip", "Sector",
                    NSTRIPS, 0.5, NSTRIPS + 0.5,
                    NSECTORS, 0.5, NSECTORS + 0.5);
            hi.setTitle(SUBDETS[i]);
            dg.addDataSet(hi, i);
        }
        return dg;
    }

    /**
     * Per-sector 1D Vz background origin for one sub-detector.
     * Layout: DataGroup(2,3) — 6 pads, one per sector.
     * Multiple particle-species histograms overlaid per pad.
     */
    public DataGroup sectorBG(int sdIdx) {
        DataGroup dg = new DataGroup(2, 3);
        for (int is = 0; is < NSECTORS; is++) {
            int sector = is + 1;
            for (int ip = 0; ip < PNAMES.length; ip++) {
                H1F hi = histo1D(
                        "hi_bg_sec_" + SUBDETS[sdIdx] + "_s" + sector + "_" + PNAMES[ip],
                        SUBDETS[sdIdx] + " S" + sector + "-" + PNAMES[ip],
                        "Vz(m)", "Rate [kHz]",
                        200, ORIG_VZ_MIN, ORIG_VZ_MAX, 0);
                this.setHistoAttr(hi, ip < 5 ? ip + 1 : ip + 3);
                dg.addDataSet(hi, is);
            }
        }
        return dg;
    }

    /**
     * Origin of background: 2D vertex maps plus 1D Vz species distributions.
     * Layout (3 columns × 3 rows):
     * Row 0 (pads 0-2): Vx vs Vy per sub-detector
     * Row 1 (pads 3-5): Vz vs r per sub-detector
     * Row 2 (pads 6-8): 1D Vz by particle species per sub-detector
     */
    public DataGroup origin_bg() {
        DataGroup dg = new DataGroup(3, 3);
        for (int i = 0; i < 3; i++) {
            H2F hi_xy = histo2D("hi_bg_origin_xy_" + SUBDETS[i],
                    "Vx(m)", "Vy(m)",
                    200, -ORIG_DR, ORIG_DR,
                    200, -ORIG_DR, ORIG_DR);
            H2F hi_rz = histo2D("hi_bg_origin_rz_" + SUBDETS[i],
                    "Vz(m)", "r(m)",
                    200, ORIG_VZ_MIN, ORIG_VZ_MAX,
                    200, ORIG_R_MIN, ORIG_DR);
            dg.addDataSet(hi_xy, 0 + i);
            dg.addDataSet(hi_rz, 3 + i);

            for (int ip = 0; ip < PNAMES.length; ip++) {
                H1F hi_bg = histo1D("hi_bg_" + SUBDETS[i] + "_" + PNAMES[ip],
                        PNAMES[ip], "Vz(m)", "Rate [kHz]",
                        200, ORIG_VZ_MIN, ORIG_VZ_MAX, 0);
                this.setHistoAttr(hi_bg, ip < 5 ? ip + 1 : ip + 3);
                dg.addDataSet(hi_bg, 6 + i);
            }
        }
        return dg;
    }

    // -----------------------------------------------------------------------
    // Fill methods
    // -----------------------------------------------------------------------

    public void fillPositionMap(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            if (h.getTrue() == null)
                continue;
            int isd = subDetIdx(h.getLayer());
            if (isd < 0)
                continue;
            double x = h.getTrue().getPosition().x() / 1000.0;
            double y = h.getTrue().getPosition().y() / 1000.0;
            double z = h.getTrue().getPosition().z() / 1000.0;
            double r = Math.sqrt(x * x + y * y);
            dg.getH2F("hi_posRZ_" + SUBDETS[isd]).fill(z, r);
            dg.getH2F("hi_posXY_" + SUBDETS[isd]).fill(x, y);
        }
    }

    public void fillRateMap(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            if (h.getTrue() == null || h.getTrue().getEdep() <= EDEP_THRESHOLD)
                continue;
            int isd = subDetIdx(h.getLayer());
            if (isd < 0)
                continue;
            dg.getH2F("hi_rate_map_" + SUBDETS[isd]).fill(h.getComponent(), h.getSector());
        }
    }

    public void fillViewRateMap(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            if (h.getTrue() == null || h.getTrue().getEdep() <= EDEP_THRESHOLD)
                continue;
            int isd = subDetIdx(h.getLayer());
            if (isd < 0)
                continue;
            dg.getH2F("hi_view_rate_" + SUBDETS[isd])
                    .fill(h.getComponent(), localView(h.getLayer()));
        }
    }

    public void fillCurrentMap(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            if (h.getTrue() == null)
                continue;
            double edep = h.getTrue().getEdep();
            if (edep <= EDEP_THRESHOLD)
                continue;
            int isd = subDetIdx(h.getLayer());
            if (isd < 0)
                continue;
            double current = edep * YIELD[isd] * 1.6e-7 / TWO_H[isd];
            dg.getH2F("hi_curr_map_" + SUBDETS[isd]).fill(h.getComponent(), h.getSector(), current);
        }
    }

    public void fillSectorBG(DataGroup dg, List<Hit> hits, int targetSdIdx) {
        for (Hit h : hits) {
            if (h.getTrue() == null)
                continue;
            int isd = subDetIdx(h.getLayer());
            if (isd != targetSdIdx)
                continue;
            int sector = h.getSector();
            if (sector < 1 || sector > NSECTORS)
                continue;
            True t = h.getTrue();
            double vz = t.getVertex().z() / 1000.0;
            dg.getH1F("hi_bg_sec_" + SUBDETS[isd] + "_s" + sector + "_all").fill(vz);
            String pname = this.pidToName(Math.abs(t.getPid()));
            if (pname != null)
                dg.getH1F("hi_bg_sec_" + SUBDETS[isd] + "_s" + sector + "_" + pname).fill(vz);
            else
                dg.getH1F("hi_bg_sec_" + SUBDETS[isd] + "_s" + sector + "_other").fill(vz);
        }
    }

    public void fillOrigin(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            if (h.getTrue() == null)
                continue;
            True t = h.getTrue();
            int isd = subDetIdx(h.getLayer());
            if (isd < 0)
                continue;

            double vx = t.getVertex().x() / 1000.0;
            double vy = t.getVertex().y() / 1000.0;
            double vz = t.getVertex().z() / 1000.0;
            double r = Math.sqrt(vx * vx + vy * vy);

            dg.getH2F("hi_bg_origin_xy_" + SUBDETS[isd]).fill(vx, vy);
            dg.getH2F("hi_bg_origin_rz_" + SUBDETS[isd]).fill(vz, r);

            dg.getH1F("hi_bg_" + SUBDETS[isd] + "_all").fill(vz);
            String pname = this.pidToName(Math.abs(t.getPid()));
            if (pname != null)
                dg.getH1F("hi_bg_" + SUBDETS[isd] + "_" + pname).fill(vz);
            else
                dg.getH1F("hi_bg_" + SUBDETS[isd] + "_other").fill(vz);
        }
    }

    // -----------------------------------------------------------------------
    // Module interface
    // -----------------------------------------------------------------------

    @Override
    public void createHistos() {
        this.getHistos().put("Sector-Layer Map", this.sectorLayerMap());
        this.getHistos().put("Rate by Sector", this.rateBySector());
        this.getHistos().put("PCAL Strip Map", this.sectorStripMap(0));
        this.getHistos().put("ECIN Strip Map", this.sectorStripMap(1));
        this.getHistos().put("ECOUT Strip Map", this.sectorStripMap(2));
        this.getHistos().put("Position Map", this.positionMap());
        this.getHistos().put("Rate Map", this.rateMap());
        this.getHistos().put("View Rate Map", this.viewRateMap());
        this.getHistos().put("PMT Current Map", this.pmtCurrentMap());
        this.getHistos().put("Origin of Bg", this.origin_bg());
        for (int i = 0; i < SUBDETS.length; i++)
            this.getHistos().put("Origin of Bg - Sector " + SUBDETS[i], this.sectorBG(i));
    }

    @Override
    public void fillHistos(Event event) {
        List<Hit> hits = event.getHits(DetectorType.ECAL);
        if (hits == null)
            return;

        DataGroup dgMap = this.getHistos().get("Sector-Layer Map");
        DataGroup dgSec = this.getHistos().get("Rate by Sector");

        for (Hit h : hits) {
            if (h.getTrue() == null || h.getTrue().getEdep() <= EDEP_THRESHOLD)
                continue;
            int sector = h.getSector();
            int layer = h.getLayer();
            int comp = h.getComponent();
            int isd = subDetIdx(layer);
            if (isd < 0)
                continue;

            dgMap.getH2F("hi_sec_lay").fill(sector, layer);
            dgSec.getH1F("hi_rate_sec_" + SUBDETS[isd]).fill(sector);

            DataGroup dgStrip = this.getHistos().get(SUBDETS[isd] + " Strip Map");
            dgStrip.getH2F("hi_strip_" + SUBDETS[isd] + "_s" + sector)
                    .fill(comp, localView(layer));
        }

        this.fillPositionMap(this.getHistos().get("Position Map"), hits);
        this.fillRateMap(this.getHistos().get("Rate Map"), hits);
        this.fillViewRateMap(this.getHistos().get("View Rate Map"), hits);
        this.fillCurrentMap(this.getHistos().get("PMT Current Map"), hits);
        this.fillOrigin(this.getHistos().get("Origin of Bg"), hits);
        for (int i = 0; i < SUBDETS.length; i++)
            this.fillSectorBG(this.getHistos().get("Origin of Bg - Sector " + SUBDETS[i]), hits, i);
    }

    @Override
    public void analyzeHistos() {
        double x = Double.parseDouble(System.getProperty("lumi", "117000.0"));
        double lumiScale = x / 117000.0;

        // Sector-Layer Map → occupancy [%]
        this.normalize(this.getHistos().get("Sector-Layer Map"), 0.01 * this.getNevents());
        this.normalize(this.getHistos().get("Sector-Layer Map"), lumiScale);

        // Rate by Sector → kHz, lumi-scaled
        this.normalizeToTime(this.getHistos().get("Rate by Sector"));
        this.normalize(this.getHistos().get("Rate by Sector"), lumiScale);

        // Strip maps → occupancy [%]
        for (String sd : SUBDETS) {
            this.normalize(this.getHistos().get(sd + " Strip Map"), 0.01 * this.getNevents());
            this.normalize(this.getHistos().get(sd + " Strip Map"), lumiScale);
        }

        // Position Map → rate [kHz], lumi-scaled
        this.normalizeToTime(this.getHistos().get("Position Map"));
        this.normalize(this.getHistos().get("Position Map"), lumiScale);

        // Rate Map → kHz (1 MeV threshold keeps max ~140 kHz, well below GROOT's ×10³
        // prefix)
        this.normalizeToTime(this.getHistos().get("Rate Map")); // kHz
        this.normalize(this.getHistos().get("Rate Map"), lumiScale);

        // View Rate Map → kHz per sector (matches RGA reference plot)
        // Divide by nActiveSectors to convert total-over-sectors → per-sector average
        int nActive = NSECTORS;
        this.normalizeToTime(this.getHistos().get("View Rate Map")); // kHz total
        this.normalize(this.getHistos().get("View Rate Map"), (double) nActive); // → per sector
        this.normalize(this.getHistos().get("View Rate Map"), lumiScale);

        // PMT Current Map → µA (normalizeToTime with units=1 → Hz, so weight × Hz = µA)
        for (String sd : SUBDETS)
            this.normalizeToTime(this.getHistos().get("PMT Current Map").getH2F("hi_curr_map_" + sd), 1);
        this.normalize(this.getHistos().get("PMT Current Map"), lumiScale);

        // Origin of Bg → kHz, lumi-scaled
        this.normalizeToTime(this.getHistos().get("Origin of Bg"));
        this.normalize(this.getHistos().get("Origin of Bg"), lumiScale);

        // Origin of Bg - Sector → kHz, lumi-scaled
        for (String sd : SUBDETS) {
            this.normalizeToTime(this.getHistos().get("Origin of Bg - Sector " + sd));
            this.normalize(this.getHistos().get("Origin of Bg - Sector " + sd), lumiScale);
        }

        System.out.println("=====================================================\n");

    }

    @Override
    public void setPlottingOptions(String name) {
        GStyle.getAxisAttributesX().setTitleFontSize(18);
        GStyle.getAxisAttributesY().setTitleFontSize(18);
        GStyle.getAxisAttributesX().setLabelFontSize(18);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        GStyle.getAxisAttributesZ().setTitleFontSize(18);
        GStyle.getAxisAttributesZ().setLabelFontSize(18);

        for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
            pad.setTitle("");

        if (name.equals("Sector-Layer Map")) {
            this.setupColorAxis(name, "Occupancy [%]");
            // this.setLogZ(name);
        }

        if (name.contains("Strip Map")) {
            this.setupColorAxis(name, "Occupancy [%]");
            // this.setLogZ(name);
            // ECIN and ECOUT only have 36 strips — clip X to 40
            if (name.startsWith("ECIN") || name.startsWith("ECOUT")) {
                for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
                    pad.getAxisX().setRange(0.5, 40.5);
            }
        }

        if (name.equals("Position Map")) {
            this.setupColorAxis(name, "Rate [kHz]");
            // this.setLogZ(name);
            // // All 6 pads (R vs Z top row + Y vs X bottom row) share same Z range.
            // for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
            // pad.getAxisZ().setRange(1, 10000);
        }

        if (name.equals("Rate by Sector"))
            this.setLegend(name, 250, 50);

        if (name.equals("View Rate Map")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
            this.getCanvas(name).setPreferredSize(new java.awt.Dimension(700, 1100));
            // for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
            // pad.getAxisY().setRange(0.5, 3.5);
            // pad.getAxisZ().setRange(0, 500);
            // }
        }

        if (name.equals("Rate Map")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
            this.getCanvas(name).setPreferredSize(new java.awt.Dimension(1600, 550));
            // double[] zMaxKHz = { 900, 700, 500 }; // PCAL, ECIN, ECOUT
            java.util.List<EmbeddedPad> pads = this.getCanvas(name).getCanvasPads();
            for (int i = 0; i < pads.size(); i++) {
                // pads.get(i).getAxisZ().setRange(0, zMaxKHz[i]);
                if (i > 0) // ECIN and ECOUT: only 36 strips, show up to 40
                    pads.get(i).getAxisX().setRange(0.5, 40.5);
            }
        }

        if (name.equals("PMT Current Map")) {
            this.setupColorAxis(name, "Current [#muA]");
            this.setWideColorBar(name);
            java.util.List<EmbeddedPad> pads = this.getCanvas(name).getCanvasPads();
            for (int i = 1; i < pads.size(); i++) // ECIN (1) and ECOUT (2): 36 strips
                pads.get(i).getAxisX().setRange(0.5, 40.5);
        }

        if (name.equals("Origin of Bg")) {
            this.setupColorAxis(name, "Rate [kHz]");
            this.setWideColorBar(name);
            this.setLogZ(name);
            this.setLegend(name, 200, 70);
            // Fix Y max on the 1D Vz pads (row 2, pad indices 6-8).
            java.util.List<EmbeddedPad> pads = this.getCanvas(name).getCanvasPads();
            // double[] vzYMax = { 200000, 30000, 30000 }; // PCAL, ECIN, ECOUT
            for (int i = 6; i < 9 && i < pads.size(); i++)
                // pads.get(i).getAxisY().setRange(0, vzYMax[i - 6]);
                this.addDetectorLabels(name);
        }

        if (name.startsWith("Origin of Bg - Sector"))
            this.setLegend(name, 450, 150);
    }
}
