package modules;

import java.util.List;
import java.util.stream.Collectors;
import objects.Hit;
import analysis.Module;
import objects.Event;
import objects.True;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.base.PadMargins;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.LatexText;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;

/**
 * Background module for the HTCC (High Threshold Cherenkov Counter).
 *
 * Studies the impact of adding wing shielding in front of the HTCC by
 * measuring hit rates, NPhe spectra, and background particle origins.
 *
 * HTCC bank column mapping (from htcc_hitprocess.cpp integrateDgt):
 * bank "layer" = idhalf → half (1 or 2, left/right PMT within sector-ring)
 * bank "component" = idring → ring (1–4, polar angle layer)
 * A unique PMT = (sector, half, ring), giving 6×2×4 = 48 PMTs total.
 *
 * ADC → NPhe conversion:
 * GEMC stores ADC = 100 × NPhe × mc_gain.
 * mc_gain = 0.5 for all channels (from /calibration/htcc/mc_gain CCDB table).
 * Therefore: NPhe = ADC / (100 × 0.5) = ADC / 50.
 *
 * For optical-photon simulation MC::True pid=0 (optical photon);
 * the module falls back to mpid to identify the Cherenkov-producing particle.
 */
public class HTCCmodule extends Module {

    private static final int NSECTORS = 6;
    private static final int NRINGS = 4; // component 1–4
    private static final int NHALVES = 2; // layer 1–2
    private static final double ADC_TO_NPHE = 50.0; // 100 * mc_gain(=0.5)

    private static final double PMT_GAIN = Double.parseDouble(System.getProperty("htcc.gain", "1.0e6"));
    private static final double ELECTRON_C = 1.6e-19; // coulombs

    // Raw hit counts per ring bin (index 0–5 maps to bins 0.5–5.5)
    // kept separately so we can compute Poisson errors after normalization
    private final long[] rawNphe0 = new long[6];
    private final long[] rawNphe10 = new long[6];

    // Total photoelectrons accumulated per ring (1-indexed, [0] unused)
    // Used to compute anode current: I = R_NPhe × e × G
    private final long[] totalNpheByRing = new long[NRINGS + 1];

    // Hit count per ring (1-indexed) for mean-NPhe-per-hit calculation
    private final long[] countHitsByRing = new long[NRINGS + 1];

    public HTCCmodule() {
        super(DetectorType.HTCC);
    }

    ///////////// debug helper //////
    private void printHTCCSummary(double lumiScale,
            int activeSectors,
            int nAllPmts,
            int nPerRing) {

        String config = System.getenv().getOrDefault("HTCC_CONFIG", "unknown");
        String deltaTns = System.getenv().getOrDefault("HTCC_DT_NS", "252.0");
        String outFile = System.getenv().getOrDefault("HTCC_CSV", "htcc_rate_probability.csv");

        System.out.println("\n================ HTCC MODULE SUMMARY ================");
        System.out.printf("Config name              = %s%n", config);
        System.out.printf("Number of events         = %d%n", this.getNevents());
        System.out.printf("Luminosity scale         = %.6f%n", lumiScale);
        System.out.printf("Probability time window  = %s ns%n", deltaTns);

        System.out.printf("Active sectors           = %d%n", activeSectors);
        System.out.printf("HTCC rings               = %d%n", NRINGS);
        System.out.printf("HTCC halves per ring     = %d%n", NHALVES);
        System.out.printf("Total PMTs used          = %d%n", nAllPmts);
        System.out.printf("PMTs per ring used       = %d%n", nPerRing);
        System.out.printf("ADC to NPhe factor       = %.3f%n", ADC_TO_NPHE);
        System.out.printf("CSV output file          = %s%n", outFile);
        System.out.println("Probability formula      = 100 * (1 - exp(-rate_Hz * deltaT_s))");
        System.out.println("=====================================================\n");
    }

    // ////////// For csv file /////////
    // private void writeRateProbabilityCSV(double lumiScale,
    // double[] npheRateKHz,
    // double[] anodeCurrentUa) {

    // String config = System.getenv().getOrDefault("HTCC_CONFIG", "unknown");
    // double deltaTns =
    // Double.parseDouble(System.getenv().getOrDefault("HTCC_DT_NS", "252.0"));
    // String outFile = System.getenv().getOrDefault("HTCC_CSV",
    // "htcc_rate_probability_no_sec4.csv");
    // double deltaTs = deltaTns * 1.0e-9;

    // System.out.printf("%n>>>>> timeWindow_ns = %.3f ns%n", deltaTns);

    // DataGroup dg = this.getHistos().get("Rate vs Ring (Threshold)");
    // H1F h0 = dg.getH1F("hi_rvr_nphe0");
    // H1F h10 = dg.getH1F("hi_rvr_nphe3");

    // boolean append = true;
    // boolean writeHeader = !(new java.io.File(outFile).exists());

    // try (PrintWriter out = new PrintWriter(new FileWriter(outFile, append))) {

    // if (writeHeader) {
    // out.println("config,lumi_scale,delta_t_ns,ring,"
    // + "threshold,rate_kHz_per_PMT,rate_Hz_per_PMT,probability_percent,"
    // + "nphe_rate_kHz_per_PMT,anode_current_uA");
    // }

    // for (int ring = 1; ring <= NRINGS; ring++) {

    // int bin = h0.getAxis().getBin(ring + 0.5);

    // double rate0_kHz = h0.getBinContent(bin);
    // double rate3_kHz = h10.getBinContent(bin);
    // double rate0_Hz = rate0_kHz * 1000.0;
    // double rate3_Hz = rate3_kHz * 1000.0;
    // double prob0 = 100.0 * (1.0 - Math.exp(-rate0_Hz * deltaTs));
    // double prob3 = 100.0 * (1.0 - Math.exp(-rate3_Hz * deltaTs));

    // // NPhe rate and current are per-ring (not per-threshold), written once for
    // // NPhe>0 row
    // out.printf("%s,%.6f,%.3f,%d,NPhe>0,%.6f,%.6f,%.6f,%.4f,%.6f%n",
    // config, lumiScale, deltaTns, ring,
    // rate0_kHz, rate0_Hz, prob0,
    // npheRateKHz[ring], anodeCurrentUa[ring]);

    // out.printf("%s,%.6f,%.3f,%d,NPhe>3,%.6f,%.6f,%.6f,%.4f,%.6f%n",
    // config, lumiScale, deltaTns, ring,
    // rate3_kHz, rate3_Hz, prob3,
    // npheRateKHz[ring], anodeCurrentUa[ring]);
    // }

    // } catch (IOException e) {
    // System.err.println("[HTCCmodule] ERROR writing CSV file: " + outFile);
    // e.printStackTrace();
    // }
    // }

    // -----------------------------------------------------------------------
    // Histogram definitions
    // -----------------------------------------------------------------------

    /** 2D hit-rate map: sector (x) vs ring (y), in % occupancy per PMT pair. */
    public DataGroup occupancy() {
        DataGroup dg = new DataGroup(1, 1);
        H2F hi = histo2D("hi_occ_htcc", "Sector", "Ring",
                NSECTORS, 0.5, NSECTORS + 0.5,
                NRINGS, 0.5, NRINGS + 0.5);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /**
     * Hit rate [kHz] vs sector, Stacked by particle species (uses mpid for
     * optical-photon hits).
     */
    public DataGroup rateBySector() {
        DataGroup dg = new DataGroup(1, 1);
        for (int ip = 0; ip < PNAMES.length; ip++) {
            H1F hi = histo1D("hi_rate_sector_" + PNAMES[ip], PNAMES[ip],
                    "Sector", "Rate [kHz]",
                    NSECTORS, 0.5, NSECTORS + 0.5, 0);
            this.setHistoAttr(hi, ip < 5 ? ip + 1 : ip + 3);
            dg.addDataSet(hi, 0);
        }
        return dg;
    }

    /**
     * Hit rate [kHz] vs ring — shows which polar-angle ring is most affected
     * by the wing shadow (ring 1 is most forward, closest to Moller cone).
     * Ring = bank "component" (1–4); bank "layer" is the half (1–2).
     */
    public DataGroup rateByRing() {
        DataGroup dg = new DataGroup(1, 1);
        for (int ip = 0; ip < PNAMES.length; ip++) {
            H1F hi = histo1D("hi_rate_ring_" + PNAMES[ip], PNAMES[ip],
                    "Ring", "Rate [kHz]",
                    NRINGS, 0.5, NRINGS + 0.5, 0);
            this.setHistoAttr(hi, ip < 5 ? ip + 1 : ip + 3);
            dg.addDataSet(hi, 0);
        }
        return dg;
    }

    /**
     * NPhe spectrum [kHz] — distribution of photoelectrons per HTCC hit.
     * NPhe = ADC / ADC_TO_NPHE = ADC / 50.
     * Background hits from Moller electrons typically give NPhe ~ 5–30.
     */
    public DataGroup npheSpectrum() {
        DataGroup dg = new DataGroup(1, 1);
        H1F hi = histo1D("hi_nphe", "HTCC NPhe", "NPhe", "Rate [kHz]",
                100, 0, 50, 4);
        hi.setLineWidth(2);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /**
     * NPhe spectrum per PMT [kHz/PMT], log-Y — mirrors Figure 6 left:
     * black = all 48 PMTs averaged, then rings 1–4 individually.
     * Ring is identified by bank "component" (1–4).
     */
    public DataGroup nphePerPMT() {
        DataGroup dg = new DataGroup(1, 1);
        H1F hi_all = histo1D("hi_nphe_pmt_all", "All PMTs", "NPhe", "Rate [kHz/PMT]", 100, 0, 100, 0);
        H1F hi_ring1 = histo1D("hi_nphe_pmt_ring1", "Ring 1", "NPhe", "Rate [kHz/PMT]", 100, 0, 100, 0);
        H1F hi_ring2 = histo1D("hi_nphe_pmt_ring2", "Ring 2", "NPhe", "Rate [kHz/PMT]", 100, 0, 100, 0);
        H1F hi_ring3 = histo1D("hi_nphe_pmt_ring3", "Ring 3", "NPhe", "Rate [kHz/PMT]", 100, 0, 100, 0);
        H1F hi_ring4 = histo1D("hi_nphe_pmt_ring4", "Ring 4", "NPhe", "Rate [kHz/PMT]", 100, 0, 100, 0);
        hi_all.setLineColor(1);
        hi_all.setLineWidth(2);
        hi_ring1.setLineColor(4);
        hi_ring1.setLineWidth(2);
        hi_ring2.setLineColor(3);
        hi_ring2.setLineWidth(2);
        hi_ring3.setLineColor(44);
        hi_ring3.setLineWidth(2);
        hi_ring4.setLineColor(2);
        hi_ring4.setLineWidth(2);
        dg.addDataSet(hi_all, 0);
        dg.addDataSet(hi_ring1, 0);
        dg.addDataSet(hi_ring2, 0);
        dg.addDataSet(hi_ring3, 0);
        dg.addDataSet(hi_ring4, 0);
        return dg;
    }

    /**
     * PMT rate [kHz/PMT] vs ring — mirrors Figure 6 right:
     * blue = NPhe > 0 (all hits), red = NPhe > 3 (quality threshold).
     * NPhe computed as ADC / 50 before threshold comparison.
     * hi_10 added last so it draws on top; hi_0 (wider line) stays visible
     * underneath.
     */
    public DataGroup rateVsRingThreshold() {
        DataGroup dg = new DataGroup(1, 1);
        H1F hi0 = histo1D(
                "hi_rvr_nphe0",
                "NPhe > 0",
                "Ring",
                "Rate [kHz/PMT]",
                6, 0.0, 6.0, 0);

        H1F hi10 = histo1D(
                "hi_rvr_nphe3",
                "NPhe > 3",
                "Ring",
                "Rate [kHz/PMT]",
                6, 0.0, 6.0, 0);

        hi0.setLineColor(4);
        hi0.setLineWidth(3);
        hi10.setLineColor(2);
        hi10.setLineWidth(2);
        dg.addDataSet(hi0, 0);
        dg.addDataSet(hi10, 0);
        return dg;
    }

    /**
     * Mean NPhe per hit vs ring — placeholder H1F replaced by GraphErrors in
     * analyzeHistos after accumulation. No time/lumi normalization needed
     * because it is a ratio (sum_nphe / n_hits).
     */
    public DataGroup meanNpheByRing() {
        DataGroup dg = new DataGroup(1, 1);
        H1F hi = histo1D("hi_mean_nphe_placeholder", "Mean NPhe", "Ring", "<NPhe> / hit",
                NRINGS, 0.5, NRINGS + 0.5, 0);
        dg.addDataSet(hi, 0);
        return dg;
    }

    /**
     * Background-origin plots:
     * pad 0 — 2D vertex (Vz, Vr)
     * pad 1 — 2D vertex (Vx, Vy)
     * pad 2 — Vz stacked by species [kHz]
     * pad 3 — kinetic energy stacked by species [kHz]
     */
    public DataGroup originBG() {
        DataGroup dg = new DataGroup(2, 2);

        H2F hi_rz = histo2D("hi_bg_rz", "Vz (mm)", "Vr (mm)",
                200, -500, 1500, 200, 0, 700);
        H2F hi_xy = histo2D("hi_bg_xy", "Vx (mm)", "Vy (mm)",
                200, -700, 700, 200, -700, 700);
        dg.addDataSet(hi_rz, 0);
        dg.addDataSet(hi_xy, 1);

        for (int ip = 0; ip < PNAMES.length; ip++) {
            H1F hi_vz = histo1D("hi_bg_vz_" + PNAMES[ip], PNAMES[ip],
                    "Vz (mm)", "Rate [kHz]",
                    200, -500, 1500, 0);
            this.setHistoAttr(hi_vz, ip < 5 ? ip + 1 : ip + 3);
            dg.addDataSet(hi_vz, 2);
        }

        for (int ip = 0; ip < PNAMES.length; ip++) {
            H1F hi_e = histo1D("hi_bg_energy_" + PNAMES[ip], PNAMES[ip],
                    "E (MeV)", "Rate [kHz]",
                    200, 0, 500, 0);
            this.setHistoAttr(hi_e, ip < 5 ? ip + 1 : ip + 3);
            dg.addDataSet(hi_e, 3);
        }

        return dg;
    }

    /**
     * Per-sector background vertex Vz [kHz].
     * Comparing sector 4 (with wings) to sectors 1–3, 5–6 (without wings)
     * is the main diagnostic for the shielding impact.
     */
    public DataGroup sectorBG() {
        DataGroup dg = new DataGroup(3, 2);
        for (int is = 0; is < NSECTORS; is++) {
            int s = is + 1;
            for (int ip = 0; ip < PNAMES.length; ip++) {
                H1F hi = histo1D("hi_bg_s" + s + "_" + PNAMES[ip], PNAMES[ip],
                        "Vz (mm)", "Rate [kHz]",
                        200, -500, 1500, 0);
                this.setHistoAttr(hi, ip < 5 ? ip + 1 : ip + 3);
                dg.addDataSet(hi, is);
            }
        }
        return dg;
    }

    // -----------------------------------------------------------------------
    // Module interface
    // -----------------------------------------------------------------------

    @Override
    public void createHistos() {
        this.getHistos().put("Occupancy", this.occupancy());
        this.getHistos().put("Rate by Sector", this.rateBySector());
        this.getHistos().put("Rate by Ring", this.rateByRing());
        this.getHistos().put("NPhe Spectrum", this.npheSpectrum());
        this.getHistos().put("NPhe per PMT", this.nphePerPMT());
        this.getHistos().put("Rate vs Ring (Threshold)", this.rateVsRingThreshold());
        this.getHistos().put("Mean NPhe by Ring", this.meanNpheByRing());
        this.getHistos().put("Origin of Bg", this.originBG());
        this.getHistos().put("Sector BG", this.sectorBG());
    }

    @Override
    public void fillHistos(Event event) {
        List<Hit> hits = event.getHits(DetectorType.HTCC);
        if (hits == null)
            return;
        fillOccupancy(this.getHistos().get("Occupancy"), hits);
        fillRateBySector(this.getHistos().get("Rate by Sector"), hits);
        fillRateByRing(this.getHistos().get("Rate by Ring"), hits);
        fillNPhe(this.getHistos().get("NPhe Spectrum"), hits);
        fillNphePerPMT(this.getHistos().get("NPhe per PMT"), hits);
        fillRateVsRingThreshold(this.getHistos().get("Rate vs Ring (Threshold)"), hits);
        fillOrigin(this.getHistos().get("Origin of Bg"), hits);
        fillSectorBG(this.getHistos().get("Sector BG"), hits);
    }

    @Override
    public void analyzeHistos() {
        final int activeSectors = NSECTORS;
        final int N_ALL_PMTS = activeSectors * NRINGS * NHALVES;
        final int N_PER_RING = activeSectors * NHALVES;

        // Luminosity scaling:
        double x = Double.parseDouble(System.getProperty("lumi", "117000.0"));
        double lumiScale = x / 117000.0;

        printHTCCSummary(lumiScale, activeSectors, N_ALL_PMTS, N_PER_RING);

        // Occupancy [%] per PMT pair: hits / (0.01 × nevents × NHALVES)
        this.normalize(this.getHistos().get("Occupancy"),
                0.01 * this.getNevents() * NHALVES);
        this.normalize(this.getHistos().get("Occupancy"), lumiScale);

        this.normalizeToTime(this.getHistos().get("Rate by Sector"));
        this.normalize(this.getHistos().get("Rate by Sector"), lumiScale);

        this.normalizeToTime(this.getHistos().get("Rate by Ring"));
        this.normalize(this.getHistos().get("Rate by Ring"), lumiScale);

        this.normalizeToTime(this.getHistos().get("NPhe Spectrum"));
        this.normalize(this.getHistos().get("NPhe Spectrum"), lumiScale);

        this.normalizeToTime(this.getHistos().get("Origin of Bg"));
        this.normalize(this.getHistos().get("Origin of Bg"), lumiScale);

        this.normalizeToTime(this.getHistos().get("Sector BG"));
        this.normalize(this.getHistos().get("Sector BG"), lumiScale);

        // NPhe per PMT: kHz, lumi-scaled, then divide by number of PMTs
        DataGroup npheGroup = this.getHistos().get("NPhe per PMT");
        this.normalizeToTime(npheGroup);
        this.normalize(npheGroup, lumiScale);
        this.normalize(npheGroup.getH1F("hi_nphe_pmt_all"), (double) N_ALL_PMTS);
        for (int ir = 1; ir <= NRINGS; ir++)
            this.normalize(npheGroup.getH1F("hi_nphe_pmt_ring" + ir), (double) N_PER_RING);

        // Rate vs ring by threshold: kHz, lumi-scaled, then per-PMT (12 per ring)
        DataGroup rvrGroup = this.getHistos().get("Rate vs Ring (Threshold)");
        this.normalizeToTime(rvrGroup);
        this.normalize(rvrGroup, lumiScale);
        this.normalize(rvrGroup, (double) N_PER_RING);

        // Set Poisson errors: relative error = 1/sqrt(N_raw), so |err| =
        // rate/sqrt(N_raw)
        H1F h0 = rvrGroup.getH1F("hi_rvr_nphe0");
        H1F h10 = rvrGroup.getH1F("hi_rvr_nphe3");
        for (int ir = 0; ir < NRINGS; ir++) {
            int bin = h0.getAxis().getBin(ir + 1.5);
            if (rawNphe0[ir] > 0)
                h0.setBinError(bin, h0.getBinContent(bin) / Math.sqrt(rawNphe0[ir]));
            if (rawNphe10[ir] > 0)
                h10.setBinError(bin, h10.getBinContent(bin) / Math.sqrt(rawNphe10[ir]));
        }

        // Compute total NPhe rate per PMT per ring → anode current estimate
        double deltaTns = Double.parseDouble(System.getenv().getOrDefault("HTCC_DT_NS", "252.0"));
        double deltaTs = deltaTns * 1.0e-9;
        double[] npheRateKHz = new double[NRINGS + 1]; // [kHz/PMT]
        double[] anodeCurrentUa = new double[NRINGS + 1]; // [μA]
        for (int ir = 1; ir <= NRINGS; ir++) {
            // kHz/PMT = total_NPhe × lumiScale / (N_events × T_s × N_PMTs/ring × 1000)
            npheRateKHz[ir] = totalNpheByRing[ir] * lumiScale
                    / (this.getNevents() * deltaTs * N_PER_RING * 1000.0);
            anodeCurrentUa[ir] = npheRateKHz[ir] * 1000.0 * ELECTRON_C * PMT_GAIN * 1.0e6;
        }

        System.out.println("\n--- HTCC Anode Current Estimate (G=" + (long) PMT_GAIN + ") ---");
        System.out.printf("%-8s  %-18s  %-18s%n", "Ring", "NPhe Rate [kHz/PMT]", "I_anode [μA]");
        for (int ir = 1; ir <= NRINGS; ir++)
            System.out.printf("%-8d  %-18.2f  %-18.4f%n", ir, npheRateKHz[ir], anodeCurrentUa[ir]);
        System.out.println("----------------------------------------------------");

        // Write rate and probability table after all normalizations
        // writeRateProbabilityCSV(lumiScale, npheRateKHz, anodeCurrentUa);

        // Mean NPhe per hit vs ring (ratio: no time/lumi normalization needed)
        GraphErrors gMean = new GraphErrors("Mean NPhe");
        gMean.setTitleX("Ring");
        gMean.setTitleY("<NPhe> / hit");
        gMean.setMarkerColor(4);
        gMean.setMarkerStyle(20);
        gMean.setMarkerSize(8);
        gMean.setLineColor(4);
        for (int ir = 1; ir <= NRINGS; ir++) {
            if (countHitsByRing[ir] > 0) {
                double mean = (double) totalNpheByRing[ir] / countHitsByRing[ir];
                double err = mean / Math.sqrt(countHitsByRing[ir]);
                gMean.addPoint(ir, mean, 0.4, err);
            }
        }
        DataGroup meanGroup = new DataGroup(1, 1);
        meanGroup.addDataSet(gMean, 0);
        this.getHistos().put("Mean NPhe by Ring", meanGroup);

        // Replace the H1F DataGroup with GraphErrors BEFORE drawHistos() runs.
        // GROOT H1F silently ignores the "PE" draw option; GraphErrors always
        // renders as PE-style markers+error bars, so swapping here is the only
        // reliable way to get the ROOT-style PE plot.
        GraphErrors g0 = new GraphErrors("NPhe > 0");
        GraphErrors g10 = new GraphErrors("NPhe > 3");
        g0.setTitle("NPhe > 0");
        g10.setTitle("NPhe > 3");
        // Use bin centers + half-bin-width x error so each point spans its full bin
        for (int ib = 0; ib < h0.getAxis().getNBins(); ib++) {
            double xc = h0.getAxis().getBinCenter(ib);
            double ex = h0.getAxis().getBinWidth(ib) / 2.0;
            double y0 = h0.getBinContent(ib);
            double y10 = h10.getBinContent(ib);
            if (y0 > 0)
                g0.addPoint(xc, y0, ex, h0.getBinError(ib));
            if (y10 > 0)
                g10.addPoint(xc, y10, ex, h10.getBinError(ib));
        }
        g0.setTitleX("Ring");
        g0.setTitleY("Rate [kHz/PMT]");
        g0.setMarkerColor(4);
        g0.setMarkerStyle(20);
        g0.setMarkerSize(10);
        g0.setLineColor(4);
        g0.setLineThickness(3);
        g10.setMarkerColor(2);
        g10.setMarkerStyle(20);
        g10.setMarkerSize(10);
        g10.setLineColor(2);
        g10.setLineThickness(3);
        DataGroup rvrGroupPE = new DataGroup(1, 1);
        rvrGroupPE.addDataSet(g0, 0);
        rvrGroupPE.addDataSet(g10, 0);
        this.getHistos().put("Rate vs Ring (Threshold)", rvrGroupPE);
    }

    @Override
    public void setPlottingOptions(String name) {
        GStyle.getAxisAttributesX().setTitleFontSize(52);
        GStyle.getAxisAttributesY().setTitleFontSize(52);
        GStyle.getAxisAttributesX().setLabelFontSize(48);
        GStyle.getAxisAttributesY().setLabelFontSize(48);
        GStyle.getAxisAttributesZ().setTitleFontSize(46);
        GStyle.getAxisAttributesZ().setLabelFontSize(42);

        // Not available in your GROOT version
        // GStyle.getLegendAttributes().setFontSize(42);

        for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
            pad.setTitle("");

        if (name.equals("Occupancy")) {
            this.getCanvas(name).draw(this.getHistos().get(name));
            this.setupColorAxis(name, "Occupancy [%]");
        }

        if (name.equals("NPhe per PMT")) {
            for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads())
                pad.getAxisY().setLog(true);
            this.setLegend(name, 1100, 150);
        }

        if (name.equals("Mean NPhe by Ring")) {
            for (EmbeddedPad pad : this.getCanvas(name).getCanvasPads()) {
                pad.getAxisX().setRange(0.5, 4.5);
                pad.getAxisY().setRange(0, 50);
            }
        }

        if (name.equals("Rate vs Ring (Threshold)")) {
            EmbeddedCanvas canvas = this.getCanvas(name);
            canvas.setAxisTitleSize(52);
            canvas.setAxisLabelSize(48);
            PadMargins margins = new PadMargins()
                    .setLeftMargin(120)
                    .setBottomMargin(100)
                    .setTopMargin(30)
                    .setRightMargin(30);
            margins.setFixed(true);
            for (EmbeddedPad pad : canvas.getCanvasPads()) {
                pad.getAxisX().setRange(0.5, 5.8);
                pad.getAxisY().setLog(true);
                pad.setLegend(false); // replaced by manual LatexText below
                pad.setMargins(margins);
            }
            // GStyle.getLegendAttributes().setFontSize() is not available in this
            // GROOT version, so draw labels manually at 3× default size (~42 pt).
            canvas.cd(0);
            LatexText lbl0 = new LatexText("NPhe > 0", 1100, 140);
            lbl0.setFontSize(42);
            lbl0.setColor(4); // blue — matches g0 marker color
            canvas.draw(lbl0);
            LatexText lbl3 = new LatexText("NPhe > 3", 1100, 200);
            lbl3.setFontSize(42);
            lbl3.setColor(2); // red — matches g10 marker color
            canvas.draw(lbl3);
        }
    }

    // -----------------------------------------------------------------------
    // Fill methods
    // -----------------------------------------------------------------------

    public void fillOccupancy(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits)
            // sector = sector (1-6), ring = component (1-4)
            fillH2(dg, "hi_occ_htcc", h.getSector(), h.getComponent());
    }

    public void fillRateBySector(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            String pname = particleName(h);
            fillH1(dg, "hi_rate_sector_all", h.getSector());
            if (pname != null)
                fillH1(dg, "hi_rate_sector_" + pname, h.getSector());
            else
                fillH1(dg, "hi_rate_sector_other", h.getSector());
        }
    }

    public void fillRateByRing(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            // ring = component (1-4); layer is the half (1-2)
            String pname = particleName(h);
            fillH1(dg, "hi_rate_ring_all", h.getComponent());
            if (pname != null)
                fillH1(dg, "hi_rate_ring_" + pname, h.getComponent());
            else
                fillH1(dg, "hi_rate_ring_other", h.getComponent());
        }
    }

    public void fillNPhe(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            double nphe = h.getADC() / ADC_TO_NPHE;
            if (nphe > 0)
                fillH1(dg, "hi_nphe", nphe);
        }
    }

    public void fillOrigin(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            True t = h.getTrue();
            if (t == null)
                continue;

            double vx = t.getVertex().x();
            double vy = t.getVertex().y();
            double vz = t.getVertex().z();
            double r = Math.sqrt(vx * vx + vy * vy);

            fillH2(dg, "hi_bg_rz", vz, r);
            fillH2(dg, "hi_bg_xy", vx, vy);

            fillH1(dg, "hi_bg_vz_all", vz);
            fillH1(dg, "hi_bg_energy_all", t.getKinEnergy());

            String pname = particleName(h);
            if (pname != null) {
                fillH1(dg, "hi_bg_vz_" + pname, vz);
                fillH1(dg, "hi_bg_energy_" + pname, t.getKinEnergy());
            } else {
                fillH1(dg, "hi_bg_vz_other", vz);
                fillH1(dg, "hi_bg_energy_other", t.getKinEnergy());
            }
        }
    }

    public void fillSectorBG(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            True t = h.getTrue();
            if (t == null)
                continue;
            int s = h.getSector();

            fillH1(dg, "hi_bg_s" + s + "_all", t.getVertex().z());
            String pname = particleName(h);
            if (pname != null)
                fillH1(dg, "hi_bg_s" + s + "_" + pname, t.getVertex().z());
            else
                fillH1(dg, "hi_bg_s" + s + "_other", t.getVertex().z());
        }
    }

    public void fillNphePerPMT(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            double nphe = h.getADC() / ADC_TO_NPHE;
            if (nphe <= 0)
                continue;
            fillH1(dg, "hi_nphe_pmt_all", nphe);
            // ring = component (1-4)
            int ring = h.getComponent();
            if (ring >= 1 && ring <= NRINGS) {
                fillH1(dg, "hi_nphe_pmt_ring" + ring, nphe);
                totalNpheByRing[ring] += (long) nphe;
                countHitsByRing[ring]++;
            }
        }
    }

    public void fillRateVsRingThreshold(DataGroup dg, List<Hit> hits) {
        for (Hit h : hits) {
            double nphe = h.getADC() / ADC_TO_NPHE;

            int ring = h.getComponent(); // 1..4
            if (ring < 1 || ring > NRINGS)
                continue;

            double x = ring + 0.5; // fills bins 1-2, 2-3, 3-4, 4-5

            if (nphe > 0) {
                fillH1(dg, "hi_rvr_nphe0", x);
                rawNphe0[ring - 1]++;
            }
            if (nphe > 3) {
                fillH1(dg, "hi_rvr_nphe3", x);
                rawNphe10[ring - 1]++;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /** Fills a 1D histogram by name, skipping silently if it wasn't booked. */
    private void fillH1(DataGroup dg, String name, double val) {
        H1F h = dg.getH1F(name);
        if (h != null)
            h.fill(val);
    }

    /** Fills a 2D histogram by name, skipping silently if it wasn't booked. */
    private void fillH2(DataGroup dg, String name, double x, double y) {
        H2F h = dg.getH2F(name);
        if (h != null)
            h.fill(x, y);
    }

    /**
     * Returns the named particle category for a hit.
     * For optical-photon hits (pid == 0) the mother particle pid (mpid) is used
     * so that the Cherenkov-producing track species is reported.
     */
    private String particleName(Hit h) {
        True t = h.getTrue();
        if (t == null)
            return null;
        int pid = Math.abs(t.getPid());
        if (pid == 0)
            pid = Math.abs(t.getMPID());
        return this.pidToName(pid);
    }
}
