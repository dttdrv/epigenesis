import AppKit
import SwiftUI

// MARK: - 1. Develop Pipeline Workbench

struct DevelopView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        let req: JSONDict = [
            "mode": "develop",
            "preset": engine.devPreset,
            "accessibility": engine.devAccessibility,
            "samples": engine.devSamples,
            "include_coding_truncation": engine.devIncludeCoding,
            "run_cx3d": engine.devRunCx3d,
        ]
        engine.runScore(tab: .develop, request: req, summaryLabel: "\(engine.devPreset) (acc=\(String(format: "%.2f", engine.devAccessibility)))")
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                // Parameter Inspector Ribbon
                LabSectionCard(title: engine.t("tab.develop"), subtitle: engine.t("dev.sub"), trailing: AnyView(
                    Button(action: triggerRun) {
                        Label(engine.t("act.run"), systemImage: "play.fill")
                    }
                    .buttonStyle(.borderedProminent)
                    .controlSize(.small)
                    .keyboardShortcut(.return, modifiers: .command)
                    .disabled(engine.isBusy || !engine.isReady)
                )) {
                    VStack(spacing: 10) {
                        HStack(spacing: 16) {
                            VStack(alignment: .leading, spacing: 3) {
                                Text(engine.t("lbl.preset")).font(.system(size: 11, weight: .medium)).foregroundStyle(.secondary)
                                Picker("", selection: $engine.devPreset) {
                                    Text("neural-tube-multi-locus (Nkx2.2 / Olig2 / Pax6 / Irx3)").tag("neural-tube-multi-locus")
                                    Text("tri-lineage-chromatin (Excitatory / Inhibitory / Progenitor)").tag("tri-lineage-chromatin")
                                }
                                .labelsHidden()
                                .frame(minWidth: 280)
                            }

                            VStack(alignment: .leading, spacing: 3) {
                                HStack {
                                    Text(engine.t("lbl.accessibility")).font(.system(size: 11, weight: .medium)).foregroundStyle(.secondary)
                                    Spacer()
                                    Text(String(format: "%.2f", engine.devAccessibility)).font(.system(size: 11, weight: .semibold, design: .monospaced))
                                }
                                Slider(value: $engine.devAccessibility, in: 0...1, step: 0.05)
                                    .frame(width: 170)
                            }

                            VStack(alignment: .leading, spacing: 3) {
                                Text(engine.t("lbl.samples")).font(.system(size: 11, weight: .medium)).foregroundStyle(.secondary)
                                Stepper(value: $engine.devSamples, in: 4...64, step: 4) {
                                    Text("\(engine.devSamples)").font(.system(size: 12, design: .monospaced))
                                }
                                .frame(width: 100)
                            }

                            Spacer()
                        }

                        Divider()

                        HStack(spacing: 20) {
                            Toggle(engine.t("lbl.coding"), isOn: $engine.devIncludeCoding)
                                .toggleStyle(.checkbox)
                                .font(.system(size: 12))
                            Toggle(engine.t("lbl.cx3d"), isOn: $engine.devRunCx3d)
                                .toggleStyle(.checkbox)
                                .font(.system(size: 12))
                            Spacer()
                        }
                    }
                }

                if let receipt = engine.receipts[.develop], let res = receipt.dict("result") {
                    let nom = res.dict("nominal") ?? [:]
                    let ctrl = res.dict("controls") ?? [:]
                    let ens = res.dict("uncertainty_ensemble") ?? [:]
                    let prov = res.dict("provenance") ?? [:]
                    let rust = prov.dict("rust_cross_check") ?? [:]

                    // KPI Strip
                    HStack(spacing: 10) {
                        MetricTile(
                            label: "Cell Fate Transition",
                            value: "\(nom.str("wt_fate")) → \(nom.str("variant_fate"))",
                            detail: nom.flag("fate_flipped") ? "Bifurcation crossed" : "Fate preserved",
                            tone: nom.flag("fate_flipped") ? .orange : .green
                        )
                        MetricTile(
                            label: "Ensemble Flip Prob.",
                            value: String(format: "%.1f%%", ens.dbl("fate_flip_probability") * 100),
                            detail: "\(engine.devSamples) parameter draws"
                        )
                        MetricTile(
                            label: "Causal Controls",
                            value: ctrl.flag("reversion_exact") && ctrl.flag("closed_chromatin_cis_exact") ? "All Exact (Δ=0)" : "Check",
                            detail: "Reversion & closed chromatin",
                            tone: .green
                        )
                        MetricTile(
                            label: "Cross-Language Replay",
                            value: (prov.flag("development_valid") && rust.flag("development_valid")) ? "Python + Rust OK" : (prov.flag("development_valid") ? "Python OK" : "Failed"),
                            detail: "GA4GH VRS 2.x + brainc",
                            tone: prov.flag("development_valid") ? .green : .red
                        )
                    }

                    // GRN Steady-State & 90% CI Chart + Table
                    let genes = res.strArr("genes")
                    let wtEnd = nom.dict("wt_end_state") ?? [:]
                    let varEnd = nom.dict("variant_end_state") ?? [:]
                    let deltaSum = ens.dict("gene_delta_summary") ?? [:]
                    let chartRows: [GRNUncertaintyChartView.Row] = genes.map { g in
                        let ds = deltaSum.dict(g) ?? [:]
                        return GRNUncertaintyChartView.Row(
                            id: g,
                            wt: wtEnd.dbl(g),
                            variant: varEnd.dbl(g),
                            deltaMean: ds.dbl("mean"),
                            p05: ds.dbl("p05"),
                            p95: ds.dbl("p95")
                        )
                    }

                    LabSectionCard(title: "Gene Regulatory Network Steady State & 90% Uncertainty Ensemble", subtitle: "Top bar: wild-type (blue) · Bottom bar: edited state with p05–p95 whisker") {
                        GRNUncertaintyChartView(rows: chartRows)
                    }

                    // 3D CX3D Morphogenesis & Provenance Inspector
                    HStack(alignment: .top, spacing: 12) {
                        if let morph = res.dict("morphogenesis"), morph.flag("executed"), let branches = morph.dict("branches") {
                            LabSectionCard(title: "3D CX3D Morphogenesis & LIF Spiking") {
                                VStack(alignment: .leading, spacing: 6) {
                                    ForEach(["reference", "variant", "allele_rescue"], id: \.self) { bKey in
                                        if let b = branches.dict(bKey), let cx = b.dict("construction") {
                                            HStack {
                                                Text(bKey).font(.system(size: 11.5, weight: .semibold, design: .monospaced)).frame(width: 96, alignment: .leading)
                                                Text("Cells: \(cx.int("cells"))").font(.system(size: 11, design: .monospaced))
                                                Text(String(format: "Neurite: %.2f", cx.dbl("length"))).font(.system(size: 11, design: .monospaced))
                                                Text("Synapses: \(cx.int("synapses"))").font(.system(size: 11, design: .monospaced))
                                                Spacer()
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        LabSectionCard(title: "Cryptographic Provenance & Controls") {
                            VStack(alignment: .leading, spacing: 6) {
                                HStack(spacing: 8) {
                                    BadgePill(text: "VRS 2.x Valid", ok: prov.flag("variants_valid"))
                                    BadgePill(text: "Bundle Valid", ok: prov.flag("development_valid"))
                                    if rust.flag("available") {
                                        BadgePill(text: "Rust Verified", ok: rust.flag("development_valid"))
                                    }
                                }
                                Text("Bundle SHA-256: \(prov.str("development_bundle_sha256"))")
                                    .font(.system(size: 10.5, design: .monospaced))
                                    .foregroundStyle(.secondary)
                                    .textSelection(.enabled)
                                Text("Receipt SHA-256: \(receipt.str("receipt_sha256"))")
                                    .font(.system(size: 10.5, design: .monospaced))
                                    .foregroundStyle(.secondary)
                                    .textSelection(.enabled)
                            }
                        }
                    }
                } else {
                    Text(engine.t("lbl.noData"))
                        .font(.system(size: 12.5))
                        .foregroundStyle(.secondary)
                        .padding(24)
                        .frame(maxWidth: .infinity)
                }
            }
            .padding(14)
        }
    }
}

// MARK: - 2. Inverse Compiler & 532-TF Discovery Workbench

struct InverseCompilerView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        if engine.invGoal == "discover" {
            engine.runScore(tab: .inverse, request: ["mode": "discover", "bootstraps": 200], summaryLabel: "532-TF GSE188264 Discovery")
        } else if engine.invGoal == "rescue" {
            let req: JSONDict = [
                "mode": "synthesize",
                "goal": "rescue",
                "locus": engine.invLocus,
                "breaking_edit": engine.invBreakingEdit,
                "signal": engine.invSignal,
                "run_cx3d": engine.invRunCx3d,
            ]
            engine.runScore(tab: .inverse, request: req, summaryLabel: "Rescue \(engine.invLocus) (\(engine.invBreakingEdit))")
        } else {
            let req: JSONDict = [
                "mode": "synthesize",
                "goal": "reprogram",
                "locus": engine.invLocus,
                "target_fate": engine.invTargetFate,
                "signal": engine.invSignal,
                "run_cx3d": engine.invRunCx3d,
            ]
            engine.runScore(tab: .inverse, request: req, summaryLabel: "Reprogram \(engine.invLocus) → \(engine.invTargetFate)")
        }
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                LabSectionCard(title: engine.t("tab.inverse"), subtitle: engine.t("inv.sub"), trailing: AnyView(
                    Button(action: triggerRun) {
                        Label(engine.t("act.run"), systemImage: "play.fill")
                    }
                    .buttonStyle(.borderedProminent)
                    .controlSize(.small)
                    .keyboardShortcut(.return, modifiers: .command)
                    .disabled(engine.isBusy || !engine.isReady)
                )) {
                    VStack(alignment: .leading, spacing: 10) {
                        Picker("Goal", selection: $engine.invGoal) {
                            Text("Epistatic 3D Rescue (Primary Lesion Intact)").tag("rescue")
                            Text("Non-Coding Fate Reprogramming").tag("reprogram")
                            Text("532-TF Human Enhancer Discovery (GSE188264)").tag("discover")
                        }
                        .pickerStyle(.segmented)

                        if engine.invGoal != "discover" {
                            HStack(spacing: 14) {
                                VStack(alignment: .leading, spacing: 3) {
                                    Text(engine.t("lbl.locus")).font(.system(size: 11)).foregroundStyle(.secondary)
                                    Picker("", selection: $engine.invLocus) {
                                        Text("chr18:56475717-56475888").tag("chr18:56475717-56475888")
                                        Text("chr1:175099416-175099587").tag("chr1:175099416-175099587")
                                        Text("chr15:90936956-90937127").tag("chr15:90936956-90937127")
                                    }
                                    .labelsHidden()
                                    .frame(width: 220)
                                }

                                if engine.invGoal == "rescue" {
                                    VStack(alignment: .leading, spacing: 3) {
                                        Text(engine.t("lbl.breakingEdit")).font(.system(size: 11)).foregroundStyle(.secondary)
                                        TextField("165:G>T", text: $engine.invBreakingEdit)
                                            .font(.system(size: 12, design: .monospaced))
                                            .textFieldStyle(.roundedBorder)
                                            .frame(width: 110)
                                    }
                                } else {
                                    VStack(alignment: .leading, spacing: 3) {
                                        Text(engine.t("lbl.targetFate")).font(.system(size: 11)).foregroundStyle(.secondary)
                                        Picker("", selection: $engine.invTargetFate) {
                                            Text("Olig2 (pMN)").tag("Olig2")
                                            Text("Nkx2.2 (p3)").tag("Nkx2.2")
                                            Text("Pax6").tag("Pax6")
                                        }
                                        .labelsHidden()
                                        .frame(width: 130)
                                    }
                                }

                                VStack(alignment: .leading, spacing: 3) {
                                    Text("\(engine.t("lbl.signal")): \(String(format: "%.2f", engine.invSignal))").font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary)
                                    Slider(value: $engine.invSignal, in: 0.05...0.80, step: 0.05).frame(width: 140)
                                }

                                Toggle(engine.t("lbl.cx3d"), isOn: $engine.invRunCx3d)
                                    .toggleStyle(.checkbox)
                                    .font(.system(size: 12))
                                    .padding(.top, 14)

                                Spacer()
                            }
                        }
                    }
                }

                if let receipt = engine.receipts[.inverse], let res = receipt.dict("result") {
                    if let heldout = res.dict("heldout_test") {
                        // Discovery result view
                        HStack(spacing: 10) {
                            MetricTile(label: "Held-Out Pearson r", value: String(format: "%.4f", heldout.dbl("pearson_r")), detail: "1,895 unseen human edits", tone: .accentColor)
                            MetricTile(label: "Top-10% Direction", value: String(format: "%.1f%%", heldout.dbl("top_decile_direction_accuracy") * 100), detail: "vs 50.0% chance", tone: .green)
                            MetricTile(label: "MSE Gain vs Zero", value: String(format: "%.2f%%", heldout.dbl("mse_reduction_vs_zero_pct")), detail: "Pure 532-TF Boltzmann")
                            MetricTile(label: "Sparse Top-25 TF r", value: String(format: "%.4f", heldout.dbl("sparse_top25_pearson_r")), detail: "Interpretable TF subset")
                        }
                        if let tfs = res.dict("top_tfs") {
                            HStack(alignment: .top, spacing: 12) {
                                LabSectionCard(title: "Top Activating TFs (JASPAR 2024)") {
                                    VStack(alignment: .leading, spacing: 4) {
                                        ForEach(Array(tfs.dictArr("activators").prefix(8).enumerated()), id: \.offset) { _, row in
                                            HStack {
                                                Text(row.str("motif_id")).font(.system(size: 11.5, weight: .semibold, design: .monospaced))
                                                Spacer()
                                                Text(String(format: "%+.4f", row.dbl("weight"))).font(.system(size: 11.5, design: .monospaced)).foregroundStyle(.orange)
                                            }
                                        }
                                    }
                                }
                                LabSectionCard(title: "Top Repressing TFs (JASPAR 2024)") {
                                    VStack(alignment: .leading, spacing: 4) {
                                        ForEach(Array(tfs.dictArr("repressors").prefix(8).enumerated()), id: \.offset) { _, row in
                                            HStack {
                                                Text(row.str("motif_id")).font(.system(size: 11.5, weight: .semibold, design: .monospaced))
                                                Spacer()
                                                Text(String(format: "%+.4f", row.dbl("weight"))).font(.system(size: 11.5, design: .monospaced)).foregroundStyle(.blue)
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        // Synthesize Rescue / Reprogram view
                        let syn = res.dict("synthesized_edit") ?? res.dict("synthesized_rescue") ?? [:]
                        let prov = res.dict("provenance") ?? [:]
                        let rust = prov.dict("rust_cross_check") ?? [:]
                        HStack(spacing: 10) {
                            MetricTile(label: "Compiled DNA Edit", value: syn.str("edit_id", syn.str("rescue_edit", "—")), detail: syn.str("vrs_allele_id", "").prefix(22).description, tone: .accentColor)
                            MetricTile(label: "Affinity Ratio q", value: String(format: "%.5f", syn.dbl("affinity_ratio", syn.dbl("double_mutant_q"))), detail: "Strand-symmetric Boltzmann")
                            MetricTile(label: "Outcome", value: "\(syn.str("nominal_wt_fate", res.str("wt_fate"))) → \(syn.str("nominal_variant_fate", res.str("target_fate", "Restored")))", detail: "Primary constraint satisfied", tone: .green)
                            MetricTile(label: "Cross-Language Proof", value: (prov.flag("development_valid") && rust.flag("development_valid")) ? "Python + Rust OK" : "Python OK", detail: "RFC 8785 SHA-256", tone: .green)
                        }

                        let pool = res.dictArr("top_pool_evaluated")
                        if !pool.isEmpty {
                            LabSectionCard(title: "Ranked Candidate Designs Evaluated", subtitle: "Sorted by biophysical target fit & ensemble stability") {
                                VStack(spacing: 4) {
                                    HStack {
                                        Text("Edit").frame(width: 90, alignment: .leading)
                                        Text("Category").frame(width: 140, alignment: .leading)
                                        Text("Affinity q").frame(width: 100, alignment: .trailing)
                                        Text("Ensemble Prob.").frame(width: 110, alignment: .trailing)
                                        Spacer()
                                    }
                                    .font(.system(size: 10.5, weight: .semibold))
                                    .foregroundStyle(.secondary)
                                    Divider()
                                    ForEach(Array(pool.prefix(10).enumerated()), id: \.offset) { _, c in
                                        HStack {
                                            Text(c.str("edit_id")).font(.system(size: 11.5, weight: .semibold, design: .monospaced)).frame(width: 90, alignment: .leading)
                                            Text(c.str("category", "candidate")).font(.system(size: 11)).frame(width: 140, alignment: .leading)
                                            Text(String(format: "%.5f", c.dbl("affinity_ratio"))).font(.system(size: 11, design: .monospaced)).frame(width: 100, alignment: .trailing)
                                            Text(String(format: "%.1f%%", c.dbl("target_fate_ensemble_prob", 1.0) * 100)).font(.system(size: 11, design: .monospaced)).frame(width: 110, alignment: .trailing)
                                            Spacer()
                                        }
                                        .padding(.vertical, 2)
                                    }
                                }
                            }
                        }
                    }
                } else {
                    Text(engine.t("lbl.noData")).font(.system(size: 12.5)).foregroundStyle(.secondary).padding(24).frame(maxWidth: .infinity)
                }
            }
            .padding(14)
        }
    }
}

// MARK: - 3. Regulatory Site Simulation (9-bp GLI)

struct SimulateSiteView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        let req: JSONDict = [
            "mode": "simulate",
            "site": engine.simSite,
            "protein": engine.simProtein,
            "phases": [[engine.simDuration, engine.simSignal]],
        ]
        engine.runScore(tab: .simulate, request: req, summaryLabel: "\(engine.simSite) (\(engine.simProtein), S=\(String(format: "%.2f", engine.simSignal)))")
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                LabSectionCard(title: engine.t("tab.simulate"), subtitle: engine.t("sim.sub"), trailing: AnyView(
                    Button(action: triggerRun) {
                        Label(engine.t("act.run"), systemImage: "play.fill")
                    }
                    .buttonStyle(.borderedProminent)
                    .controlSize(.small)
                    .keyboardShortcut(.return, modifiers: .command)
                    .disabled(engine.isBusy || !engine.isReady)
                )) {
                    VStack(alignment: .leading, spacing: 10) {
                        HStack(spacing: 16) {
                            VStack(alignment: .leading, spacing: 3) {
                                Text("9-bp GLI Site (Ref: GACCACCCA)").font(.system(size: 11)).foregroundStyle(.secondary)
                                TextField("GAACACCCA", text: $engine.simSite)
                                    .font(.system(size: 13, weight: .semibold, design: .monospaced))
                                    .textFieldStyle(.roundedBorder)
                                    .frame(width: 150)
                            }
                            VStack(alignment: .leading, spacing: 3) {
                                Text(engine.t("lbl.protein")).font(.system(size: 11)).foregroundStyle(.secondary)
                                Picker("", selection: $engine.simProtein) {
                                    Text("GLI2 (Primary)").tag("GLI2")
                                    Text("GLI3 (Sensitivity)").tag("GLI3")
                                }
                                .labelsHidden()
                                .frame(width: 150)
                            }
                            VStack(alignment: .leading, spacing: 3) {
                                Text("Duration: \(Int(engine.simDuration)) u").font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary)
                                Slider(value: $engine.simDuration, in: 5...60, step: 5).frame(width: 120)
                            }
                            VStack(alignment: .leading, spacing: 3) {
                                Text("Signal S: \(String(format: "%.2f", engine.simSignal))").font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary)
                                Slider(value: $engine.simSignal, in: 0.05...0.95, step: 0.05).frame(width: 140)
                            }
                            Spacer()
                        }
                        DNASequenceTrackView(sequence: engine.simSite, highlightPos1Based: 3, altBase: String(engine.simSite.dropFirst(2).prefix(1)))
                    }
                }

                if let receipt = engine.receipts[.simulate], let res = receipt.dict("result") {
                    let genes = res.strArr("genes")
                    let times = res.dblArr("time")
                    let runs = res.dict("runs") ?? [:]
                    let refRun = runs.dict("reference") ?? [:]
                    let varRun = runs.dict("variant") ?? [:]
                    let repRun = runs.dict("reporter") ?? [:]
                    let refStates = refRun.arr("states").compactMap { ($0 as? [Any])?.compactMap { ($0 as? Double) ?? ($0 as? Int).map(Double.init) } }
                    let varStates = varRun.arr("states").compactMap { ($0 as? [Any])?.compactMap { ($0 as? Double) ?? ($0 as? Int).map(Double.init) } }
                    let sInfo = res.dict("site") ?? [:]
                    let ctrl = res.dict("controls") ?? [:]

                    HStack(spacing: 10) {
                        MetricTile(label: "Measured Affinity Ratio q", value: String(format: "%.5f", sInfo.dbl("ratio")), detail: "Hallikas \(engine.simProtein) Table S1")
                        MetricTile(label: "Endogenous Nkx2.2", value: String(format: "%.4f → %.4f", refStates.last?[2] ?? 0, varStates.last?[2] ?? 0), detail: "Coupled 4-gene feedback", tone: .orange)
                        MetricTile(label: "Passive Reporter Occ.", value: String(format: "%.4f → %.4f", refRun.dblArr("occupancy").last ?? 0, repRun.dblArr("occupancy").last ?? 0), detail: "Unchanged reference network")
                        MetricTile(label: "Reversion Control", value: ctrl.dbl("reversion_max_difference") == 0 ? "Exact (0.0)" : "Check", detail: "Identical brainc module", tone: .green)
                    }

                    LabSectionCard(title: "4-Gene Neural Tube ODE Trajectories (Hover to inspect concentrations)") {
                        ODETrajectoryPlotView(times: times, genes: genes, refStates: refStates, varStates: varStates, hoverIdx: $engine.simHoverIdx)
                    }
                } else {
                    Text(engine.t("lbl.noData")).font(.system(size: 12.5)).foregroundStyle(.secondary).padding(24).frame(maxWidth: .infinity)
                }
            }
            .padding(14)
        }
    }
}

// MARK: - 4. Enhancer Circuit Scanner

struct EnhancerCircuitView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        let req: JSONDict = [
            "mode": "circuit",
            "sequence": engine.circuitSequence,
            "circuit": engine.cirName,
            "protein": engine.cirProtein,
            "phases": [[20.0, engine.cirSignal]],
        ]
        engine.runScore(tab: .circuit, request: req, summaryLabel: "Circuit \(engine.cirName) (\(engine.cirProtein), \(engine.circuitSequence.count) bp)")
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                LabSectionCard(title: engine.t("tab.circuit"), subtitle: engine.t("cir.sub"), trailing: AnyView(
                    Button(action: triggerRun) {
                        Label(engine.t("act.run"), systemImage: "play.fill")
                    }
                    .buttonStyle(.borderedProminent)
                    .controlSize(.small)
                    .keyboardShortcut(.return, modifiers: .command)
                    .disabled(engine.isBusy || !engine.isReady)
                )) {
                    VStack(alignment: .leading, spacing: 8) {
                        HStack(spacing: 14) {
                            Picker(engine.t("lbl.circuit"), selection: $engine.cirName) {
                                Text("neural-tube (4 genes)").tag("neural-tube")
                                Text("mutual-toggle (2 genes)").tag("mutual-toggle")
                                Text("tri-lineage-switch (3 genes)").tag("tri-lineage-switch")
                            }
                            .frame(width: 240)

                            Picker(engine.t("lbl.protein"), selection: $engine.cirProtein) {
                                Text("GLI2").tag("GLI2")
                                Text("GLI3").tag("GLI3")
                                Text("MA0069.1 (PAX6)").tag("MA0069.1")
                                Text("MA0463.3 (BCL6)").tag("MA0463.3")
                            }
                            .frame(width: 190)

                            Text("S = \(String(format: "%.2f", engine.cirSignal))").font(.system(size: 11, design: .monospaced))
                            Slider(value: $engine.cirSignal, in: 0.05...0.80, step: 0.05).frame(width: 120)
                            Spacer()
                        }

                        TextEditor(text: $engine.circuitSequence)
                            .font(.system(size: 11.5, design: .monospaced))
                            .frame(height: 58)
                            .overlay(RoundedRectangle(cornerRadius: 4).stroke(Color(nsColor: .separatorColor)))
                    }
                }

                if let receipt = engine.receipts[.circuit], let res = receipt.dict("result") {
                    let counts = res.dict("category_counts") ?? [:]
                    let vars = res.dictArr("variants")
                    let cryptic = vars.filter { $0.str("category") == "cryptic_bifurcation" }

                    HStack(spacing: 10) {
                        MetricTile(label: "Variants Scanned", value: "\(vars.count)", detail: "\(res.int("length")) bp enhancer")
                        MetricTile(label: "Cryptic Bifurcations", value: "\(counts.int("cryptic_bifurcation"))", detail: "|Δreporter|<0.02 yet fate flips", tone: counts.int("cryptic_bifurcation") > 0 ? .orange : .primary)
                        MetricTile(label: "Concordant Switches", value: "\(counts.int("concordant_switch"))", detail: "Reporter & ODE agree")
                        MetricTile(label: "Neutral / Buffered", value: "\(counts.int("neutral") + counts.int("canalized_buffered"))", detail: "Within stability basin")
                    }

                    LabSectionCard(title: "Top Circuit-Perturbing Variants", subtitle: "Cryptic bifurcations shown first") {
                        let shown = Array((cryptic + vars.filter { $0.str("category") != "cryptic_bifurcation" }).prefix(16))
                        VStack(spacing: 3) {
                            HStack {
                                Text("Edit").frame(width: 80, alignment: .leading)
                                Text("Classification").frame(width: 150, alignment: .leading)
                                Text("Affinity q").frame(width: 90, alignment: .trailing)
                                Text("Δ Reporter").frame(width: 90, alignment: .trailing)
                                Text("Δ Endogenous").frame(width: 100, alignment: .trailing)
                                Spacer()
                            }
                            .font(.system(size: 10.5, weight: .semibold)).foregroundStyle(.secondary)
                            Divider()
                            ForEach(Array(shown.enumerated()), id: \.offset) { _, v in
                                HStack {
                                    Text(v.str("id")).font(.system(size: 11, weight: .semibold, design: .monospaced)).frame(width: 80, alignment: .leading)
                                    Text(v.str("category")).font(.system(size: 11)).foregroundStyle(v.str("category") == "cryptic_bifurcation" ? .orange : .secondary).frame(width: 150, alignment: .leading)
                                    Text(String(format: "%.5f", v.dbl("affinity_ratio"))).font(.system(size: 11, design: .monospaced)).frame(width: 90, alignment: .trailing)
                                    Text(String(format: "%+.4f", v.dbl("delta_reporter"))).font(.system(size: 11, design: .monospaced)).frame(width: 90, alignment: .trailing)
                                    Text(String(format: "%+.4f", v.dbl("delta_endogenous"))).font(.system(size: 11, design: .monospaced)).frame(width: 100, alignment: .trailing)
                                    Spacer()
                                }
                                .padding(.vertical, 2)
                            }
                        }
                    }
                } else {
                    Text(engine.t("lbl.noData")).font(.system(size: 12.5)).foregroundStyle(.secondary).padding(24).frame(maxWidth: .infinity)
                }
            }
            .padding(14)
        }
    }
}

// MARK: - 5. Positional Saturation Scan Workbench

struct PositionalScanView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        let req: JSONDict = [
            "mode": "scan",
            "model": engine.scanPredictor,
            "sequence": engine.scanSequence,
        ]
        engine.scanSelectedIdx = nil
        engine.runScore(tab: .scan, request: req, summaryLabel: "Scan (\(engine.scanPredictor), \(engine.scanSequence.count) bp)")
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                LabSectionCard(title: engine.t("tab.scan"), subtitle: engine.t("scan.sub"), trailing: AnyView(
                    HStack(spacing: 8) {
                        Picker("", selection: $engine.scanPredictor) {
                            Text("Biophysical 532-TF+ v2 (r=0.51, 99% Top-100)").tag("biophys")
                            Text("Kreimer 0–72h Motif-Temporal").tag("kreimer")
                            Text("Siraj SK-N-SH Huber-Linear").tag("siraj")
                        }
                        .labelsHidden()
                        .frame(width: 290)

                        Button(action: triggerRun) {
                            Label(engine.t("act.run"), systemImage: "play.fill")
                        }
                        .buttonStyle(.borderedProminent)
                        .controlSize(.small)
                        .keyboardShortcut(.return, modifiers: .command)
                        .disabled(engine.isBusy || !engine.isReady)
                    }
                )) {
                    TextEditor(text: $engine.scanSequence)
                        .font(.system(size: 11.5, design: .monospaced))
                        .frame(height: 62)
                        .overlay(RoundedRectangle(cornerRadius: 4).stroke(Color(nsColor: .separatorColor)))
                }

                if let receipt = engine.receipts[.scan], let res = receipt.dict("result") {
                    let preds = res.dictArr("predictions")
                    let cleanSeq = receipt.dict("request")?.str("sequence") ?? engine.scanSequence
                    let sorted = preds.sorted { abs($0.dbl("effect")) > abs($1.dbl("effect")) }
                    let topAct = preds.max { $0.dbl("effect") < $1.dbl("effect") }
                    let topRep = preds.min { $0.dbl("effect") < $1.dbl("effect") }

                    HStack(spacing: 10) {
                        MetricTile(label: "Substitutions Scored", value: "\(preds.count)", detail: "\(cleanSeq.count) bp · \(engine.scanPredictor)")
                        MetricTile(label: "Strongest Activation", value: String(format: "%+.3f", topAct?.dbl("effect") ?? 0), detail: topAct?.str("id") ?? "—", tone: .orange)
                        MetricTile(label: "Strongest Repression", value: String(format: "%+.3f", topRep?.dbl("effect") ?? 0), detail: topRep?.str("id") ?? "—", tone: .blue)
                        MetricTile(label: "Receipt Seal", value: String(receipt.str("receipt_sha256").prefix(10)), detail: "RFC 8785 SHA-256", tone: .green)
                    }

                    LabSectionCard(title: "Positional Saturation Mutagenesis Heatmap") {
                        SaturationHeatmapView(sequence: cleanSeq, predictions: preds, selectedIndex: $engine.scanSelectedIdx)
                    }

                    let activePred: JSONDict? = engine.scanSelectedIdx.flatMap { $0 < preds.count ? preds[$0] : nil } ?? sorted.first
                    if let sel = activePred {
                        let idParts = sel.str("id").split(separator: ":")
                        let selPos = idParts.count == 2 ? (Int(idParts[0]) ?? 0) : sel.int("position")
                        let chParts = idParts.count == 2 ? idParts[1].split(separator: ">") : []
                        let selAlt = chParts.count == 2 ? String(chParts[1]) : sel.str("alt")
                        HStack(alignment: .top, spacing: 12) {
                            LabSectionCard(title: "Selected Variant Inspector — \(sel.str("id"))") {
                                VStack(alignment: .leading, spacing: 8) {
                                    DNASequenceTrackView(sequence: cleanSeq, highlightPos1Based: selPos > 0 ? selPos : nil, altBase: selAlt.isEmpty ? nil : selAlt)
                                    HStack(spacing: 16) {
                                        Text(String(format: "Total Effect: %+.4f", sel.dbl("effect"))).font(.system(size: 12, weight: .bold, design: .monospaced))
                                        Text(String(format: "Sequence k-mer: %+.4f", sel.dbl("word_effect"))).font(.system(size: 11.5, design: .monospaced)).foregroundStyle(.secondary)
                                        Text(String(format: "TF Motif / Boltzmann: %+.4f", sel.dbl("motif_effect"))).font(.system(size: 11.5, design: .monospaced)).foregroundStyle(.secondary)
                                    }
                                    let topM = sel.dictArr("top_motifs")
                                    if !topM.isEmpty {
                                        Divider()
                                        Text("Top Contributing JASPAR 2024 TFs").font(.system(size: 11, weight: .semibold)).foregroundStyle(.secondary)
                                        ForEach(Array(topM.enumerated()), id: \.offset) { _, m in
                                            HStack {
                                                Text(m.str("id")).font(.system(size: 11, weight: .medium, design: .monospaced))
                                                Spacer()
                                                Text(String(format: "Δbind %+.3f   contrib %+.4f", m.dbl("delta"), m.dbl("effect")))
                                                    .font(.system(size: 11, design: .monospaced))
                                                    .foregroundStyle(m.dbl("effect") >= 0 ? .orange : .blue)
                                            }
                                        }
                                    }
                                }
                            }

                            let temp = sel.dblArr("temporal")
                            if temp.count == 7 {
                                LabSectionCard(title: "0–72 h Neural Induction Kinetics") {
                                    TemporalCurvePlotView(times: [0, 3, 6, 12, 24, 48, 72], values: temp)
                                }
                                .frame(width: 340)
                            }
                        }
                    }

                    // Ranked Data Table with Filter
                    LabSectionCard(title: "Ranked Variant Effects", trailing: AnyView(
                        TextField("Filter by position or base…", text: $engine.scanSearchText)
                            .textFieldStyle(.roundedBorder)
                            .controlSize(.small)
                            .frame(width: 180)
                    )) {
                        let filtered = sorted.filter { engine.scanSearchText.isEmpty || $0.str("id").localizedCaseInsensitiveContains(engine.scanSearchText) }
                        VStack(spacing: 2) {
                            HStack {
                                Text("Variant").frame(width: 90, alignment: .leading)
                                Text("Effect").frame(width: 90, alignment: .trailing)
                                Text("k-mer").frame(width: 90, alignment: .trailing)
                                Text("TF / Motif").frame(width: 90, alignment: .trailing)
                                Text("Top TF Driver").frame(maxWidth: .infinity, alignment: .leading).padding(.leading, 12)
                            }
                            .font(.system(size: 10.5, weight: .semibold)).foregroundStyle(.secondary)
                            Divider()
                            ForEach(Array(filtered.prefix(25).enumerated()), id: \.offset) { _, p in
                                HStack {
                                    Text(p.str("id")).font(.system(size: 11.5, weight: .semibold, design: .monospaced)).frame(width: 90, alignment: .leading)
                                    Text(String(format: "%+.4f", p.dbl("effect"))).font(.system(size: 11.5, weight: .semibold, design: .monospaced)).foregroundStyle(p.dbl("effect") >= 0 ? .orange : .blue).frame(width: 90, alignment: .trailing)
                                    Text(String(format: "%+.4f", p.dbl("word_effect"))).font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(width: 90, alignment: .trailing)
                                    Text(String(format: "%+.4f", p.dbl("motif_effect"))).font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(width: 90, alignment: .trailing)
                                    Text(p.dictArr("top_motifs").first?.str("id") ?? "—").font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(maxWidth: .infinity, alignment: .leading).padding(.leading, 12)
                                }
                                .padding(.vertical, 2)
                            }
                        }
                    }
                } else {
                    Text(engine.t("lbl.noData")).font(.system(size: 12.5)).foregroundStyle(.secondary).padding(24).frame(maxWidth: .infinity)
                }
            }
            .padding(14)
        }
    }
}

// MARK: - 6. Batch Variant Table Workbench

struct BatchEditsView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        guard let data = engine.editsText.data(using: .utf8),
              let rows = (try? JSONSerialization.jsonObject(with: data)) as? [JSONDict] else {
            engine.lastError = "Input must be a JSON array of {id, reference, sequence, family?} objects (or open a CSV/JSON via File → Open)."
            return
        }
        let req: JSONDict = [
            "mode": "pairs",
            "model": engine.editsPredictor,
            "rows": rows,
            "budget": min(engine.editsBudget, max(1, rows.count)),
        ]
        engine.runScore(tab: .edits, request: req, summaryLabel: "Batch \(rows.count) edits (\(engine.editsPredictor))")
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                LabSectionCard(title: engine.t("tab.edits"), subtitle: engine.t("edits.sub"), trailing: AnyView(
                    HStack(spacing: 8) {
                        Picker("", selection: $engine.editsPredictor) {
                            Text("Biophysical 532-TF+ v2").tag("biophys")
                            Text("Kreimer 0–72h").tag("kreimer")
                            Text("Siraj Huber-Linear").tag("siraj")
                        }
                        .labelsHidden()
                        .frame(width: 190)

                        Stepper("Panel: \(engine.editsBudget)", value: $engine.editsBudget, in: 1...50).font(.system(size: 11.5)).frame(width: 110)

                        Button(action: triggerRun) {
                            Label(engine.t("act.run"), systemImage: "play.fill")
                        }
                        .buttonStyle(.borderedProminent)
                        .controlSize(.small)
                        .keyboardShortcut(.return, modifiers: .command)
                        .disabled(engine.isBusy || !engine.isReady)
                    }
                )) {
                    TextEditor(text: $engine.editsText)
                        .font(.system(size: 11, design: .monospaced))
                        .frame(height: 110)
                        .overlay(RoundedRectangle(cornerRadius: 4).stroke(Color(nsColor: .separatorColor)))
                }

                if let receipt = engine.receipts[.edits], let res = receipt.dict("result") {
                    let sel = res.dictArr("selection")
                    let preds = res.dictArr("predictions")

                    LabSectionCard(title: "Selected Family-Diverse Panel (\(sel.count) edits from \(preds.count) scored)") {
                        VStack(spacing: 3) {
                            HStack {
                                Text("#").frame(width: 28, alignment: .leading)
                                Text("Variant ID").frame(width: 180, alignment: .leading)
                                Text("Family").frame(width: 110, alignment: .leading)
                                Text("Predicted Effect").frame(width: 110, alignment: .trailing)
                                Text("k-mer").frame(width: 85, alignment: .trailing)
                                Text("TF / Motif").frame(width: 85, alignment: .trailing)
                                Spacer()
                            }
                            .font(.system(size: 10.5, weight: .semibold)).foregroundStyle(.secondary)
                            Divider()
                            ForEach(Array(sel.enumerated()), id: \.offset) { idx, row in
                                HStack {
                                    Text("\(idx + 1)").font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(width: 28, alignment: .leading)
                                    Text(row.str("id")).font(.system(size: 11.5, weight: .semibold, design: .monospaced)).frame(width: 180, alignment: .leading)
                                    Text(row.str("family", "—")).font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(width: 110, alignment: .leading)
                                    Text(String(format: "%+.4f", row.dbl("effect"))).font(.system(size: 11.5, weight: .bold, design: .monospaced)).foregroundStyle(row.dbl("effect") >= 0 ? .orange : .blue).frame(width: 110, alignment: .trailing)
                                    Text(String(format: "%+.4f", row.dbl("word_effect"))).font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(width: 85, alignment: .trailing)
                                    Text(String(format: "%+.4f", row.dbl("motif_effect"))).font(.system(size: 11, design: .monospaced)).foregroundStyle(.secondary).frame(width: 85, alignment: .trailing)
                                    Spacer()
                                }
                                .padding(.vertical, 2)
                            }
                        }
                    }
                } else {
                    Text(engine.t("lbl.noData")).font(.system(size: 12.5)).foregroundStyle(.secondary).padding(24).frame(maxWidth: .infinity)
                }
            }
            .padding(14)
        }
    }
}

// MARK: - 7. Cryptographic Receipt Verifier

struct VerifyReceiptView: View {
    @ObservedObject var engine: EngineClient

    func triggerRun() {
        guard let data = engine.verifyText.data(using: .utf8),
              let obj = (try? JSONSerialization.jsonObject(with: data)) as? JSONDict else {
            engine.lastError = "Paste a valid JSON receipt or open one via File → Open (⌘O)."
            return
        }
        engine.runVerify(receiptObj: obj)
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 12) {
                LabSectionCard(title: engine.t("tab.verify"), subtitle: engine.t("ver.sub"), trailing: AnyView(
                    HStack(spacing: 8) {
                        Button("Open Receipt…") { engine.chooseFilePanel() }.controlSize(.small)
                        Button(action: triggerRun) {
                            Label("Verify Replay", systemImage: "checkmark.shield.fill")
                        }
                        .buttonStyle(.borderedProminent)
                        .controlSize(.small)
                        .keyboardShortcut(.return, modifiers: .command)
                        .disabled(engine.isBusy || !engine.isReady)
                    }
                )) {
                    TextEditor(text: $engine.verifyText)
                        .font(.system(size: 11, design: .monospaced))
                        .frame(height: 150)
                        .overlay(RoundedRectangle(cornerRadius: 4).stroke(Color(nsColor: .separatorColor)))
                }

                if let outcome = engine.receipts[.verify] {
                    let ok = outcome.flag("ok")
                    let checks = outcome.dictArr("checks")
                    LabSectionCard(title: ok ? "Verification Passed — Exact Cryptographic Replay" : "Verification Failed — Tampering or Model Mismatch Detected") {
                        VStack(alignment: .leading, spacing: 8) {
                            BadgePill(text: ok ? "ALL CHECKS PASSED" : "VERIFICATION FAILED", ok: ok)
                            Divider()
                            ForEach(Array(checks.enumerated()), id: \.offset) { _, ch in
                                HStack(spacing: 10) {
                                    Image(systemName: ch.flag("ok") ? "checkmark.circle.fill" : "xmark.octagon.fill")
                                        .foregroundStyle(ch.flag("ok") ? .green : .red)
                                    Text(ch.str("check")).font(.system(size: 12, weight: .medium, design: .monospaced))
                                    Spacer()
                                    Text(ch.flag("ok") ? "PASS" : "FAIL")
                                        .font(.system(size: 11, weight: .bold, design: .monospaced))
                                        .foregroundStyle(ch.flag("ok") ? .green : .red)
                                }
                            }
                        }
                    }
                }
            }
            .padding(14)
        }
    }
}
