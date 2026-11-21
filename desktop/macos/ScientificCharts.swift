import AppKit
import SwiftUI

// MARK: - Laboratory UI Building Blocks (SnapGene / GraphPad Prism / Xcode style)

struct LabSectionCard<Content: View>: View {
    let title: String
    var subtitle: String? = nil
    var trailing: AnyView? = nil
    @ViewBuilder let content: () -> Content

    var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            HStack(alignment: .firstTextBaseline, spacing: 8) {
                VStack(alignment: .leading, spacing: 2) {
                    Text(title)
                        .font(.system(size: 12.5, weight: .semibold))
                        .foregroundStyle(.primary)
                    if let sub = subtitle, !sub.isEmpty {
                        Text(sub)
                            .font(.system(size: 11))
                            .foregroundStyle(.secondary)
                    }
                }
                Spacer()
                if let trailing { trailing }
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 8)
            .background(Color(nsColor: .controlBackgroundColor))

            Divider()

            content()
                .padding(12)
        }
        .background(Color(nsColor: .textBackgroundColor))
        .clipShape(RoundedRectangle(cornerRadius: 6, style: .continuous))
        .overlay(
            RoundedRectangle(cornerRadius: 6, style: .continuous)
                .strokeBorder(Color(nsColor: .separatorColor), lineWidth: 1)
        )
    }
}

struct MetricTile: View {
    let label: String
    let value: String
    var detail: String = ""
    var tone: Color = .primary

    var body: some View {
        VStack(alignment: .leading, spacing: 3) {
            Text(label)
                .font(.system(size: 10.5, weight: .medium))
                .foregroundStyle(.secondary)
                .textCase(.uppercase)
            Text(value)
                .font(.system(size: 18, weight: .semibold, design: .monospaced))
                .monospacedDigit()
                .foregroundStyle(tone)
                .lineLimit(1)
                .minimumScaleFactor(0.75)
            if !detail.isEmpty {
                Text(detail)
                    .font(.system(size: 10.5))
                    .foregroundStyle(.secondary)
                    .lineLimit(1)
            }
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(.horizontal, 10)
        .padding(.vertical, 8)
        .background(Color(nsColor: .controlBackgroundColor).opacity(0.6))
        .clipShape(RoundedRectangle(cornerRadius: 5, style: .continuous))
        .overlay(
            RoundedRectangle(cornerRadius: 5, style: .continuous)
                .strokeBorder(Color(nsColor: .separatorColor).opacity(0.8), lineWidth: 1)
        )
    }
}

struct BadgePill: View {
    let text: String
    let ok: Bool

    var body: some View {
        HStack(spacing: 4) {
            Image(systemName: ok ? "checkmark.circle.fill" : "exclamationmark.triangle.fill")
                .font(.system(size: 10))
            Text(text)
                .font(.system(size: 10.5, weight: .medium, design: .monospaced))
        }
        .padding(.horizontal, 7)
        .padding(.vertical, 3)
        .foregroundStyle(ok ? Color.green : Color.orange)
        .background((ok ? Color.green : Color.orange).opacity(0.12))
        .clipShape(RoundedRectangle(cornerRadius: 4, style: .continuous))
    }
}

// MARK: - SnapGene / IGV Style Interactive DNA Sequence Ruler & Track

struct DNASequenceTrackView: View {
    let sequence: String
    var highlightPos1Based: Int? = nil
    var altBase: String? = nil
    var onSelectPos: ((Int) -> Void)? = nil

    private func baseColor(_ ch: Character) -> Color {
        switch ch {
        case "A": return Color(red: 0.13, green: 0.56, blue: 0.28)
        case "C": return Color(red: 0.15, green: 0.42, blue: 0.80)
        case "G": return Color(red: 0.78, green: 0.50, blue: 0.06)
        case "T": return Color(red: 0.78, green: 0.20, blue: 0.20)
        default: return .secondary
        }
    }

    var body: some View {
        let chars = Array(sequence.uppercased().filter { "ACGTN".contains($0) })
        ScrollView(.horizontal, showsIndicators: true) {
            VStack(alignment: .leading, spacing: 2) {
                // 10-bp coordinate ruler
                HStack(spacing: 1) {
                    ForEach(chars.indices, id: \.self) { i in
                        let pos = i + 1
                        Text(pos == 1 || pos % 10 == 0 ? "\(pos)" : (pos % 5 == 0 ? "·" : ""))
                            .font(.system(size: 8.5, design: .monospaced))
                            .foregroundStyle(.secondary)
                            .frame(width: 14, alignment: .leading)
                    }
                }
                // Nucleotide track
                HStack(spacing: 1) {
                    ForEach(chars.indices, id: \.self) { i in
                        let pos = i + 1
                        let ch = chars[i]
                        let isHi = (highlightPos1Based == pos)
                        Button {
                            onSelectPos?(pos)
                        } label: {
                            Text(String(ch))
                                .font(.system(size: 11, weight: isHi ? .bold : .medium, design: .monospaced))
                                .frame(width: 14, height: 20)
                                .foregroundStyle(isHi ? .white : baseColor(ch))
                                .background(isHi ? Color.accentColor : baseColor(ch).opacity(0.11))
                                .clipShape(RoundedRectangle(cornerRadius: 2))
                        }
                        .buttonStyle(.plain)
                        .help("Position \(pos): \(String(ch))")
                    }
                }
                if let hp = highlightPos1Based, let alt = altBase, hp >= 1, hp <= chars.count {
                    HStack(spacing: 6) {
                        Text("Selected coordinate \(hp): \(String(chars[hp - 1])) → \(alt)")
                            .font(.system(size: 10.5, weight: .medium, design: .monospaced))
                            .foregroundStyle(Color.accentColor)
                    }
                    .padding(.top, 3)
                }
            }
            .padding(.vertical, 4)
        }
    }
}

// MARK: - Publication-Grade Multi-Gene ODE Trajectory Plot (with Hover Crosshair)

struct ODETrajectoryPlotView: View {
    let times: [Double]
    let genes: [String]
    let refStates: [[Double]]
    let varStates: [[Double]]
    var refOcc: [Double] = []
    var varOcc: [Double] = []
    var repOcc: [Double] = []
    @Binding var hoverIdx: Int?

    private let palette: [Color] = [
        Color(red: 0.14, green: 0.45, blue: 0.78),
        Color(red: 0.85, green: 0.37, blue: 0.05),
        Color(red: 0.08, green: 0.58, blue: 0.42),
        Color(red: 0.56, green: 0.27, blue: 0.68),
    ]

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            // Legend + hover readout bar
            HStack(spacing: 14) {
                ForEach(Array(genes.enumerated()), id: \.offset) { gi, g in
                    HStack(spacing: 5) {
                        RoundedRectangle(cornerRadius: 1.5)
                            .fill(palette[gi % palette.count])
                            .frame(width: 12, height: 3)
                        Text(g)
                            .font(.system(size: 11, weight: .medium, design: .monospaced))
                    }
                }
                Spacer()
                Text("Solid: WT   Dashed: Variant")
                    .font(.system(size: 10.5))
                    .foregroundStyle(.secondary)
            }

            GeometryReader { geo in
                let padL: CGFloat = 38
                let padR: CGFloat = 12
                let padT: CGFloat = 10
                let padB: CGFloat = 24
                let w = max(10, geo.size.width - padL - padR)
                let h = max(10, geo.size.height - padT - padB)
                let maxT = max(times.last ?? 1.0, 1e-6)
                let maxY = max(1.0, (refStates.flatMap { $0 } + varStates.flatMap { $0 }).max() ?? 1.0)

                ZStack(alignment: .topLeading) {
                    Canvas { ctx, size in
                        // Gridlines & Y ticks
                        for tick in stride(from: 0.0, through: 1.0, by: 0.25) {
                            let y = padT + h * CGFloat(1.0 - tick)
                            var grid = Path()
                            grid.move(to: CGPoint(x: padL, y: y))
                            grid.addLine(to: CGPoint(x: padL + w, y: y))
                            ctx.stroke(grid, with: .color(.gray.opacity(0.22)), lineWidth: 0.8)
                            let label = Text(String(format: "%.2f", tick * maxY))
                                .font(.system(size: 9.5, design: .monospaced))
                                .foregroundStyle(.secondary)
                            ctx.draw(label, at: CGPoint(x: padL - 6, y: y), anchor: .trailing)
                        }
                        // X axis & ticks
                        for frac in stride(from: 0.0, through: 1.0, by: 0.25) {
                            let x = padL + w * CGFloat(frac)
                            let label = Text(String(format: "%.1f", frac * maxT))
                                .font(.system(size: 9.5, design: .monospaced))
                                .foregroundStyle(.secondary)
                            ctx.draw(label, at: CGPoint(x: x, y: padT + h + 11), anchor: .center)
                        }

                        // Plot each gene (WT solid, Variant dashed)
                        for gi in 0..<genes.count {
                            let col = palette[gi % palette.count]
                            var pRef = Path()
                            var pVar = Path()
                            for (ti, t) in times.enumerated() {
                                let x = padL + w * CGFloat(t / maxT)
                                if ti < refStates.count, gi < refStates[ti].count {
                                    let y = padT + h * CGFloat(1.0 - min(maxY, max(0, refStates[ti][gi])) / maxY)
                                    if ti == 0 { pRef.move(to: CGPoint(x: x, y: y)) } else { pRef.addLine(to: CGPoint(x: x, y: y)) }
                                }
                                if ti < varStates.count, gi < varStates[ti].count {
                                    let y = padT + h * CGFloat(1.0 - min(maxY, max(0, varStates[ti][gi])) / maxY)
                                    if ti == 0 { pVar.move(to: CGPoint(x: x, y: y)) } else { pVar.addLine(to: CGPoint(x: x, y: y)) }
                                }
                            }
                            ctx.stroke(pRef, with: .color(col.opacity(0.45)), lineWidth: 1.6)
                            ctx.stroke(pVar, with: .color(col), style: StrokeStyle(lineWidth: 2.0, dash: [5, 3]))
                        }

                        // Hover crosshair
                        if let hi = hoverIdx, hi < times.count {
                            let x = padL + w * CGFloat(times[hi] / maxT)
                            var ch = Path()
                            ch.move(to: CGPoint(x: x, y: padT))
                            ch.addLine(to: CGPoint(x: x, y: padT + h))
                            ctx.stroke(ch, with: .color(.primary.opacity(0.35)), lineWidth: 1)
                        }
                    }
                    .onContinuousHover { phase in
                        switch phase {
                        case .active(let loc):
                            let frac = max(0, min(1, Double((loc.x - padL) / w)))
                            let targetT = frac * maxT
                            if let idx = times.enumerated().min(by: { abs($0.element - targetT) < abs($1.element - targetT) })?.offset {
                                hoverIdx = idx
                            }
                        case .ended:
                            hoverIdx = nil
                        }
                    }
                }
            }
            .frame(height: 210)

            if let hi = hoverIdx, hi < times.count, hi < varStates.count {
                HStack(spacing: 12) {
                    Text(String(format: "t = %.2f", times[hi]))
                        .font(.system(size: 10.5, weight: .semibold, design: .monospaced))
                    ForEach(Array(genes.enumerated()), id: \.offset) { gi, g in
                        let r = gi < refStates[hi].count ? refStates[hi][gi] : 0
                        let v = gi < varStates[hi].count ? varStates[hi][gi] : 0
                        Text("\(g): \(String(format: "%.3f", r)) → \(String(format: "%.3f", v))")
                            .font(.system(size: 10.5, design: .monospaced))
                            .foregroundStyle(palette[gi % palette.count])
                    }
                }
                .padding(.horizontal, 8)
                .padding(.vertical, 4)
                .background(Color(nsColor: .controlBackgroundColor))
                .clipShape(RoundedRectangle(cornerRadius: 4))
            }
        }
    }
}

// MARK: - GRN Steady-State + 90% Parameter Uncertainty Interval Plot

struct GRNUncertaintyChartView: View {
    struct Row: Identifiable {
        let id: String
        let wt: Double
        let variant: Double
        let deltaMean: Double
        let p05: Double
        let p95: Double
    }
    let rows: [Row]

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            ForEach(rows) { row in
                HStack(spacing: 10) {
                    Text(row.id)
                        .font(.system(size: 11.5, weight: .semibold, design: .monospaced))
                        .frame(width: 82, alignment: .leading)

                    GeometryReader { geo in
                        let w = geo.size.width
                        ZStack(alignment: .leading) {
                            RoundedRectangle(cornerRadius: 3)
                                .fill(Color.gray.opacity(0.14))
                                .frame(height: 16)
                            // WT bar (muted slate)
                            RoundedRectangle(cornerRadius: 2)
                                .fill(Color.blue.opacity(0.32))
                                .frame(width: max(2, w * CGFloat(min(1, max(0, row.wt)))), height: 7)
                                .offset(y: -4)
                            // Variant bar (crisp orange/blue)
                            RoundedRectangle(cornerRadius: 2)
                                .fill(row.variant >= row.wt ? Color.orange : Color.blue)
                                .frame(width: max(2, w * CGFloat(min(1, max(0, row.variant)))), height: 7)
                                .offset(y: 4)
                            // 90% CI whisker around WT + delta
                            let lo = max(0, min(1, row.wt + row.p05))
                            let hi = max(0, min(1, row.wt + row.p95))
                            if hi > lo {
                                Rectangle()
                                    .fill(Color.primary.opacity(0.65))
                                    .frame(width: max(2, w * CGFloat(hi - lo)), height: 1.5)
                                    .offset(x: w * CGFloat(lo), y: 4)
                            }
                        }
                    }
                    .frame(height: 18)

                    Text(String(format: "WT %.3f → Var %.3f  (90%% CI [%+.3f, %+.3f])", row.wt, row.variant, row.p05, row.p95))
                        .font(.system(size: 10.5, design: .monospaced))
                        .monospacedDigit()
                        .foregroundStyle(.secondary)
                        .frame(width: 280, alignment: .trailing)
                }
            }
        }
    }
}

// MARK: - 0-72h Neural Induction Temporal Kinetics Curve

struct TemporalCurvePlotView: View {
    let times: [Double]
    let values: [Double]

    var body: some View {
        GeometryReader { geo in
            let padL: CGFloat = 38
            let padR: CGFloat = 14
            let padT: CGFloat = 10
            let padB: CGFloat = 22
            let w = max(10, geo.size.width - padL - padR)
            let h = max(10, geo.size.height - padT - padB)
            let maxAbs = max(0.25, (values.map { abs($0) }.max() ?? 0.5) * 1.15)
            let maxT = max(times.last ?? 72.0, 1.0)

            Canvas { ctx, _ in
                let yZero = padT + h * 0.5
                var zeroLine = Path()
                zeroLine.move(to: CGPoint(x: padL, y: yZero))
                zeroLine.addLine(to: CGPoint(x: padL + w, y: yZero))
                ctx.stroke(zeroLine, with: .color(.secondary.opacity(0.45)), style: StrokeStyle(lineWidth: 1, dash: [3, 3]))

                for (lbl, yPos) in [(String(format: "+%.2f", maxAbs), padT), ("0.00", yZero), (String(format: "-%.2f", maxAbs), padT + h)] {
                    ctx.draw(Text(lbl).font(.system(size: 9, design: .monospaced)).foregroundStyle(.secondary),
                             at: CGPoint(x: padL - 5, y: yPos), anchor: .trailing)
                }

                var curve = Path()
                for (i, t) in times.enumerated() where i < values.count {
                    let x = padL + w * CGFloat(t / maxT)
                    let y = padT + h * CGFloat(0.5 - (values[i] / (2.0 * maxAbs)))
                    if i == 0 { curve.move(to: CGPoint(x: x, y: y)) } else { curve.addLine(to: CGPoint(x: x, y: y)) }
                    let dot = CGRect(x: x - 3.5, y: y - 3.5, width: 7, height: 7)
                    ctx.fill(Path(ellipseIn: dot), with: .color(values[i] >= 0 ? .orange : .blue))
                    ctx.draw(Text("\(Int(t))h").font(.system(size: 9, design: .monospaced)).foregroundStyle(.secondary),
                             at: CGPoint(x: x, y: padT + h + 10), anchor: .center)
                }
                ctx.stroke(curve, with: .color(.accentColor), lineWidth: 2)
            }
        }
        .frame(height: 150)
    }
}

// MARK: - Clickable Retina Saturation Mutagenesis Heatmap

struct SaturationHeatmapView: View {
    let sequence: String
    let predictions: [JSONDict] // [{id, position, ref, alt, effect, ...}]
    @Binding var selectedIndex: Int?

    private let bases = ["A", "C", "G", "T"]

    private func parseEditId(_ id: String) -> (pos: Int, alt: String)? {
        let parts = id.split(separator: ":")
        guard parts.count == 2, let pos = Int(parts[0]) else { return nil }
        let change = parts[1].split(separator: ">")
        guard change.count == 2 else { return nil }
        return (pos, String(change[1]))
    }

    var body: some View {
        let len = sequence.count
        let scale = max(0.2, predictions.map { abs($0.dbl("effect")) }.max() ?? 0.5)
        let lookup: [String: (Int, Double)] = Dictionary(uniqueKeysWithValues: predictions.enumerated().compactMap { idx, p in
            let parsed = parseEditId(p.str("id"))
            let pos = p.int("position", parsed?.pos ?? 0)
            let alt = p.str("alt", parsed?.alt ?? "")
            guard pos > 0, !alt.isEmpty else { return nil }
            return ("\(pos):\(alt)", (idx, p.dbl("effect")))
        })
        let cellW: CGFloat = 10
        let cellH: CGFloat = 18

        VStack(alignment: .leading, spacing: 6) {
            ScrollView(.horizontal, showsIndicators: true) {
                HStack(alignment: .top, spacing: 4) {
                    VStack(spacing: 2) {
                        Text(" ").font(.system(size: 9)).frame(height: 14)
                        ForEach(bases, id: \.self) { b in
                            Text(b)
                                .font(.system(size: 10, weight: .semibold, design: .monospaced))
                                .frame(width: 14, height: cellH)
                        }
                    }

                    Canvas { ctx, _ in
                        // Ruler ticks
                        for pos in stride(from: 1, through: len, by: 10) {
                            let x = CGFloat(pos - 1) * (cellW + 1)
                            ctx.draw(Text("\(pos)").font(.system(size: 8.5, design: .monospaced)).foregroundStyle(.secondary),
                                     at: CGPoint(x: x + 2, y: 6), anchor: .leading)
                        }
                        // Cells
                        for pos in 1...max(1, len) {
                            let x = CGFloat(pos - 1) * (cellW + 1)
                            for (bi, alt) in bases.enumerated() {
                                let y = 16 + CGFloat(bi) * (cellH + 2)
                                let rect = CGRect(x: x, y: y, width: cellW, height: cellH)
                                if let (idx, eff) = lookup["\(pos):\(alt)"] {
                                    let mag = min(1.0, abs(eff) / scale)
                                    let col: Color = eff >= 0
                                        ? Color(red: 0.90, green: 0.42, blue: 0.08).opacity(0.12 + 0.85 * mag)
                                        : Color(red: 0.16, green: 0.44, blue: 0.86).opacity(0.12 + 0.85 * mag)
                                    ctx.fill(Path(roundedRect: rect, cornerRadius: 1.5), with: .color(col))
                                    if selectedIndex == idx {
                                        ctx.stroke(Path(roundedRect: rect, cornerRadius: 1.5), with: .color(.primary), lineWidth: 1.8)
                                    }
                                } else {
                                    ctx.fill(Path(roundedRect: rect, cornerRadius: 1.5), with: .color(.gray.opacity(0.12)))
                                }
                            }
                        }
                    }
                    .frame(width: max(200, CGFloat(len) * (cellW + 1)), height: 16 + 4 * (cellH + 2))
                    .contentShape(Rectangle())
                    .gesture(
                        DragGesture(minimumDistance: 0).onEnded { value in
                            let pos = Int(value.location.x / (cellW + 1)) + 1
                            let bi = Int((value.location.y - 16) / (cellH + 2))
                            guard pos >= 1, pos <= len, bi >= 0, bi < 4 else { return }
                            if let (idx, _) = lookup["\(pos):\(bases[bi])"] {
                                selectedIndex = idx
                            }
                        }
                    )
                }
            }
            HStack(spacing: 10) {
                Text("Repress (−\(String(format: "%.2f", scale)))").font(.system(size: 10, design: .monospaced)).foregroundStyle(.blue)
                LinearGradient(colors: [.blue, .gray.opacity(0.2), .orange], startPoint: .leading, endPoint: .trailing)
                    .frame(width: 110, height: 7)
                    .clipShape(Capsule())
                Text("Activate (+\(String(format: "%.2f", scale)))").font(.system(size: 10, design: .monospaced)).foregroundStyle(.orange)
                Spacer()
                Text("Click any base cell to inspect motif & TF drivers").font(.system(size: 10.5)).foregroundStyle(.secondary)
            }
        }
    }
}
