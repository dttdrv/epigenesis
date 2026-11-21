import AppKit
import Combine
import Foundation
import SwiftUI

// MARK: - Navigation & Bilingual Strings (EN / BG)

enum WorkbenchTab: String, CaseIterable, Identifiable {
    case develop
    case inverse
    case simulate
    case circuit
    case scan
    case edits
    case verify

    var id: String { rawValue }

    var symbol: String {
        switch self {
        case .develop: return "point.3.connected.trianglepath.dotted"
        case .inverse: return "wand.and.stars"
        case .simulate: return "waveform.path.ecg"
        case .circuit: return "cpu"
        case .scan: return "square.grid.3x3.topleft.filled"
        case .edits: return "tablecells"
        case .verify: return "checkmark.shield"
        }
    }

    var section: String {
        switch self {
        case .develop, .simulate, .circuit: return "sec.dev"
        case .inverse: return "sec.design"
        case .scan, .edits: return "sec.screen"
        case .verify: return "sec.audit"
        }
    }
}

struct L10n {
    static let en: [String: String] = [
        "sec.dev": "DEVELOPMENTAL MODELS",
        "sec.design": "INVERSE DESIGN",
        "sec.screen": "VARIANT SCREENING",
        "sec.audit": "VERIFICATION",
        "sec.files": "WORKSPACE FILES",
        "sec.history": "SESSION RUNS",
        "tab.develop": "Develop Pipeline",
        "tab.inverse": "Inverse Compiler",
        "tab.simulate": "Regulatory Site",
        "tab.circuit": "Enhancer Circuit",
        "tab.scan": "Positional Scan",
        "tab.edits": "Batch Edits",
        "tab.verify": "Verify Receipt",
        "act.run": "Run Analysis",
        "act.open": "Open File…",
        "act.folder": "Workspace…",
        "act.exportCsv": "Export CSV…",
        "act.exportReceipt": "Save Receipt…",
        "act.copy": "Copy",
        "status.ready": "Engine Ready",
        "status.busy": "Computing…",
        "status.offline": "Engine Offline",
        "dev.sub": "Multi-locus DNA + GA4GH VRS 2.x → chromatin gating & coding translation → N-gene GRN uncertainty → 3D CX3D morphogenesis",
        "inv.sub": "Synthesize minimal compensatory 3D rescues, single-base fate edits, or run out-of-sample 532-TF human enhancer discovery",
        "sim.sub": "9-bp Nkx2.2 GLI regulatory site perturbation with Hallikas affinity profiles, multi-phase Shh signal protocol, and exact reversion",
        "cir.sub": "Strand-symmetric Boltzmann partition scanning across arbitrary human enhancers and N-gene circuits for cryptic fate bifurcations",
        "scan.sub": "Saturation mutagenesis across all 3L single-base substitutions with positional heatmap, 0–72 h temporal kinetics, and JASPAR TF attributions",
        "edits.sub": "Score arbitrary SNV, MNV, insertion, and deletion pairs with family-diverse panel selection and motif decomposition",
        "ver.sub": "Offline RFC 8785 SHA-256 digest check, model identity verification, and deterministic replay audit for saved receipts",
        "lbl.preset": "Program Preset",
        "lbl.accessibility": "Chromatin Accessibility",
        "lbl.samples": "Uncertainty Draws",
        "lbl.coding": "Include Coding Truncation",
        "lbl.cx3d": "Run 3D CX3D + LIF Assay",
        "lbl.model": "Predictor Model",
        "lbl.locus": "Human Enhancer Locus",
        "lbl.protein": "Transcription Factor",
        "lbl.circuit": "Regulatory Circuit",
        "lbl.targetGene": "Target Gene",
        "lbl.targetFate": "Target Cell Fate",
        "lbl.breakingEdit": "Locked Primary Lesion",
        "lbl.signal": "Morphogen Signal (S)",
        "lbl.budget": "Panel Budget",
        "lbl.noData": "Configure parameters above and press Run Analysis (⌘↩) to compute results.",
    ]

    static let bg: [String: String] = [
        "sec.dev": "МОДЕЛИ НА РАЗВИТИЕ",
        "sec.design": "ОБРАТЕН ДИЗАЙН",
        "sec.screen": "СКРИНИНГ НА ВАРИАНТИ",
        "sec.audit": "ВЕРИФИКАЦИЯ",
        "sec.files": "РАБОТНИ ФАЙЛОВЕ",
        "sec.history": "ИСТОРИЯ НА СЕСИЯТА",
        "tab.develop": "Конвейер за развитие",
        "tab.inverse": "Обратен компилатор",
        "tab.simulate": "Регулаторен сайт",
        "tab.circuit": "Верига на енхансер",
        "tab.scan": "Позиционно сканиране",
        "tab.edits": "Таблица с редакции",
        "tab.verify": "Проверка на удостоверение",
        "act.run": "Изпълни анализ",
        "act.open": "Отвори файл…",
        "act.folder": "Папка…",
        "act.exportCsv": "Експорт CSV…",
        "act.exportReceipt": "Запази удостоверение…",
        "act.copy": "Копирай",
        "status.ready": "Готовност",
        "status.busy": "Изчисляване…",
        "status.offline": "Извън линия",
        "dev.sub": "Многолокусна ДНК + GA4GH VRS 2.x → хроматинов вентил и кодираща транслация → N-генна ГРМ с 90% ДИ → 3D CX3D морфогенеза",
        "inv.sub": "Синтез на минимални епистатични 3D спасявания, еднобазово препрограмиране на съдба и 532-ТФ биофизично откриване",
        "sim.sub": "Пертурбация на 9-базов GLI сайт на Nkx2.2 с профили на Hallikas, многофазов Shh протокол и точно възстановяване",
        "cir.sub": "Болцманово сканиране на произволни човешки енхансери и N-генни вериги за скрити бифуркации и второсайтови спасявания",
        "scan.sub": "Насищащ мутагенезис по всички 3L единични замени с топлинна карта, 0–72 ч кинетика и приноси на JASPAR ТФ",
        "edits.sub": "Оценка на произволни SNV, MNV, инсерции и делеции с подбор на семейно-разнообразен панел и мотивен анализ",
        "ver.sub": "Офлайн проверка на RFC 8785 SHA-256 контролна сума, идентичност на модела и детерминистично преизчисляване",
        "lbl.preset": "Програмен профил",
        "lbl.accessibility": "Хроматинова достъпност",
        "lbl.samples": "Ансамблови проби",
        "lbl.coding": "Кодиращо скъсяване",
        "lbl.cx3d": "3D CX3D + LIF симулация",
        "lbl.model": "Модел предсказвател",
        "lbl.locus": "Човешки енхансер",
        "lbl.protein": "Транскрипционен фактор",
        "lbl.circuit": "Регулаторна верига",
        "lbl.targetGene": "Целеви ген",
        "lbl.targetFate": "Целева клетъчна съдба",
        "lbl.breakingEdit": "Заключена мутация",
        "lbl.signal": "Морфогенен сигнал (S)",
        "lbl.budget": "Размер на панела",
        "lbl.noData": "Задайте параметрите по-горе и натиснете Изпълни анализ (⌘↩).",
    ]

    static func t(_ key: String, lang: String) -> String {
        (lang == "bg" ? bg[key] : en[key]) ?? en[key] ?? key
    }
}

// MARK: - Dynamic JSON Access Helpers

typealias JSONDict = [String: Any]

extension Dictionary where Key == String, Value == Any {
    func str(_ k: String, _ fallback: String = "") -> String { (self[k] as? String) ?? fallback }
    func dbl(_ k: String, _ fallback: Double = 0) -> Double {
        if let d = self[k] as? Double { return d }
        if let i = self[k] as? Int { return Double(i) }
        return fallback
    }
    func int(_ k: String, _ fallback: Int = 0) -> Int {
        if let i = self[k] as? Int { return i }
        if let d = self[k] as? Double { return Int(d) }
        return fallback
    }
    func flag(_ k: String, _ fallback: Bool = false) -> Bool { (self[k] as? Bool) ?? fallback }
    func dict(_ k: String) -> JSONDict? { self[k] as? JSONDict }
    func arr(_ k: String) -> [Any] { (self[k] as? [Any]) ?? [] }
    func dictArr(_ k: String) -> [JSONDict] { (self[k] as? [JSONDict]) ?? [] }
    func dblArr(_ k: String) -> [Double] {
        arr(k).compactMap { ($0 as? Double) ?? ($0 as? Int).map(Double.init) }
    }
    func strArr(_ k: String) -> [String] { (self[k] as? [String]) ?? [] }
}

struct WorkspaceFile: Identifiable, Hashable {
    let id = UUID()
    let name: String
    let size: Int
    let kind: String
}

struct RunLogEntry: Identifiable {
    let id = UUID()
    let timestamp: Date
    let tab: WorkbenchTab
    let summary: String
    let latencyMs: Int
    let receipt: JSONDict
}

// MARK: - Native Stdio Pipe Engine Client (Zero Localhost / Zero HTTP)

final class EngineClient: ObservableObject {
    @Published var lang: String = "en"
    @Published var activeTab: WorkbenchTab = .develop
    @Published var isReady: Bool = false
    @Published var isBusy: Bool = false
    @Published var statusMessage: String = "Starting local compute kernel…"
    @Published var lastError: String? = nil
    @Published var lastLatencyMs: Int = 0
    @Published var engineVersion: String = "1.7.0"
    @Published var workspaceFolder: String = ""
    @Published var workspaceFiles: [WorkspaceFile] = []
    @Published var runHistory: [RunLogEntry] = []

    // Loaded metadata from kernel
    @Published var modelInfo: JSONDict = [:]
    @Published var siteTable: JSONDict = [:]
    @Published var circuitExample: JSONDict = [:]
    @Published var pipelinePresets: JSONDict = [:]
    @Published var exampleEdits: [JSONDict] = []

    // Latest receipt per tab
    @Published var receipts: [WorkbenchTab: JSONDict] = [:]

    // Active input buffers shared across workbenches (persisted across tab switches, zero SwiftUI macro dependency)
    @Published var scanSequence: String = ""
    @Published var circuitSequence: String = ""
    @Published var editsText: String = ""
    @Published var verifyText: String = ""

    // Develop Pipeline state
    @Published var devPreset: String = "neural-tube-multi-locus"
    @Published var devAccessibility: Double = 1.0
    @Published var devSamples: Int = 16
    @Published var devIncludeCoding: Bool = false
    @Published var devRunCx3d: Bool = false

    // Inverse Compiler state
    @Published var invGoal: String = "rescue"
    @Published var invLocus: String = "chr18:56475717-56475888"
    @Published var invBreakingEdit: String = "165:G>T"
    @Published var invTargetFate: String = "Olig2"
    @Published var invSignal: Double = 0.20
    @Published var invRunCx3d: Bool = false

    // Simulate Site state
    @Published var simSite: String = "GAACACCCA"
    @Published var simProtein: String = "GLI2"
    @Published var simDuration: Double = 20.0
    @Published var simSignal: Double = 0.20
    @Published var simHoverIdx: Int? = nil

    // Enhancer Circuit state
    @Published var cirName: String = "neural-tube"
    @Published var cirProtein: String = "GLI2"
    @Published var cirSignal: Double = 0.20

    // Positional Scan & Batch Edits state
    @Published var scanPredictor: String = "biophys"
    @Published var scanSelectedIdx: Int? = nil
    @Published var scanSearchText: String = ""
    @Published var editsPredictor: String = "biophys"
    @Published var editsBudget: Int = 5

    private var proc: Process?
    private var stdinPipe: Pipe?
    private var stdoutPipe: Pipe?
    private let ioQueue = DispatchQueue(label: "org.epigenesis.engine.io", qos: .userInitiated)
    private var nextId: Int = 1
    private var pending: [Int: (Result<JSONDict, Error>) -> Void] = [:]
    private let lock = NSLock()
    private var repoPath: String = ""

    init() {
        startKernel()
    }

    deinit {
        stopKernel()
    }

    func t(_ key: String) -> String { L10n.t(key, lang: lang) }

    func currentReceipt() -> JSONDict? { receipts[activeTab] }

    private func loadRuntimeConfig() -> (python: String, repo: String) {
        let fm = FileManager.default
        let env = ProcessInfo.processInfo.environment
        let home = fm.homeDirectoryForCurrentUser.path
        var py = ""
        var repo = ""

        if let cfgURL = Bundle.main.url(forResource: "runtime", withExtension: "json"),
           let data = try? Data(contentsOf: cfgURL),
           let obj = try? JSONSerialization.jsonObject(with: data) as? [String: String] {
            if let p = obj["python"], fm.isExecutableFile(atPath: p) { py = p }
            if let r = obj["repo"], fm.fileExists(atPath: "\(r)/desktop/epigenesis_engine.py") { repo = r }
        }
        if py.isEmpty {
            let candidates = [
                "\(home)/.local/share/uv/tools/epigenesis-studio/bin/python",
                repo.isEmpty ? "" : "\(repo)/.venv/bin/python",
                "/opt/homebrew/bin/python3",
                "/usr/local/bin/python3",
                "/usr/bin/python3",
            ]
            py = candidates.first(where: { !$0.isEmpty && fm.isExecutableFile(atPath: $0) }) ?? "/usr/bin/python3"
        }
        if repo.isEmpty {
            let repoCandidates = [
                env["EPIGENESIS_HOME"] ?? "",
                "\(home)/Projects/epigenesis",
                fm.currentDirectoryPath,
            ]
            repo = repoCandidates.first(where: { !$0.isEmpty && fm.fileExists(atPath: "\($0)/desktop/epigenesis_engine.py") }) ?? fm.currentDirectoryPath
        }
        return (py, repo)
    }

    func startKernel() {
        let cfg = loadRuntimeConfig()
        repoPath = cfg.repo
        let script = "\(cfg.repo)/desktop/epigenesis_engine.py"
        guard FileManager.default.fileExists(atPath: script) else {
            DispatchQueue.main.async {
                self.lastError = "Compute engine not found at \(script). Set EPIGENESIS_HOME or rebuild Epigenesis.app."
                self.statusMessage = "Engine missing"
            }
            return
        }

        let p = Process()
        p.executableURL = URL(fileURLWithPath: cfg.python)
        p.arguments = ["-u", script]
        p.currentDirectoryURL = URL(fileURLWithPath: cfg.repo)
        var env = ProcessInfo.processInfo.environment
        env["PYTHONUNBUFFERED"] = "1"
        env["PYTHONPATH"] = "\(cfg.repo):\(cfg.repo)/studio"
        p.environment = env

        let inPipe = Pipe()
        let outPipe = Pipe()
        let errPipe = Pipe()
        p.standardInput = inPipe
        p.standardOutput = outPipe
        p.standardError = errPipe

        do {
            try p.run()
        } catch {
            DispatchQueue.main.async {
                self.lastError = "Failed to launch engine: \(error.localizedDescription)"
            }
            return
        }
        self.proc = p
        self.stdinPipe = inPipe
        self.stdoutPipe = outPipe

        ioQueue.async { [weak self] in
            self?.readLoop(outPipe: outPipe, errPipe: errPipe)
        }
    }

    func stopKernel() {
        stdinPipe?.fileHandleForWriting.closeFile()
        if let p = proc, p.isRunning {
            p.terminate()
        }
        proc = nil
    }

    private func readLoop(outPipe: Pipe, errPipe: Pipe) {
        let handle = outPipe.fileHandleForReading
        var buffer = Data()
        while true {
            let chunk = handle.availableData
            if chunk.isEmpty { break }
            buffer.append(chunk)
            while let nl = buffer.firstIndex(of: 0x0A) {
                let lineData = buffer.subdata(in: 0..<nl)
                buffer.removeSubrange(0...nl)
                guard !lineData.isEmpty,
                      let obj = (try? JSONSerialization.jsonObject(with: lineData)) as? JSONDict else {
                    continue
                }
                if obj["ready"] as? Bool == true {
                    let ver = obj.str("version", "0.3.0")
                    DispatchQueue.main.async {
                        self.isReady = true
                        self.engineVersion = ver
                        self.statusMessage = self.t("status.ready")
                        self.loadInitialInfo()
                    }
                    continue
                }
                if let reqId = obj["id"] as? Int {
                    lock.lock()
                    let cb = pending.removeValue(forKey: reqId)
                    lock.unlock()
                    if obj.flag("ok"), let data = obj.dict("data") {
                        cb?(.success(data))
                    } else {
                        let msg = obj.str("error", "Unknown engine error")
                        cb?(.failure(NSError(domain: "EpigenesisEngine", code: 1, userInfo: [NSLocalizedDescriptionKey: msg])))
                    }
                }
            }
        }
        let errText = String(data: errPipe.fileHandleForReading.readDataToEndOfFile(), encoding: .utf8) ?? ""
        DispatchQueue.main.async {
            self.isReady = false
            if !errText.isEmpty {
                self.lastError = "Engine exited: \(errText.prefix(500))"
            }
        }
    }

    func call(op: String, payload: JSONDict = [:], completion: @escaping (Result<JSONDict, Error>) -> Void) {
        lock.lock()
        let reqId = nextId
        nextId += 1
        pending[reqId] = completion
        lock.unlock()

        var msg = payload
        msg["id"] = reqId
        msg["op"] = op
        guard let data = try? JSONSerialization.data(withJSONObject: msg),
              let handle = stdinPipe?.fileHandleForWriting else {
            completion(.failure(NSError(domain: "EpigenesisEngine", code: 2, userInfo: [NSLocalizedDescriptionKey: "Engine pipe unavailable"])))
            return
        }
        var line = data
        line.append(0x0A)
        ioQueue.async {
            handle.write(line)
        }
    }

    func loadInitialInfo() {
        call(op: "info") { [weak self] res in
            DispatchQueue.main.async {
                guard let self else { return }
                switch res {
                case .success(let info):
                    self.modelInfo = info
                    self.workspaceFolder = info.str("folder")
                    self.workspaceFiles = info.dictArr("files").map {
                        WorkspaceFile(name: $0.str("name"), size: $0.int("size"), kind: $0.str("kind"))
                    }
                    if let st = info.dict("site_table") { self.siteTable = st }
                    if let ce = info.dict("circuit_example") {
                        self.circuitExample = ce
                        if self.circuitSequence.isEmpty {
                            self.circuitSequence = ce.str("sequence")
                        }
                    }
                    if let pp = info.dict("pipeline_presets") { self.pipelinePresets = pp }
                    let edits = info.dictArr("example_edits")
                    self.exampleEdits = edits
                    if self.scanSequence.isEmpty, let first = edits.first {
                        self.scanSequence = first.str("reference")
                    }
                    if self.editsText.isEmpty, !edits.isEmpty,
                       let data = try? JSONSerialization.data(withJSONObject: edits, options: [.prettyPrinted]),
                       let str = String(data: data, encoding: .utf8) {
                        self.editsText = str
                    }
                case .failure(let err):
                    self.lastError = err.localizedDescription
                }
            }
        }
    }

    func runScore(tab: WorkbenchTab, request: JSONDict, summaryLabel: String, onDone: ((JSONDict) -> Void)? = nil) {
        let t0 = Date()
        DispatchQueue.main.async {
            self.isBusy = true
            self.lastError = nil
            self.statusMessage = self.t("status.busy")
        }
        call(op: "score", payload: ["request": request]) { [weak self] res in
            let ms = Int(Date().timeIntervalSince(t0) * 1000)
            DispatchQueue.main.async {
                guard let self else { return }
                self.isBusy = false
                self.lastLatencyMs = ms
                self.statusMessage = self.t("status.ready")
                switch res {
                case .success(let receipt):
                    self.receipts[tab] = receipt
                    self.runHistory.insert(
                        RunLogEntry(timestamp: Date(), tab: tab, summary: summaryLabel, latencyMs: ms, receipt: receipt),
                        at: 0
                    )
                    if self.runHistory.count > 40 { self.runHistory.removeLast() }
                    onDone?(receipt)
                case .failure(let err):
                    self.lastError = err.localizedDescription
                }
            }
        }
    }

    func runVerify(receiptObj: JSONDict, onDone: ((JSONDict) -> Void)? = nil) {
        let t0 = Date()
        DispatchQueue.main.async {
            self.isBusy = true
            self.lastError = nil
            self.statusMessage = self.t("status.busy")
        }
        call(op: "verify", payload: ["receipt": receiptObj]) { [weak self] res in
            let ms = Int(Date().timeIntervalSince(t0) * 1000)
            DispatchQueue.main.async {
                guard let self else { return }
                self.isBusy = false
                self.lastLatencyMs = ms
                self.statusMessage = self.t("status.ready")
                switch res {
                case .success(let outcome):
                    self.receipts[.verify] = outcome
                    onDone?(outcome)
                case .failure(let err):
                    self.lastError = err.localizedDescription
                }
            }
        }
    }

    func openDocument(path: String) {
        call(op: "read_file", payload: ["path": path]) { [weak self] res in
            DispatchQueue.main.async {
                guard let self else { return }
                switch res {
                case .success(let doc):
                    let kind = doc.str("kind")
                    let text = doc.str("text")
                    if kind == "receipt" {
                        self.verifyText = text
                        self.activeTab = .verify
                        if let data = text.data(using: .utf8),
                           let obj = (try? JSONSerialization.jsonObject(with: data)) as? JSONDict {
                            self.runVerify(receiptObj: obj)
                        }
                    } else if kind == "edits" {
                        self.editsText = text
                        self.activeTab = .edits
                    } else {
                        self.scanSequence = text
                        self.circuitSequence = text
                        self.activeTab = .scan
                    }
                case .failure(let err):
                    self.lastError = err.localizedDescription
                }
            }
        }
    }

    func chooseFilePanel() {
        let panel = NSOpenPanel()
        panel.canChooseFiles = true
        panel.canChooseDirectories = false
        panel.allowsMultipleSelection = false
        panel.begin { resp in
            guard resp == .OK, let url = panel.url else { return }
            self.openDocument(path: url.path)
        }
    }

    func chooseFolderPanel() {
        let panel = NSOpenPanel()
        panel.canChooseFiles = false
        panel.canChooseDirectories = true
        panel.begin { resp in
            guard resp == .OK, let url = panel.url else { return }
            self.call(op: "set_folder", payload: ["folder": url.path]) { res in
                DispatchQueue.main.async {
                    if case .success(let info) = res {
                        self.workspaceFolder = info.str("folder")
                        self.workspaceFiles = info.dictArr("files").map {
                            WorkspaceFile(name: $0.str("name"), size: $0.int("size"), kind: $0.str("kind"))
                        }
                    }
                }
            }
        }
    }

    func exportCurrentReceipt() {
        guard let receipt = currentReceipt(),
              let data = try? JSONSerialization.data(withJSONObject: receipt, options: [.prettyPrinted, .sortedKeys]),
              let text = String(data: data, encoding: .utf8) else { return }
        let sha = receipt.str("receipt_sha256")
        let tag = sha.isEmpty ? activeTab.rawValue : String(sha.prefix(8))
        saveTextFile(defaultName: "epigenesis-receipt-\(tag).json", text: text + "\n")
    }

    func saveTextFile(defaultName: String, text: String) {
        let panel = NSSavePanel()
        panel.nameFieldStringValue = defaultName
        panel.begin { resp in
            guard resp == .OK, let dest = panel.url else { return }
            try? text.write(to: dest, atomically: true, encoding: .utf8)
        }
    }
}
