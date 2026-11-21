import AppKit
import SwiftUI

// MARK: - Native macOS Desktop Workstation Shell (Zero WKWebView / Zero Localhost)

struct MainWorkbenchView: View {
    @ObservedObject var engine: EngineClient

    private let sections: [(String, [WorkbenchTab])] = [
        ("sec.dev", [.develop, .simulate, .circuit]),
        ("sec.design", [.inverse]),
        ("sec.screen", [.scan, .edits]),
        ("sec.audit", [.verify]),
    ]

    var body: some View {
        VStack(spacing: 0) {
            NavigationSplitView {
                List(selection: Binding<WorkbenchTab?>(
                    get: { engine.activeTab },
                    set: { if let v = $0 { engine.activeTab = v } }
                )) {
                    ForEach(sections, id: \.0) { secKey, tabs in
                        Section(engine.t(secKey)) {
                            ForEach(tabs) { tab in
                                Label(engine.t("tab.\(tab.rawValue)"), systemImage: tab.symbol)
                                    .tag(tab)
                            }
                        }
                    }

                    if !engine.workspaceFiles.isEmpty {
                        Section(engine.t("sec.files")) {
                            ForEach(engine.workspaceFiles.prefix(12)) { f in
                                Button {
                                    let full = (engine.workspaceFolder as NSString).appendingPathComponent(f.name)
                                    engine.openDocument(path: full)
                                } label: {
                                    HStack(spacing: 6) {
                                        Image(systemName: f.kind == "receipt" ? "checkmark.seal" : (f.kind == "edits" ? "tablecells" : "doc.text"))
                                            .font(.system(size: 11))
                                            .foregroundStyle(.secondary)
                                        Text(f.name)
                                            .font(.system(size: 11.5))
                                            .lineLimit(1)
                                            .truncationMode(.middle)
                                    }
                                }
                                .buttonStyle(.plain)
                            }
                        }
                    }

                    if !engine.runHistory.isEmpty {
                        Section(engine.t("sec.history")) {
                            ForEach(engine.runHistory.prefix(10)) { item in
                                Button {
                                    engine.activeTab = item.tab
                                    engine.receipts[item.tab] = item.receipt
                                } label: {
                                    VStack(alignment: .leading, spacing: 1) {
                                        Text(item.summary)
                                            .font(.system(size: 11, weight: .medium))
                                            .lineLimit(1)
                                        Text("\(item.latencyMs) ms")
                                            .font(.system(size: 9.5, design: .monospaced))
                                            .foregroundStyle(.secondary)
                                    }
                                }
                                .buttonStyle(.plain)
                            }
                        }
                    }
                }
                .listStyle(.sidebar)
                .navigationSplitViewColumnWidth(min: 215, ideal: 238, max: 290)
            } detail: {
                VStack(spacing: 0) {
                    if let err = engine.lastError {
                        HStack(spacing: 8) {
                            Image(systemName: "exclamationmark.triangle.fill").foregroundStyle(.red)
                            Text(err).font(.system(size: 11.5)).lineLimit(2)
                            Spacer()
                            Button("Dismiss") { engine.lastError = nil }
                                .buttonStyle(.borderless)
                                .font(.system(size: 11))
                        }
                        .padding(.horizontal, 14)
                        .padding(.vertical, 6)
                        .background(Color.red.opacity(0.10))
                        Divider()
                    }

                    Group {
                        switch engine.activeTab {
                        case .develop: DevelopView(engine: engine)
                        case .inverse: InverseCompilerView(engine: engine)
                        case .simulate: SimulateSiteView(engine: engine)
                        case .circuit: EnhancerCircuitView(engine: engine)
                        case .scan: PositionalScanView(engine: engine)
                        case .edits: BatchEditsView(engine: engine)
                        case .verify: VerifyReceiptView(engine: engine)
                        }
                    }
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
                    .background(Color(nsColor: .windowBackgroundColor))
                }
            }
            .toolbar {
                ToolbarItemGroup(placement: .navigation) {
                    Button { engine.chooseFilePanel() } label: {
                        Label(engine.t("act.open"), systemImage: "doc.badge.plus")
                    }
                    .help("Open DNA FASTA, CSV/JSON edits, or receipt (⌘O)")

                    Button { engine.chooseFolderPanel() } label: {
                        Label(engine.t("act.folder"), systemImage: "folder")
                    }
                    .help("Choose active workspace directory")
                }

                ToolbarItemGroup(placement: .primaryAction) {
                    if engine.isBusy {
                        ProgressView().controlSize(.small)
                    }

                    Button { engine.exportCurrentReceipt() } label: {
                        Label(engine.t("act.exportReceipt"), systemImage: "square.and.arrow.up")
                    }
                    .disabled(engine.currentReceipt() == nil)
                    .help("Save RFC 8785 cryptographic receipt JSON (⌘S)")

                    Picker("Language", selection: $engine.lang) {
                        Text("EN").tag("en")
                        Text("BG").tag("bg")
                    }
                    .pickerStyle(.segmented)
                    .frame(width: 78)
                }
            }

            Divider()

            // Classic Laboratory / Office Desktop Status Bar
            HStack(spacing: 14) {
                HStack(spacing: 5) {
                    Circle()
                        .fill(engine.isReady ? (engine.isBusy ? Color.orange : Color.green) : Color.red)
                        .frame(width: 7, height: 7)
                    Text(engine.statusMessage)
                        .font(.system(size: 11, weight: .medium))
                }
                Divider().frame(height: 12)
                Text("Native Pipe Kernel · No Network Ports")
                    .font(.system(size: 10.5, design: .monospaced))
                    .foregroundStyle(.secondary)
                if engine.lastLatencyMs > 0 {
                    Divider().frame(height: 12)
                    Text("Last run: \(engine.lastLatencyMs) ms")
                        .font(.system(size: 10.5, design: .monospaced))
                        .monospacedDigit()
                        .foregroundStyle(.secondary)
                }
                Spacer()
                Text("⌘O Open  ·  ⌘↩ Run  ·  ⌘S Save Receipt")
                    .font(.system(size: 10.5, design: .monospaced))
                    .foregroundStyle(.secondary)
            }
            .padding(.horizontal, 12)
            .frame(height: 24)
            .background(Color(nsColor: .controlBackgroundColor))
        }
        .frame(minWidth: 1080, minHeight: 700)
    }
}

// MARK: - AppKit Delegate & Native macOS Menu Bar

final class EpigenesisAppDelegate: NSObject, NSApplicationDelegate {
    let engine = EngineClient()
    var window: NSWindow!

    func applicationDidFinishLaunching(_ notification: Notification) {
        NSApp.setActivationPolicy(.regular)
        buildMainMenu()

        let rootView = MainWorkbenchView(engine: engine)
        let hosting = NSHostingController(rootView: rootView)

        window = NSWindow(
            contentRect: NSRect(x: 0, y: 0, width: 1320, height: 860),
            styleMask: [.titled, .closable, .miniaturizable, .resizable, .fullSizeContentView],
            backing: .buffered,
            defer: false
        )
        window.title = "Epigenesis"
        window.subtitle = "DNA-to-Development Workbench"
        window.toolbarStyle = .unified
        window.contentViewController = hosting
        window.center()
        window.makeKeyAndOrderFront(nil)
        NSApp.activate(ignoringOtherApps: true)
    }

    func application(_ application: NSApplication, open urls: [URL]) {
        if let first = urls.first {
            engine.openDocument(path: first.path)
        }
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { true }

    func applicationWillTerminate(_ notification: Notification) {
        engine.stopKernel()
    }

    private func buildMainMenu() {
        let mainMenu = NSMenu()

        let appItem = NSMenuItem()
        let appMenu = NSMenu()
        appMenu.addItem(withTitle: "About Epigenesis", action: #selector(showAbout(_:)), keyEquivalent: "")
        appMenu.addItem(.separator())
        appMenu.addItem(withTitle: "Hide Epigenesis", action: #selector(NSApplication.hide(_:)), keyEquivalent: "h")
        let hideOthers = appMenu.addItem(withTitle: "Hide Others", action: #selector(NSApplication.hideOtherApplications(_:)), keyEquivalent: "h")
        hideOthers.keyEquivalentModifierMask = [.command, .option]
        appMenu.addItem(withTitle: "Show All", action: #selector(NSApplication.unhideAllApplications(_:)), keyEquivalent: "")
        appMenu.addItem(.separator())
        appMenu.addItem(withTitle: "Quit Epigenesis", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
        appItem.submenu = appMenu
        mainMenu.addItem(appItem)

        let fileItem = NSMenuItem()
        let fileMenu = NSMenu(title: "File")
        fileMenu.addItem(withTitle: "Open DNA / Edit Table / Receipt…", action: #selector(openFile(_:)), keyEquivalent: "o")
        let openDir = fileMenu.addItem(withTitle: "Choose Workspace Folder…", action: #selector(openFolder(_:)), keyEquivalent: "O")
        openDir.keyEquivalentModifierMask = [.command, .shift]
        fileMenu.addItem(.separator())
        fileMenu.addItem(withTitle: "Save Cryptographic Receipt…", action: #selector(saveReceipt(_:)), keyEquivalent: "s")
        fileMenu.addItem(.separator())
        fileMenu.addItem(withTitle: "Close Window", action: #selector(NSWindow.performClose(_:)), keyEquivalent: "w")
        fileItem.submenu = fileMenu
        mainMenu.addItem(fileItem)

        let editItem = NSMenuItem()
        let editMenu = NSMenu(title: "Edit")
        editMenu.addItem(withTitle: "Undo", action: Selector(("undo:")), keyEquivalent: "z")
        editMenu.addItem(withTitle: "Redo", action: Selector(("redo:")), keyEquivalent: "Z")
        editMenu.addItem(.separator())
        editMenu.addItem(withTitle: "Cut", action: #selector(NSText.cut(_:)), keyEquivalent: "x")
        editMenu.addItem(withTitle: "Copy", action: #selector(NSText.copy(_:)), keyEquivalent: "c")
        editMenu.addItem(withTitle: "Paste", action: #selector(NSText.paste(_:)), keyEquivalent: "v")
        editMenu.addItem(withTitle: "Select All", action: #selector(NSText.selectAll(_:)), keyEquivalent: "a")
        editItem.submenu = editMenu
        mainMenu.addItem(editItem)

        let winItem = NSMenuItem()
        let winMenu = NSMenu(title: "Window")
        winMenu.addItem(withTitle: "Minimize", action: #selector(NSWindow.performMiniaturize(_:)), keyEquivalent: "m")
        winMenu.addItem(withTitle: "Zoom", action: #selector(NSWindow.performZoom(_:)), keyEquivalent: "")
        winItem.submenu = winMenu
        mainMenu.addItem(winItem)
        NSApp.windowsMenu = winMenu

        NSApp.mainMenu = mainMenu
    }

    @objc private func showAbout(_ sender: Any?) {
        let alert = NSAlert()
        alert.messageText = "Epigenesis"
        alert.informativeText = "DNA-to-Development Simulation, Inverse Regulatory Compiler, and Cryptographic Replay Workbench.\n\n100% Native macOS AppKit/SwiftUI · Direct Local Compute Engine (Zero Network Ports) · Apache-2.0."
        alert.alertStyle = .informational
        alert.runModal()
    }

    @objc private func openFile(_ sender: Any?) { engine.chooseFilePanel() }
    @objc private func openFolder(_ sender: Any?) { engine.chooseFolderPanel() }
    @objc private func saveReceipt(_ sender: Any?) { engine.exportCurrentReceipt() }
}

@main
struct EpigenesisLauncher {
    static func main() {
        let app = NSApplication.shared
        let delegate = EpigenesisAppDelegate()
        app.delegate = delegate
        app.run()
    }
}
