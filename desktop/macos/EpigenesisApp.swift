import AppKit
import Foundation
import WebKit

/// Native macOS laboratory desktop wrapper for Epigenesis Studio.
/// Spawns the local loopback Python server, renders the interface inside a native
/// AppKit/WKWebView window with standard macOS menus, native file open/save dialogs,
/// and guarantees clean background-server shutdown on exit.
final class EpigenesisAppDelegate: NSObject, NSApplicationDelegate, NSWindowDelegate, WKNavigationDelegate, WKUIDelegate, WKScriptMessageHandler {
    private var window: NSWindow!
    private var webView: WKWebView!
    private var serverProcess: Process?
    private var serverURL: URL?
    private var pendingFileToOpen: String?
    private var statusField: NSTextField?

    func applicationDidFinishLaunching(_ notification: Notification) {
        NSApp.setActivationPolicy(.regular)
        buildMainMenu()
        buildMainWindow()
        startBackendAndLoad()
    }

    func application(_ application: NSApplication, open urls: [URL]) {
        guard let first = urls.first else { return }
        if serverURL == nil {
            pendingFileToOpen = first.path
        } else {
            openFileInWebView(url: first)
        }
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool {
        return true
    }

    func applicationWillTerminate(_ notification: Notification) {
        stopBackend()
    }

    // MARK: - Window & Menu Construction

    private func buildMainWindow() {
        let rect = NSRect(x: 0, y: 0, width: 1280, height: 860)
        window = NSWindow(
            contentRect: rect,
            styleMask: [.titled, .closable, .miniaturizable, .resizable],
            backing: .buffered,
            defer: false
        )
        window.title = "Epigenesis"
        window.minSize = NSSize(width: 1040, height: 680)
        window.center()
        window.delegate = self
        window.isReleasedWhenClosed = false

        let config = WKWebViewConfiguration()
        config.websiteDataStore = .nonPersistent()
        let contentController = WKUserContentController()
        contentController.add(self, name: "epigenesisSave")
        // Bridge blob/anchor downloads in Studio's JS to a native macOS NSSavePanel
        let downloadBridgeJS = """
        (function() {
          const origClick = HTMLAnchorElement.prototype.click;
          HTMLAnchorElement.prototype.click = function() {
            if (this.hasAttribute('download') && this.href && this.href.startsWith('blob:')) {
              const name = this.getAttribute('download') || 'epigenesis-output.json';
              fetch(this.href).then(r => r.text()).then(text => {
                window.webkit.messageHandlers.epigenesisSave.postMessage({ name: name, text: text });
              }).catch(() => {});
              return;
            }
            return origClick.apply(this, arguments);
          };
        })();
        """
        contentController.addUserScript(WKUserScript(source: downloadBridgeJS, injectionTime: .atDocumentStart, forMainFrameOnly: true))
        config.userContentController = contentController

        webView = WKWebView(frame: window.contentView!.bounds, configuration: config)
        webView.autoresizingMask = [.width, .height]
        webView.navigationDelegate = self
        webView.uiDelegate = self
        window.contentView?.addSubview(webView)

        let loadingHTML = """
        <!doctype html><html><head><meta charset="utf-8"><style>
        :root { color-scheme: light dark; }
        body { margin: 0; height: 100vh; display: grid; place-items: center; font: 14px/1.5 -apple-system, "SF Pro Text", sans-serif; background: #f2f2ef; color: #15171a; }
        @media (prefers-color-scheme: dark) { body { background: #0b0c0e; color: #ececea; } }
        .box { text-align: center; }
        .sub { opacity: 0.6; font-size: 13px; margin-top: 6px; }
        </style></head><body><div class="box"><strong>Epigenesis</strong><div class="sub">Starting local analysis engine…</div></div></body></html>
        """
        webView.loadHTMLString(loadingHTML, baseURL: nil)
        window.makeKeyAndOrderFront(nil)
        NSApp.activate(ignoringOtherApps: true)
    }

    private func buildMainMenu() {
        let mainMenu = NSMenu()

        // App menu
        let appItem = NSMenuItem()
        let appMenu = NSMenu()
        appMenu.addItem(withTitle: "About Epigenesis", action: #selector(showAbout(_:)), keyEquivalent: "")
        appMenu.addItem(NSMenuItem.separator())
        appMenu.addItem(withTitle: "Hide Epigenesis", action: #selector(NSApplication.hide(_:)), keyEquivalent: "h")
        let hideOthers = appMenu.addItem(withTitle: "Hide Others", action: #selector(NSApplication.hideOtherApplications(_:)), keyEquivalent: "h")
        hideOthers.keyEquivalentModifierMask = [.command, .option]
        appMenu.addItem(withTitle: "Show All", action: #selector(NSApplication.unhideAllApplications(_:)), keyEquivalent: "")
        appMenu.addItem(NSMenuItem.separator())
        appMenu.addItem(withTitle: "Quit Epigenesis", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
        appItem.submenu = appMenu
        mainMenu.addItem(appItem)

        // File menu
        let fileItem = NSMenuItem()
        let fileMenu = NSMenu(title: "File")
        fileMenu.addItem(withTitle: "Open DNA / Edit Table / Receipt…", action: #selector(openDocument(_:)), keyEquivalent: "o")
        fileMenu.addItem(NSMenuItem.separator())
        fileMenu.addItem(withTitle: "Reload", action: #selector(reloadPage(_:)), keyEquivalent: "r")
        let browserItem = fileMenu.addItem(withTitle: "Open in Default Browser", action: #selector(openInBrowser(_:)), keyEquivalent: "b")
        browserItem.keyEquivalentModifierMask = [.command, .shift]
        fileMenu.addItem(NSMenuItem.separator())
        fileMenu.addItem(withTitle: "Close Window", action: #selector(NSWindow.performClose(_:)), keyEquivalent: "w")
        fileItem.submenu = fileMenu
        mainMenu.addItem(fileItem)

        // Edit menu (enables native Cmd+C / Cmd+V / Cmd+A inside WKWebView textareas)
        let editItem = NSMenuItem()
        let editMenu = NSMenu(title: "Edit")
        editMenu.addItem(withTitle: "Undo", action: Selector(("undo:")), keyEquivalent: "z")
        editMenu.addItem(withTitle: "Redo", action: Selector(("redo:")), keyEquivalent: "Z")
        editMenu.addItem(NSMenuItem.separator())
        editMenu.addItem(withTitle: "Cut", action: #selector(NSText.cut(_:)), keyEquivalent: "x")
        editMenu.addItem(withTitle: "Copy", action: #selector(NSText.copy(_:)), keyEquivalent: "c")
        editMenu.addItem(withTitle: "Paste", action: #selector(NSText.paste(_:)), keyEquivalent: "v")
        editMenu.addItem(withTitle: "Select All", action: #selector(NSText.selectAll(_:)), keyEquivalent: "a")
        editItem.submenu = editMenu
        mainMenu.addItem(editItem)

        // View menu
        let viewItem = NSMenuItem()
        let viewMenu = NSMenu(title: "View")
        viewMenu.addItem(withTitle: "Actual Size", action: #selector(zoomReset(_:)), keyEquivalent: "0")
        viewMenu.addItem(withTitle: "Zoom In", action: #selector(zoomIn(_:)), keyEquivalent: "+")
        viewMenu.addItem(withTitle: "Zoom Out", action: #selector(zoomOut(_:)), keyEquivalent: "-")
        viewMenu.addItem(NSMenuItem.separator())
        let fullItem = viewMenu.addItem(withTitle: "Enter Full Screen", action: #selector(NSWindow.toggleFullScreen(_:)), keyEquivalent: "f")
        fullItem.keyEquivalentModifierMask = [.command, .control]
        viewItem.submenu = viewMenu
        mainMenu.addItem(viewItem)

        // Window menu
        let winItem = NSMenuItem()
        let winMenu = NSMenu(title: "Window")
        winMenu.addItem(withTitle: "Minimize", action: #selector(NSWindow.performMiniaturize(_:)), keyEquivalent: "m")
        winMenu.addItem(withTitle: "Zoom", action: #selector(NSWindow.performZoom(_:)), keyEquivalent: "")
        winItem.submenu = winMenu
        mainMenu.addItem(winItem)
        NSApp.windowsMenu = winMenu

        NSApp.mainMenu = mainMenu
    }

    // MARK: - Menu Actions

    @objc private func showAbout(_ sender: Any?) {
        let alert = NSAlert()
        alert.messageText = "Epigenesis Studio"
        alert.informativeText = "DNA-to-Development Simulation, Inverse DNA Compiler & Cryptographic Replay Workbench.\n\nRuns 100% locally on loopback (127.0.0.1). Apache-2.0."
        alert.alertStyle = .informational
        alert.runModal()
    }

    @objc private func openDocument(_ sender: Any?) {
        let panel = NSOpenPanel()
        panel.canChooseFiles = true
        panel.canChooseDirectories = false
        panel.allowsMultipleSelection = false
        panel.beginSheetModal(for: window) { response in
            guard response == .OK, let url = panel.url else { return }
            self.openFileInWebView(url: url)
        }
    }

    @objc private func reloadPage(_ sender: Any?) {
        if let url = serverURL {
            webView.load(URLRequest(url: url))
        }
    }

    @objc private func openInBrowser(_ sender: Any?) {
        if let url = serverURL {
            NSWorkspace.shared.open(url)
        }
    }

    @objc private func zoomReset(_ sender: Any?) { webView.pageZoom = 1.0 }
    @objc private func zoomIn(_ sender: Any?) { webView.pageZoom = min(2.0, webView.pageZoom + 0.1) }
    @objc private func zoomOut(_ sender: Any?) { webView.pageZoom = max(0.7, webView.pageZoom - 0.1) }

    // MARK: - Native Download & Open Bridges

    func userContentController(_ userContentController: WKUserContentController, didReceive message: WKScriptMessage) {
        guard message.name == "epigenesisSave",
              let body = message.body as? [String: Any],
              let name = body["name"] as? String,
              let text = body["text"] as? String else { return }
        DispatchQueue.main.async {
            let panel = NSSavePanel()
            panel.nameFieldStringValue = name
            panel.beginSheetModal(for: self.window) { resp in
                guard resp == .OK, let dest = panel.url else { return }
                do {
                    try text.write(to: dest, atomically: true, encoding: .utf8)
                } catch {
                    let errAlert = NSAlert(error: error)
                    errAlert.runModal()
                }
            }
        }
    }

    private func openFileInWebView(url: URL) {
        guard let text = try? String(contentsOf: url, encoding: .utf8) else { return }
        let nameJSON = (try? JSONSerialization.data(withJSONObject: [url.lastPathComponent]))
            .flatMap { String(data: $0, encoding: .utf8) } ?? "[\"file\"]"
        let textJSON = (try? JSONSerialization.data(withJSONObject: [text]))
            .flatMap { String(data: $0, encoding: .utf8) } ?? "[\"\"]"
        let js = """
        (async function() {
          const name = (\(nameJSON))[0], text = (\(textJSON))[0];
          if (typeof api === 'function' && typeof openDoc === 'function') {
            try {
              const res = await api('/api/kind', { name, text });
              openDoc({ name, text, kind: res.kind });
            } catch (e) {}
          }
        })();
        """
        webView.evaluateJavaScript(js, completionHandler: nil)
    }

    // MARK: - Backend Python Server Lifecycle

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
            if let r = obj["repo"], fm.fileExists(atPath: "\(r)/studio/epigenesis_studio/server.py") { repo = r }
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
            repo = repoCandidates.first(where: { !$0.isEmpty && fm.fileExists(atPath: "\($0)/studio/epigenesis_studio/server.py") }) ?? fm.currentDirectoryPath
        }
        return (py, repo)
    }

    private func startBackendAndLoad() {
        let cfg = loadRuntimeConfig()
        let proc = Process()
        proc.executableURL = URL(fileURLWithPath: cfg.python)
        var args = ["-m", "epigenesis_studio.server", "--port", "0", "--no-browser"]
        if let f = pendingFileToOpen {
            args.append(f)
        }
        proc.arguments = args
        proc.currentDirectoryURL = URL(fileURLWithPath: cfg.repo)

        var env = ProcessInfo.processInfo.environment
        let pyPaths = [cfg.repo, "\(cfg.repo)/studio"]
        if let existing = env["PYTHONPATH"], !existing.isEmpty {
            env["PYTHONPATH"] = (pyPaths + [existing]).joined(separator: ":")
        } else {
            env["PYTHONPATH"] = pyPaths.joined(separator: ":")
        }
        env["PYTHONUNBUFFERED"] = "1"
        proc.environment = env

        let outPipe = Pipe()
        let errPipe = Pipe()
        proc.standardOutput = outPipe
        proc.standardError = errPipe

        do {
            try proc.run()
        } catch {
            showFatalError("Could not launch Epigenesis backend (\(cfg.python)): \(error.localizedDescription)")
            return
        }
        serverProcess = proc

        // Read stdout in background until the server announces its loopback URL
        DispatchQueue.global(qos: .userInitiated).async {
            let handle = outPipe.fileHandleForReading
            var buffer = ""
            while proc.isRunning {
                let chunk = handle.availableData
                if chunk.isEmpty { break }
                if let str = String(data: chunk, encoding: .utf8) {
                    buffer += str
                    if let range = buffer.range(of: #"http://127\.0\.0\.1:\d+/"#, options: .regularExpression),
                       let url = URL(string: String(buffer[range])) {
                        DispatchQueue.main.async {
                            self.serverURL = url
                            self.webView.load(URLRequest(url: url))
                        }
                        return
                    }
                }
            }
            let errText = String(data: errPipe.fileHandleForReading.readDataToEndOfFile(), encoding: .utf8) ?? ""
            DispatchQueue.main.async {
                if self.serverURL == nil {
                    self.showFatalError("Epigenesis backend exited before starting.\n\n\(errText)")
                }
            }
        }
    }

    private func showFatalError(_ message: String) {
        let escaped = message
            .replacingOccurrences(of: "&", with: "&amp;")
            .replacingOccurrences(of: "<", with: "&lt;")
            .replacingOccurrences(of: ">", with: "&gt;")
        let html = """
        <!doctype html><html><head><meta charset="utf-8"><style>
        body { font: 14px/1.5 -apple-system, sans-serif; padding: 36px; color: #15171a; background: #f2f2ef; }
        pre { background: #fff; border: 1px solid #d4d4d4; padding: 14px; border-radius: 8px; white-space: pre-wrap; }
        </style></head><body><h2>Epigenesis could not start</h2><pre>\(escaped)</pre></body></html>
        """
        webView.loadHTMLString(html, baseURL: nil)
    }

    private func stopBackend() {
        if let proc = serverProcess, proc.isRunning {
            proc.terminate()
            proc.waitUntilExit()
        }
        serverProcess = nil
    }
}

let app = NSApplication.shared
let delegate = EpigenesisAppDelegate()
app.delegate = delegate
app.run()
