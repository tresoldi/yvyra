/*
 * yvyra Web Worker
 * Pre-initializes the WASM module on load, then runs analyses on demand.
 */

var Module = null;
var ready = false;

try {
    importScripts("yvyra.js?v=" + Date.now());
} catch(e) {
    self.postMessage({ type: "error", message: "Failed to load yvyra.js: " + e.message });
}

/* Pre-initialize WASM module immediately */
if (typeof createYvyra === "function") {
    self.postMessage({ type: "status", text: "Loading yvyra engine (1.3 MB)..." });

    createYvyra({
        noInitialRun: true,
        locateFile: function(path) { return path + "?v=" + Date.now(); },
        print: function (text) {
            self.postMessage({ type: "output", text: text });
            var m = text.match(/^\s+(\d+)\s+--/);
            if (m) {
                self.postMessage({ type: "progress", generation: parseInt(m[1]) });
            }
        },
        printErr: function (text) {
            self.postMessage({ type: "output", text: text });
        }
    }).then(function(mod) {
        Module = mod;
        ready = true;
        self.postMessage({ type: "ready" });
    }).catch(function(err) {
        self.postMessage({ type: "error", message: "WASM init failed: " + (err.message || String(err)) });
    });
} else {
    self.postMessage({ type: "error", message: "createYvyra not found after importScripts" });
}

self.onmessage = function (e) {
    if (e.data.type === "run") {
        if (!ready || !Module) {
            self.postMessage({ type: "error", message: "Engine not ready" });
            return;
        }

        var content = e.data.yaml || e.data.nexus;
        var filename = e.data.yaml ? "input.yaml" : "input.nex";

        try {
            /* Clean virtual filesystem */
            try {
                var files = Module.FS.readdir(".");
                for (var i = 0; i < files.length; i++) {
                    var f = files[i];
                    if (f === "." || f === ".." || f === "tmp" || f === "home" || f === "dev" || f === "proc") continue;
                    try { Module.FS.unlink(f); } catch(ex) {}
                }
            } catch(ex) {}

            Module.FS.writeFile(filename, content);
            self.postMessage({ type: "status", text: "Running analysis..." });

            var exitCode = Module.callMain([filename]);

            /* Send generated NEXUS back for display */
            if (filename.match(/\.yaml$/)) {
                try {
                    var nexName = filename.replace(/\.yaml$/, ".nex");
                    var nexContent = Module.FS.readFile(nexName, { encoding: "utf8" });
                    self.postMessage({ type: "nexus", text: nexContent });
                } catch(ex) { /* no .nex generated */ }
            }

            /* Collect output files */
            var outFiles = {};
            var fsRoot = Module.FS.readdir(".");
            for (var i = 0; i < fsRoot.length; i++) {
                var f = fsRoot[i];
                if (f === "." || f === ".." || f === filename) continue;
                if (f.match(/\.(p|t|slk|con\.tre|parts|tstat|pstat|trprobs|mcmc|vstat|lstat)$/)) {
                    try { outFiles[f] = Module.FS.readFile(f, { encoding: "utf8" }); } catch(ex) {}
                }
            }

            self.postMessage({ type: "done", exitCode: exitCode, files: outFiles });
        } catch (err) {
            self.postMessage({ type: "error", message: err.message || String(err) });
        }
    }
};
