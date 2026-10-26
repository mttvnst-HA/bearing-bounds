// Load the vendored Clarabel WebAssembly module in Node (the browser app inlines the same two files).
const fs = require("fs"), path = require("path");
module.exports = function loadClarabel() {
  const dir = path.join(__dirname, "..", "vendor", "clarabel");
  const wb = new Function(fs.readFileSync(path.join(dir, "clarabel_bb.js"), "utf8") + "\nreturn wasm_bindgen;")();
  wb.initSync({ module: fs.readFileSync(path.join(dir, "clarabel_bb_bg.wasm")) });
  return wb;
};
