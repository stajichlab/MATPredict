// Applies the saved theme before first paint and wires the #theme selector. "auto" follows the system.
(function () {
  var KEY = "matpredict-theme", THEMES = ["auto", "light", "dark", "paper", "neuromancer"];
  function get() { try { var v = localStorage.getItem(KEY); return THEMES.indexOf(v) >= 0 ? v : "auto"; } catch (e) { return "auto"; } }
  function apply(v) { if (v === "auto") document.documentElement.removeAttribute("data-theme"); else document.documentElement.setAttribute("data-theme", v); }
  apply(get());
  document.addEventListener("DOMContentLoaded", function () {
    var sel = document.getElementById("theme"); if (!sel) return;
    sel.value = get();
    sel.addEventListener("change", function () { apply(sel.value); try { localStorage.setItem(KEY, sel.value); } catch (e) {} });
  });
})();
