var ALGORITHMS = ["azimuth", "celltypist", "popv"];

function _find_algorithm(obj) {
  for (var index = 0; index < ALGORITHMS.length; ++index) {
    var name = ALGORITHMS[index];
    if (typeof obj[name] === "object" && obj[name] !== null) {
      return name;
    }
  }

  return null;
}

function selectOutputDirectory(obj) {
  return obj["directory"] || _find_algorithm(obj) || ".";
}

function getDefaultSummarizeOptions(obj) {
  return {
    annotationMethod: _find_algorithm(obj) || "unknown",
  };
}
