var ALGORITHMS = ["azimuth", "celltypist", "popv"];

/**
 * Finds the algorithm selected in the option object
 *
 * @param {object} obj Options
 * @returns {string} Name of algorithm or null if no match was found
 */
function _find_algorithm(obj) {
  for (var index = 0; index < ALGORITHMS.length; ++index) {
    var name = ALGORITHMS[index];
    if (typeof obj[name] === "object") {
      return name;
    }
  }

  return null;
}

/**
 * Selects an output directory based on the provided options
 *
 * @param {object} obj Options
 * @returns The selected output directory
 */
function selectOutputDirectory(obj) {
  return obj["directory"] || _find_algorithm(obj) || ".";
}

/**
 * Creates default summarize step options based on the provided options
 *
 * @param {object} obj Options
 * @returns Default summarize options
 */
function getDefaultSummarizeOptions(obj) {
  return {
    annotationMethod: _find_algorithm(obj) || "unknown",
  };
}
