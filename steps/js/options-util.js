var ALGORITHMS = ['azimuth', 'celltypist', 'popv'];

function _find_algorithm(obj) {
  for (var spec of obj.algorithms) {
    for (var key of Object.keys(spec)) {
      if (ALGORITHMS.includes(key)) {
        return key;
      }
    }
  }
  return undefined;
}

function selectOutputDirectory(obj) {
  return obj['directory'] || _find_algorithm(obj) || '.';
}

function getSummarizeOptions(obj) {
  return {
    annotationMethod: _find_algorithm(obj) || 'unknown',
    ...(obj.algorithms.splice(-1)[0].summarize ?? {}),
  };
}
