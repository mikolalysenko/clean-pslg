'use strict'

module.exports = cleanPSLG

var UnionFind = require('union-find')
var boxIntersect = require('box-intersect')
var compareCell = require('compare-cell')
var segseg = require('robust-segment-intersect')
var rat = require('big-rat')
var ratCmp = require('big-rat/cmp')
var ratToFloat = require('big-rat/to-float')
var ratVec = require('rat-vec')
var nextafter = require('nextafter')

var solveIntersection = require('./lib/rat-seg-intersect')

//Bounds on a rational number when rounded to a float
function boundRat(r) {
  var f = ratToFloat(r)
  var cmp = ratCmp(rat(f), r)
  if(cmp < 0) {
    return [f, nextafter(f, Infinity)]
  } else if(cmp > 0) {
    return [nextafter(f, -Infinity), f]
  } else {
    return [f, f]
  }
}

//Convert a list of edges in a pslg to bounding boxes
function boundEdges(points, edges) {
  var bounds = new Array(edges.length)
  for(var i=0; i<edges.length; ++i) {
    var e = edges[i]
    var a = points[e[0]]
    var b = points[e[1]]
    bounds[i] = [
      Math.min(a[0], b[0]),
      Math.min(a[1], b[1]),
      Math.max(a[0], b[0]),
      Math.max(a[1], b[1]) ]
  }
  return bounds
}

//Convert a list of points into bounding boxes by duplicating coords
function boundPoints(points) {
  var bounds = new Array(points.length)
  for(var i=0; i<points.length; ++i) {
    var p = points[i]
    bounds[i] = [ p[0], p[1], p[0], p[1] ]
  }
  return bounds
}

//Find all pairs of crossing edges in a pslg (given edge bounds)
function getCrossings(points, edges, edgeBounds) {
  var result = []
  boxIntersect(edgeBounds, function(i, j) {
    var e = edges[i]
    var f = edges[j]
    if(e[0] === f[0] || e[0] === f[1] ||
       e[1] === f[0] || e[1] === f[1]) {
         return
    }
    var a = points[e[0]]
    var b = points[e[1]]
    var c = points[f[0]]
    var d = points[f[1]]
    if(segseg(a, b, c, d)) {
      result.push([i, j])
    }
  })
  return result
}

//Find all pairs of crossing vertices in a pslg (given edge/vert bounds)
function getTJunctions(points, edges, edgeBounds, vertBounds) {
  var result = []
  boxIntersect(edgeBounds, vertBounds, function(i, v) {
    var e = edges[i]
    if(e[0] === v || e[1] === v) {
      return
    }
    var p = points[v]
    var a = points[e[0]]
    var b = points[e[1]]
    if(segseg(a, b, p, p)) {
      result.push([i, v])
    }
  })
  return result
}


//Cut edges along crossings/tjunctions
function cutEdges(floatPoints, edges, crossings, junctions, useColor) {

  //Convert crossings into tjunctions by constructing rational points
  var ratPoints = []
  for(var i=0; i<crossings.length; ++i) {
    var crossing = crossings[i]
    var e = crossing[0]
    var f = crossing[1]
    var ee = edges[e]
    var ef = edges[f]
    var x = solveIntersection(
      ratVec(floatPoints[ee[0]]),
      ratVec(floatPoints[ee[1]]),
      ratVec(floatPoints[ef[0]]),
      ratVec(floatPoints[ef[1]]))
    if(!x) {
      //Segments are parallel, should already be handled by t-junctions
      continue
    }
    var idx = ratPoints.length + floatPoints.length
    ratPoints.push(x)
    junctions.push([e, idx], [f, idx])
  }

  //Sort tjunctions
  function getPoint(idx) {
    if(idx >= floatPoints.length) {
      return ratPoints[idx-floatPoints.length]
    }
    var p = floatPoints[idx]
    return [ rat(p[0]), rat(p[1]) ]
  }
  junctions.sort(function(a, b) {
    if(a[0] !== b[0]) {
      return a[0] - b[0]
    }
    var u = getPoint(a[1])
    var v = getPoint(b[1])
    return ratCmp(u[0], v[0]) || ratCmp(u[1], v[1])
  })

  //Split edges along junctions
  for(var i=junctions.length-1; i>=0; --i) {
    var junction = junctions[i]
    var e = junction[0]

    var edge = edges[e]
    var s = edge[0]
    var t = edge[1]

    //Check if edge is not lexicographically sorted
    var a = floatPoints[s]
    var b = floatPoints[t]
    if(((a[0] - b[0]) || (a[1] - b[1])) < 0) {
      var tmp = s
      s = t
      t = tmp
    }

    //Split leading edge
    edge[0] = s
    var last = edge[1] = junction[1]

    //If we are grouping edges by color, remember to track data
    var color
    if(useColor) {
      color = edge[2]
    }

    //Split other edges
    while(i > 0 && junctions[i-1][0] === e) {
      var junction = junctions[--i]
      var next = junction[1]
      if(useColor) {
        edges.push([last, next, color])
      } else {
        edges.push([last, next])
      }
      last = next
    }

    //Add final edge
    if(useColor) {
      edges.push([last, t, color])
    } else {
      edges.push([last, t])
    }
  }

  //Return constructed rational points
  return ratPoints
}

//Merge overlapping points
function dedupPoints(floatPoints, ratPoints, floatBounds) {
  var numPoints = floatPoints.length + ratPoints.length
  var uf        = new UnionFind(numPoints)

  //Compute rational bounds
  var bounds = floatBounds
  for(var i=0; i<ratPoints.length; ++i) {
    var p = ratPoints[i]
    var xb = boundRat(p[0])
    var yb = boundRat(p[1])
    bounds.push([ xb[0], yb[0], xb[1], yb[1] ])
    floatPoints.push([ ratToFloat(p[0]), ratToFloat(p[1]) ])
  }

  //Link all points with over lapping boxes
  boxIntersect(bounds, function(i, j) {
    uf.link(i, j)
  })

  //Call find on each point to get a relabeling
  var ptr = 0
  var noDupes = true
  var labels = new Array(numPoints)
  for(var i=0; i<numPoints; ++i) {
    var j = uf.find(i)
    if(j === i) {
      //If not a duplicate, then don't bother
      labels[i] = ptr
      floatPoints[ptr++] = floatPoints[i]
    } else {
      //Clear no-dupes flag, zero out label
      noDupes = false
      labels[i] = -1
    }
  }
  floatPoints.length = ptr

  //If no duplicates, return null to signal termination
  if(noDupes) {
    return null
  }

  //Do a second pass to fix up missing labels
  for(var i=0; i<numPoints; ++i) {
    if(labels[i] < 0) {
      labels[i] = labels[uf.find(i)]
    }
  }

  //Return resulting union-find data structure
  return labels
}

function compareLex2(a,b) { return (a[0]-b[0]) || (a[1]-b[1]) }
function compareLex3(a,b) {
  var d = (a[0] - b[0]) || (a[1] - b[1])
  if(d) {
    return d
  }
  if(a[2] < b[2]) {
    return -1
  } else if(a[2] > b[2]) {
    return 1
  }
  return 0
}

//Remove duplicate edge labels
function dedupEdges(edges, labels, useColor) {
  if(edges.length === 0) {
    return
  }
  if(labels) {
    for(var i=0; i<edges.length; ++i) {
      var e = edges[i]
      var a = labels[e[0]]
      var b = labels[e[1]]
      e[0] = Math.min(a, b)
      e[1] = Math.max(a, b)
    }
  } else {
    for(var i=0; i<edges.length; ++i) {
      var e = edges[i]
      var a = e[0]
      var b = e[1]
      e[0] = Math.min(a, b)
      e[1] = Math.max(a, b)
    }
  }
  if(useColor) {
    edges.sort(compareLex3)
  } else {
    edges.sort(compareLex2)
  }
  var ptr = 1
  for(var i=1; i<edges.length; ++i) {
    var prev = edges[i-1]
    var next = edges[i]
    if(next[0] === prev[0] && next[1] === prev[1] &&
      (!useColor || next[2] === prev[2])) {
      continue
    }
    edges[ptr++] = next
  }
  edges.length = ptr
}

//Repeat until convergence
function snapRound(points, edges, useColor) {

  // 1. find edge crossings
  var edgeBounds = boundEdges(points, edges)
  var crossings  = getCrossings(points, edges, edgeBounds)

  // 2. find t-junctions
  var vertBounds = boundPoints(points)
  var tjunctions = getTJunctions(points, edges, edgeBounds, vertBounds)

  // 3. cut edges, construct rational points
  var ratPoints  = cutEdges(points, edges, crossings, tjunctions, useColor)

  // 4. dedupe verts
  var labels     = dedupPoints(points, ratPoints, vertBounds)

  // 6. dedupe edges
  dedupEdges(edges, labels, useColor)

  // 5. check termination
  if(!labels) {
    return (crossings.length > 0 || tjunctions.length > 0)
  }

  // More iterations necessary
  return true
}

//Main loop, runs PSLG clean up until completion
function cleanPSLG(points, edges, colors) {
  var modified = false

  //If using colors, augment edges with color data
  var prevEdges
  if(colors) {
    prevEdges = edges
    var augEdges = new Array(edges.length)
    for(var i=0; i<edges.length; ++i) {
      var e = edges[i]
      augEdges[i] = [e[0], e[1], colors[i]]
    }
    edges = augEdges
  }

  //Run snap rounding until convergence
  while(snapRound(points, edges, !!colors)) {
    modified = true
  }

  //Strip color tags
  if(!!colors && modified) {
    prevEdges.length = 0
    colors.length = 0
    for(var i=0; i<edges.length; ++i) {
      var e = edges[i]
      prevEdges.push([e[0], e[1]])
      colors.push(e[2])
    }
  }

  return modified
}
