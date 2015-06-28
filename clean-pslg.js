'use strict'

module.exports = cleanPSLG

var UnionFind = require('union-find')
var boxIntersect = require('box-intersect')
var uniq = require('uniq')
var compareCell = require('compare-cell')
var segseg = require('robust-segment-intersect')
var rat = require('big-rat')
var ratCmp = require('big-rat/cmp')
var ratToFloat = require('big-rat/to-float')
var ratVec = require('rat-vec')

var round = require('./lib/rounding')
var solveIntersection = require('./lib/rat-seg-intersect')

//Bounds on a rational number when rounded to a float
//TODO: Use nextafter here
function boundRat(r) {
  var f = ratToFloat(r)
  var cmp = ratCmp(rat(f), r)
  if(cmp < 0) {
    return [f, round.up(f)]
  } else if(cmp > 0) {
    return [round.down(f), f]
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
function cutEdges(floatPoints, edges, crossings, junctions) {

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

    //Split other edges
    while(i > 0 && junctions[i-1][0] === e) {
      var junction = junctions[--i]
      var next = junction[1]
      edges.push([last, next])
      last = next
    }

    //Add final edge
    edges.push([last, t])
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
  for(var i=0; i<numPoints; ++i) {
    var j = uf.find(i)
    if(j === i) {
      //If not a duplicate, then don't bother
      floatPoints[ptr++] = floatPoints[i]
    } else {
      //Clear no-dupes flag
      noDupes = false

      //If order is not consistent, then swap in tree
      if(i < j) {
        var rankI = uf.ranks[i]
        var rankJ = uf.ranks[j]
        uf.ranks[i] = rankJ
        uf.ranks[j] = rankI
        uf.roots[j] = uf.roots[i] = i
        floatPoints[ptr++] = floatPoints[i]
      }
    }
  }
  floatPoints.length = ptr

  //If no duplicates, return null to signal termination
  if(noDupes) {
    return null
  }

  //Return resulting union-find data structure
  return uf.roots
}

//Remove duplicate edge labels
function dedupEdges(edges, labels) {
  for(var i=0; i<edges.length; ++i) {
    var e = edges[i]
    e[0] = labels[e[0]]
    e[1] = labels[e[1]]
  }
  uniq(edges, compareCell)
}

//Repeat until convergence
function snapRound(points, edges) {

  // 1. find edge crossings
  var edgeBounds = boundEdges(points, edges)
  var crossings  = getCrossings(points, edges, edgeBounds)

  // 2. find t-junctions
  var vertBounds = boundPoints(points)
  var tjunctions = getTJunctions(points, edges, edgeBounds, vertBounds)

  // 3. cut edges, construct rational points
  var ratPoints  = cutEdges(points, edges, crossings, tjunctions)

  // 4. dedupe verts
  var labels     = dedupPoints(points, ratPoints, vertBounds)

  // 5. check termination
  if(!labels) {
    return crossings.length > 0 || tjunctions.length > 0
  }

  // 6. dedupe edges
  dedupEdges(edges, labels)

  // More iterations necessary
  return true
}

//Main loop, runs PSLG clean up until completion
function cleanPSLG(points, edges) {
  var cleaned = false
  while(snapRound(points, edges)) {
    cleaned = true
  }
  return cleaned
}
