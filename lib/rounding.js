var ONE_PLUS_ULP  = 1.0 + Math.pow(2, -52)
var ONE_MINUS_ULP = 1.0 - Math.pow(2, -53)
var DENORMAL      = Math.pow(2, -1021)
var MIN_FLOAT     = Number.MIN_VALUE
var MAX_FLOAT     = Number.MAX_VALUE

function roundUp(x) {
  x = +x
  if(x > +0) {
    if(x < DENORMAL) {
      return +(x + MIN_FLOAT)
    }
    return +(x * ONE_PLUS_ULP)
  } else if(x < 0) {
    if(x > -DENORMAL) {
      return +(x + MIN_FLOAT)
    } else if(x === -Infinity) {
      return -MAX_FLOAT
    }
    return +(x * ONE_MINUS_ULP)
  } else {
    return +MIN_FLOAT
  }
}
exports.up = roundUp

function roundDown(x) {
  x = +x
  if(x > 0) {
    if(x < DENORMAL) {
      return +(x - MIN_FLOAT)
    } else if(x === Infinity) {
      return +MAX_FLOAT
    }
    return +(x * ONE_MINUS_ULP)
  } else if(x < 0) {
    if(x > -DENORMAL) {
      return +(x - MIN_FLOAT)
    }
    return +(x * ONE_PLUS_ULP)
  } else {
    return -MIN_FLOAT
  }
}
exports.down = roundDown
