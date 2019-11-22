 /* ODE solver using emscripten */
let params = []
var set1 = fs.readFileSync('./params.csv').toString()
var lines = set1.split('\n')
for (let i = 0; i < lines.length; i++) {
  params.push(lines[i].split(','))
}

let times = []
var set1 = fs.readFileSync('dataTimes.csv').toString()
var lines = set1.split('\n')
for (let i = 0; i < lines.length; i++) {
  times.push(lines[i].split(','))
}

let covarTime = []
var set1 = fs.readFileSync('./covarTimes.csv').toString()
var lines = set1.split('\n')
for (let i = 0; i < lines.length; i++) {
  covarTime.push(lines[i].split(','))
  covarTime[i] = Number(covarTime[i])
}

let covarTemperature = []
var set1 = fs.readFileSync('./covarTemperature').toString()
var lines = set1.split('\n')
for (let i = 0; i < lines.length; i++) {
  covarTemperature.push(lines[i].split(','))
}
console.log(covarTime[0])
function integrate (params, times, covarTime, covarTemperature) {
  let lsodaException = 0
  let arr = []
  let buffer
  let N = snippet.rInit(params)  
  let inputArray = Array(40).fill('number')
  let nByte = 8
  let lengthBuffer = times.length 
  lsodaTem = Module.cwrap('run_me', "number", inputArray)
  buffer = Module._malloc(lengthBuffer * nByte)

  /* Send covars' columns to C */
  let covarLength = covarTime.length;
  let covarTime_p = Module._malloc(covarLength * 8);
  let covarData_p = Module._malloc(covarLength * 8);
  for (let i = 0; i < covarLength; i++) {
      Module.setValue(covarTime_p + i * 8, covarTime[i], 'double');
      Module.setValue(covarData_p + i * 8, covarTemperature[i], 'double');
  }

  let timeAdd0 = [0].concat(times);
  let ptrTimes = Module._malloc(timeAdd0.length  * 8);
  for (let i = 0; i < timeAdd0.length; i++) {
      Module.setValue(ptrTimes + i * 8, timeAdd0[i], 'double');
  }
  lsodaException = lsodaTem(lengthBuffer,covarLength, buffer,ptrTimes, covarTime_p,covarData_p, ...N, ...params)
  if(lsodaException < 0){
    throw 'lsoda steps taken before reaching tout'
  } 
  for (var i = 0; i < lengthBuffer; i++) {
    arr.push(Module.getValue(buffer + i * nByte, 'double'))
  }
  
  Module._free(buffer)
  Module._free(covarTime_p)
  Module._free(covarData_p)
  Module._free(ptrTimes)
  return arr;
} 