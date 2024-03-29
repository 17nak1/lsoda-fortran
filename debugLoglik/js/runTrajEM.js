/**
 *  @file       runTraj.js        
 *              This function attempts to match trajectories of a model's deterministic skeleton to data.
 *              Trajectory matching is equivalent to maximum likelihood estimatedation under the assumption 
 *              that process noise is entirely absent, i.e., that all stochasticity is measurement error.
 *              Accordingly, this method uses only the skeleton and dmeasure components of a POMP model.
 *
 *  @author     Nazila Akhavan, nazila@kingsds.network
 *  @date       July 2019
 */

let snippet = require('./modelSnippet.js')
let model = require('./createModel')
let create  = require('./create.js')
let mathLib = require('./mathLib')
let Index = require('./indices')
let DetermineRun = require('./determineRun')
let fmin    = require('fmin')
let fs = require('fs')
const Module = require('./lsoda.js')
/************************************************************ Will be defined in ui ************************************************/
let dt = 0.005 // Step size only use in covar
let startTime = 1991
let endTime = 2008
let estIcstart = [0] // Are there no initial conditions given? 0-Given, 1-No, 2-TrajMatch
let run = 1; 
let t0 = 0;


// Parameters that is consider always fixed
let paramsIcFixed = []//snippet.statenames()
let paramsFixed =[Index.p, Index.delta, Index.mu_e, Index.mu_ql, Index.mu_el, Index.mu_qn, Index.mu_en, Index.mu_qa, Index.mu_ea,
          Index.mu_h, Index.beta_nh, Index.beta_hn, Index.beta_hl, Index.alpha, Index.c, Index.Tf, Index.gamma, ...paramsIcFixed]

// Parameters not to be transformed
let paramsNotrans = [].concat(paramsFixed)              
let ParamSetFile, paramProf
if (run === 1) {
  ParamSetFile = "../data/ParamSet_TBE.csv" 
  paramProf = null 
} else {
  ParamSetFile = `ParamSet_run${run}.csv`    
  paramProf = DetermineRun.type(run).paramProf
}  

paramsFixed = [...paramsFixed, paramProf]
paramsNoic = [Index.p, Index.omega, Index.delta, Index.mu_e, Index.mu_ql, Index.mu_el, Index.mu_qn, Index.mu_en, Index.mu_qa,
              Index.mu_ea, Index.mu_h, Index.beta_nh, Index.beta_hl,Index.beta_hn, Index.lambda_l, Index.lambda_n, Index.lambda_a,
              Index.alpha, Index.f_l, Index.f_n, Index.f_a,  Index.kappa, Index.c,Index.Tf, Index.obsprob, Index.T_min_l,Index.gamma]
paramsIc  = snippet.statenames()

// paramsFit =  paramsNoic - paramfixed
let paramsFit = []
let flag = 0
for (let i = 0; i < paramsNoic.length; i++) {
  for ( let j = 0; j < paramsFixed.length; j++) {
    if ( paramsNoic[i] === paramsFixed[j]) {
      flag = 1
      break
    }
  }
  if(flag === 0) {
    paramsFit.push(paramsNoic[i])
  }
flag = 0
}

// paramsIcFit = paramsIc - paramfixed
let paramsIcFit = []
flag = 0
for (let i = 0; i < paramsIc.length; i++) {
  for ( let j = 0; j < paramsFixed.length; j++) {
    if ( paramsIc[i] === paramsFixed[j]) {
      flag = 1
      break
    }
  }
  if(flag === 0) {
    paramsIcFit.push(paramsIc[i])
  }
flag = 0
}


let index = Array(40).fill(0)
let estimatedIndex = [...paramsFit, ...paramsIcFit]
for ( let i = 0; i < estimatedIndex.length; i++) {
  index[estimatedIndex[i]] = 1
}
let place = estimatedIndex
let fullset = []
var set1 = fs.readFileSync(ParamSetFile).toString()
var lines = set1.split('\n')
for (let i = 0; i < lines.length; i++) {
  fullset.push(lines[i].split(','))
}

// Generate covars and data  
<<<<<<< HEAD
let covars = create.covars(startTime, endTime, dt);
=======
let covars = create.covars(startTime, endTime, dt)
>>>>>>> 18ff24eef4e4e14d6b6c36a601beb3fd644922b2
let covarTime = [], covarTemperature = []
for (let i = 0; i < covars.length; i++) {
  covarTime.push(covars[i][0])
  covarTemperature.push(covars[i][1])
}
let data = create.dataset(startTime, endTime)
let times = []
for (let i = 0; i < data.length; i++) {
  times.push(data[i][0])
}
/**************************************************************************************************************************************************/
function traj_match (data, covarTime, covarTemperature, params, times, t0, index, place) {
  let deltaT = (1 / 52) * 365
  var estimated = []
  // Index of parameters that need to be transfered
  let temp = model.createPompModel(data, covars, t0 = 0, dt = 0.005, paramsFixed)
  let logTrans = temp[0]
  let logitTrans = temp[1]
  
  // Change the parameters' scale 
  model.toEstimationScale(params, logTrans, logitTrans)
 
  // Choose those that should be estimated.
  for (let i = 0; i < index.length; i++) {
    if (index[i] === 1 ) {
      estimated.push(params[i])
    }
  }

  let lo = logLik (estimated)
  
  //* calculate log likelihood
  function logLik (estimated) {
    var likvalue = 0
    var loglik = 0
    var rho 
    var psi
    var simHarranged = []
    for (let i = 0; i < estimated.length; i++) {
      params[place[i]] = estimated[i]
    }

    // Return parameters' scale to original
    model.fromEstimationScale(params, logTrans, logitTrans)
    // console.log(covars[0])
    var simH = integrate(params, times, deltaT,covarTime, covarTemperature)
    simHarranged[0] = simH[t0]
    var aa = [[0,0]]
    for ( let i = 1; i < simH.length; i++) {
      simHarranged[i] = simH[i] - simH[i - 1]
      aa.push([i , simHarranged[i]])
    }
    
    for (let i = 0; i < simHarranged.length; i++) {
      likvalue = snippet.dObs(params[Index.obsprob], simHarranged[i], data[i][1], 1)//;ar.push([likvalue])
      loglik = loglik + likvalue
    }
    console.log(params, loglik)
    return [-(loglik).toFixed(6)]
  }
  return[...params, -lo]
}

 /* ODE solver using emscripten */
function integrate (params, times, deltaT, covarTime, covarTemperature) {
  let lsodaException = 0
  let arr = []
  let buffer
  let N = snippet.rInit(params)  
  let inputArray = Array(40).fill('number')
  let nByte = 8
<<<<<<< HEAD
  let lengthBuffer = times.length 
=======
  let lengthBuffer = times.length  

>>>>>>> 18ff24eef4e4e14d6b6c36a601beb3fd644922b2
  lsodaTem = Module.cwrap('run_me', "number", inputArray)
  buffer = Module._malloc(lengthBuffer * nByte)

  /* Send covars' columns to C */
  let covarLength = covarTime.length;
  let covarTime_p = Module._malloc(covarLength * 8);
  let covarData_p = Module._malloc(covarLength * 8);
  for (let i = 0; i < covarLength; i++) {
      Module.setValue(covarTime_p + i * 8, covarTime[i], 'double');
      Module.setValue(covarData_p + i * 8, covarTemperature[i], 'double');
<<<<<<< HEAD
  }

  let timeAdd0 = [0].concat(times);
  let ptrTimes = Module._malloc(timeAdd0.length  * 8);
  for (let i = 0; i < timeAdd0.length; i++) {
      Module.setValue(ptrTimes + i * 8, timeAdd0[i], 'double');
  }
=======
  }

  let timeAdd0 = [0].concat(times);
  let ptrTimes = Module._malloc(timeAdd0.length  * 8);
  for (let i = 0; i < timeAdd0.length; i++) {
      Module.setValue(ptrTimes + i * 8, timeAdd0[i], 'double');
  }

>>>>>>> 18ff24eef4e4e14d6b6c36a601beb3fd644922b2
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
<<<<<<< HEAD
} 
=======
}
 
>>>>>>> 18ff24eef4e4e14d6b6c36a601beb3fd644922b2


/** Main program entry point */
function main() {
  let resultSet = [], result

<<<<<<< HEAD
  for(let count = 1; count <= 300; count++) {
=======
  for(let count = 1; count <= 2; count++) {
>>>>>>> 18ff24eef4e4e14d6b6c36a601beb3fd644922b2
    var params = []
    for ( let i = 0; i < fullset[0].length; i++) {
      params.push(Number(fullset[count][i]))
    }
    try{
      result = traj_match (data, covarTime, covarTemperature, params, times, t0, index, place);   
      resultSet.push(result);
    }
    catch(e) {
      resultSet.push([...Array(params.length).fill(0), NaN]);
      console.error(e);
    }
  }
  const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['p' , 'omega' , 'delta' , 'mu_e' , 'mu_ql' , 'mu_el' , 'mu_qn' , 'mu_en' , 'mu_qa' , 'mu_ea' , 'mu_h' ,
     'beta_nh' , 'beta_hl' , 'beta_hn' , 'lambda_l' , 'lambda_n' , 'lambda_a' , 'alpha' , 'f_l' , 'f_n' , 'f_a' , 'kappa' , 
     'c' , 'Tf' , 'obsprob' , 'T_min_l' , 'gamma' , 'E0' , 'QL0' , 'EL_s0' , 'EL_i0' , 'QN_s0' , 'QN_i0' , 'EN_s0' , 'EN_i0' ,
      'QA_s0' , 'QA_i0' , 'EA0' , 'H_s0' , 'H_i0' , 'LogLik'],
<<<<<<< HEAD
    path: './oct15.csv'
=======
    path: './test.csv'
>>>>>>> 18ff24eef4e4e14d6b6c36a601beb3fd644922b2
  })   
  csvWriter.writeRecords(resultSet)
    .then(() => {
    console.log('...Done')
  })
}

/* Run main only when emscripten is ready */
Module.onRuntimeInitialized = main

/*emcc lsoda.c -o lsoda.js -s  EXPORTED_FUNCTIONS='["_run_me"]' -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap","getValue", "setValue"]' -s EXIT_RUNTIME=1
*/