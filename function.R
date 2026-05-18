#Function definition
.systemInfo = function(eachThreadDir, spid, type, idx) {
  tempInfo = sprintf('%s%s.%s.', eachThreadDir, spid, type)
  return(paste0(tempInfo, idx, '.trds'))
}

function.remove = function(paths, check = FALSE) {
  stopifnot(is.character(paths))
  stopifnot(is.logical(check))
  
  if (check) {
    existsVector = file.exists(paths) | dir.exists(paths)
    invisible(suppressWarnings(file.remove(paths[existsVector])))
  } else {
    invisible(suppressWarnings(file.remove(paths)))
  }
}

function.doGC = function() {
  # silence garbage collection (gc)
  invisible(gc())
}

function.getPID = function() {
  # Get the process
  NODE = unlist(strsplit(Sys.info()['nodename'], '\\.'))[1]
  PID = Sys.getpid()
  tempString = sprintf('%s-%s', NODE, PID)
  return(tempString)
}


#temporary directory
function.getTempDir = function(useRAM = TRUE, usecache = TRUE) {
  
  currentNode = unlist(strsplit(Sys.info()['nodename'], '\\.'), use.names = FALSE)[1]
  
  
  if (usecache) {
    tempDirPath = '/'
    return(tempDirPath)
  }
  
  if (useRAM) {
    tempDirPath = '/the_temp_path/'
    return(tempDirPath)
    
  } else {
    
    tempDirPath = '/'
    return(tempDirPath)
  }
}

function.freadXZ = function(file, sep = '\t', header = TRUE) {
  stopifnot(is.character(file))
  require(data.table)
  
  tempOutDir = sprintf('%sfunction.freadXZ_%s/', function.getTempDir(useRAM = TRUE, usecache = TRUE), function.getPID())
  dir.create(tempOutDir, showWarnings = FALSE, recursive = TRUE)
  
  secondStamp = format(Sys.time(), format = '%s')
  fileID = .FileNameScramble(file)
  tempFilePath = sprintf('%s%s_%s.TMP', tempOutDir, secondStamp, fileID)
  xzCommand = sprintf('xz -d -c %s > %s', file, tempFilePath)
  system(command = xzCommand, ignore.stdout = FALSE, ignore.stderr = TRUE)

  tempDT = fread(file = tempFilePath, sep = sep, header = header, showProgress = FALSE)
  
  file.remove(tempFilePath)
  
  return(tempDT)
}

.FileNameScramble = function(charVec) {
  require(stringi)
  
  charVec = stri_trans_tolower(charVec)
  
  charVec = tail(unlist(strsplit(charVec, split = '/'), use.names = FALSE), n = 1)
  
  charVec = unlist(strsplit(charVec, split = ''), use.names = FALSE)
  
  charVec = paste0(na.omit(match(charVec, letters)), collapse = '')
  
  return(charVec)
}


function.XZSaveRDS = function(obj, file, threads = 32, compression = 6) {
  
  stopifnot(is.character(file))
  
  function.remove(paths = file, check = TRUE)
  
  xzCommand = sprintf('xz -z -T %s -%s > %s', threads, compression, file)
  
  xzConnection = pipe(description = xzCommand, open = 'wb')
  
  saveRDS(object = obj, file = xzConnection)
  
  close(xzConnection)
}