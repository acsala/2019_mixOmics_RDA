

check_specs <-
function(blocks, scaling, modes, scheme, scaled, tol, maxiter, plscomp)
{
  check_scale = check_scaling(scaling, scaled, blocks)
  
  # output
  list(scaling = check_scale$scaling,
       modes = check_modes(modes, blocks),
       scheme = check_scheme(scheme),
       scaled = check_scale$scaled,
       tol = check_tol(tol),
       maxiter = check_maxiter(maxiter),
       plscomp = check_plscomp(plscomp, scaling, modes))  
}


check_scaling <- function(scaling, scaled, blocks)
{
  # if scaling is present
  if (!is.null(scaling)) 
  {
    # make sure scaling is a list
    if (!is.list(scaling))
      stop("\nInvalid 'scaling'. Must be a list.")
    # compatibility between blocks and scaling
    if (length(blocks) != length(scaling)) {
      stop("\nLength of 'scaling' differs from length of 'blocks'.")
    }
    if (!identical(lengths(blocks), lengths(scaling))) {
      stop("\nLengths of 'scaling' differs from lengths of 'blocks'.")
    }
    
    # string manipulation of elements in 'scaling'
    scaling_aux = tolower(unlist(scaling))
    scaling_aux = substr(scaling_aux, start=1, stop=3)
    
    # are there any unrecognized scaling types?
    bad_scale <- !(scaling_aux %in% c("num", "raw", "ord", "nom"))
    if (any(bad_scale)) {
      bad = unlist(scaling)[bad_scale]
      stop(sprintf("\nSorry. Unrecognized scaling type: '%s'", bad))
    }
    
    # set all numeric when mixing only 'num' and 'raw'
    if (!any(scaling_aux %in% c('ord', 'nom'))) {
      num_and_raw = intersect(scaling_aux, c("num", "raw"))
      if (length(num_and_raw) > 1) {
        scaling = lapply(blocks, function(x) rep("num", length(x)))
      }
      if (all(unlist(scaling) == "num")) scaled = TRUE
      if (all(unlist(scaling) == "raw")) scaled = FALSE    
    }
    # final scaling
    scaling = lapply(scaling, function(x) substr(tolower(x), start=1, stop=3))
  } else {
    if (!is.logical(scaled)) scaled = TRUE
  }
  
  # output
  list(scaling = scaling, scaled = scaled)
}


check_modes <- function(modes, blocks)
{
  # default modes    
  if (is.null(modes)) {
    modes = rep("predictive", length(blocks))
  } 
  
  if (length(blocks) != length(modes)) {
    warning("Warning: length of 'modes' different from length of 'blocks'")
    message("Default modes 'predictive' is used")
    modes = rep("predictive", length(blocks))
  }
  
  # are there any unrecognized modes?
  bad_modes <- !(modes %in% c("predictive", "predicted", "NEWA", "PLSCOW", "PLSCORE"))
  if (any(bad_modes)) {
    bad = modes[bad_modes]
    stop(sprintf("\nSorry. Unrecognized mode: '%s'", bad))
  }
    
  # cannot mix modes "A" and "newA"
  mixed_modes = intersect(modes, c("predictive", "NEWA"))
  if (length(mixed_modes) > 1) {
    stop("\nSorry. Can't work with both modes 'A' and 'newA'")
  }  

  # output
  modes
}


check_scheme <- function(scheme)
{
  # some string manipulations
  if (!is.character(scheme)) scheme = as.character(scheme)
  scheme = tolower(scheme)
    
  SCHEMES = c("centroid", "factorial", "path")
  scheme_match = pmatch(scheme, SCHEMES)
  if (is.na(scheme_match)) {
    warning("\nInvalid 'scheme'. Default 'scheme=centroid' is used.")   
    scheme = "centroid"
  } else {
    scheme = SCHEMES[scheme_match]
  }
  
  # output
  scheme
}


check_tol <- function(tol)
{
  if (mode(tol) != "numeric" || length(tol) != 1L || 
        tol <= 0 || tol > 0.001) {
    warning("Warning: Invalid 'tol'. Default 'tol=0.000001' is used.")   
    tol = 0.000001
  } 
  
  # output
  tol
}


check_maxiter <- function(maxiter)
{
  if (!is_positive_integer(maxiter) || maxiter < 100) {
    warning("Warning: Invalid 'maxiter'. Default 'maxiter=100' is used.")   
    maxiter = 100
  } 
  
  # output
  maxiter
}


check_plscomp <- function(plscomp, scaling, modes)
{
  if (is.null(scaling)) plscomp = NULL
  
  if (!is.null(scaling)) {
    if (is.null(plscomp)) {
      plscomp = rep(1, length(scaling))      
    } else {
      if (any(!is_positive_integer(plscomp)))
        stop("\n'plscomp' must be an integer vector")
      if (length(scaling) != length(plscomp))
        stop("\nlength of 'plscomp' differs from number of blocks")
      
      plscomp = as.integer(plscomp)
      for (j in 1:length(plscomp)) 
      {
        # plscomp's cannot exceed number of variables in block j
        if (plscomp[j] > length(scaling[[j]]))
        {
          stop(sprintf("%s %d %s", "element", j, 
                       "in 'plscomp' exceeds number of variables"))          
        }
        # make sure plscomp[j]=1 when mode "NEWA"
        if (modes[j] == "NEWA") plscomp[j] = 1
      }
    }
  }
  
  # output
  plscomp
}
