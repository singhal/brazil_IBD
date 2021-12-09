plot.bammdata_sonal <- function (x, tau = 0.01, method = "phylogram", xlim = NULL, ylim = NULL, vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, legend = FALSE, spex = "s", lwd = 1, cex = 1, pal = "RdYlBu", mask = integer(0), mask.color = gray(0.5), colorbreaks = NULL, logcolor = FALSE, breaksmethod = "linear", color.interval = NULL, JenksSubset = 20000, par.reset = FALSE, direction = "rightwards", ...) {
  if (inherits(x, "bammdata")) {
    if (attributes(x)$order != "cladewise") {
      stop("Function requires tree in 'cladewise' order");
    }
    phy <- BAMMtools:::as.phylo.bammdata(x);
  }
  else stop("Object ephy must be of class bammdata");
  
  if (!spex %in% c('s','e','netdiv')) {
    stop("spex must be 's', 'e' or 'netdiv'.");
  }
  # 
  # if (length(pal) == 1 && !pal %in% names(get("palettes", envir=.colorEnv)) && pal != "temperature" && pal != "terrain")
  #   pal <- rep(pal, 3)
  # else if (length(pal) == 2)
  #   pal <- c(pal, pal[2]);
  
  if (breaksmethod == 'linear' & !is.null(color.interval)) {
    if (length(color.interval) != 2) {
      stop("color.interval must be a vector of 2 numeric values.");
    }
  }
  
  if (!is.binary.phylo(phy)) {
    stop("Function requires fully bifurcating tree");
  }
  if (any(phy$edge.length == 0)) {
    warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero");
  }
  if (!("dtrates" %in% names(x))) {
    x <- BAMMtools:::dtRates(x, tau);
  }
  
  NCOLORS <- 64;
  
  if (!is.null(color.interval)) {
    # change the number of breaks such that the range of color.interval 
    # is equivalent in terms of number of colors to the full range
    # this way we preserve good resolution
    # Here we will ensure that NCOLORS bins occur within the color.interval
    
    if (x$type == "trait") {
      ratesRange <- range(x$dtrates$rates);
    } else if (x$type == "diversification") {
      if (tolower(spex) == "s") {
        ratesRange <- range(x$dtrates$rates[[1]]);
      } else if (tolower(spex) == "e") {
        ratesRange <- range(x$dtrates$rates[[2]]);
      } else if (tolower(spex) == "netdiv") {
        ratesRange <- range(x$dtrates$rates[[1]] - x$dtrates$rates[[2]]);
      }
    }	
    if (all(!is.na(color.interval))) {
      brks <- seq(min(color.interval[1], ratesRange[1]), max(color.interval[2], ratesRange[2]), length.out = (NCOLORS+1));
      intervalLength <- length(which.min(abs(color.interval[1] - brks)) : which.min(abs(color.interval[2] - brks)));
    } else if (is.na(color.interval[1])) {
      brks <- seq(ratesRange[1], max(color.interval[2], ratesRange[2]), length.out = (NCOLORS+1));
      intervalLength <- length(1 : which.min(abs(color.interval[2] - brks)));
    } else if (is.na(color.interval[2])) {
      brks <- seq(min(color.interval[1], ratesRange[1]), ratesRange[2], length.out = (NCOLORS+1));
      intervalLength <- length(which.min(abs(color.interval[1] - brks)) : length(brks));
    }
    NCOLORS <- round((NCOLORS ^ 2) / intervalLength)    	
  }
  
  if (is.null(colorbreaks)) {
    colorbreaks <- BAMMtools:::assignColorBreaks(x$dtrates$rates, NCOLORS, spex, logcolor, breaksmethod, JenksSubset);
  }
  if (x$type == "trait") {
    colorobj <- BAMMtools:::colorMap(x$dtrates$rates, pal, colorbreaks, logcolor, color.interval);
  }
  else if (x$type == "diversification") {
    if (tolower(spex) == "s") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[1]], pal, colorbreaks, logcolor, color.interval);
    }
    else if (tolower(spex) == "e") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[2]], pal, colorbreaks, logcolor, color.interval);
    }
    else if (tolower(spex) == "netdiv") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], pal, colorbreaks, logcolor, color.interval);
    }
  }
  else {
    stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'");	
  }
  edge.color <- colorobj$cols;
  keep = list()
  #    if (is.ultrametric(phy))    
  #    	tH <- max(branching.times(phy))
  #    else
  #    	tH <- max(NU.branching.times(phy));
  tH <- max(x$end);
  phy$begin <- x$begin;
  phy$end <- x$end;
  tau <- x$dtrates$tau;
  if (method == "polar") {
    ret <- BAMMtools:::setPolarTreeCoords(phy, vtheta, rbf);
    rb <- tH * rbf;
    p <- BAMMtools:::mkdtsegsPolar(ret$segs[-1,], tau, x$edge);
  }
  else if (method == "phylogram") {
    ret <- BAMMtools:::setPhyloTreeCoords(phy);
    p <- BAMMtools:::mkdtsegsPhylo(ret$segs[-1,], tau, x$edge);
  }
  else {
    stop("Unimplemented method");
  }
  x0 <- c(ret$segs[1,1], p[, 1]);
  x1 <- c(ret$segs[1,3], p[, 2]);
  y0 <- c(ret$segs[1,2], p[, 3]);
  y1 <- c(ret$segs[1,4], p[, 4]);
  offset <- table(p[, 5])[as.character(unique(p[, 5]))];
  if (length(mask)) {
    edge.color[p[,5] %in% mask] <- mask.color;
  }
  arc.color <- c(edge.color[1], edge.color[match(unique(p[, 5]), p[, 5]) + offset]);
  edge.color <- c(edge.color[1], edge.color);
  
  xy = data.frame(x0, x1, y0, y1)
  keep = list(edge.color, ret, rb, phy$end, tH, arc.color, xy)
  
  if (show) {
    op <- par(no.readonly = TRUE);
    if (length(list(...))) {
      par(...);
    }
    if (legend) {
      #par(fig=c(0,0.9,0,1));
      par(mar = c(5, 4, 4, 5))
    }
    plot.new();
    ofs <- 0;
    if (labels) {
      if (method == "phylogram")
        ofs <- max(nchar(phy$tip.label) * 0.03 * cex * tH)
      else
        ofs <- max(nchar(phy$tip.label) * 0.03 * cex);
    }
    if (method == "polar") {
      if (is.null(xlim) || is.null(ylim)) {
        if (is.null(xlim))
          xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs)
        if (is.null(ylim))
          ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs) 
      }
      plot.window(xlim = xlim, ylim = ylim, asp = 1);
      segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, lend = 2);
      arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + phy$end/tH), border = arc.color, lwd = lwd);
      
      if (labels) {
        for (k in 1:length(phy$tip.label)) {
          text(ret$segs[-1, ][phy$edge[, 2] == k, 3],ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k],cex = cex, srt = (180/pi) * ret$arcs[-1,][phy$edge[, 2] == k, 1], adj = c(0, NA));
        }
      }
    }
    if (method == "phylogram") {
      direction <- match.arg(direction, c("rightwards","leftwards","downwards","upwards"));
      if (direction == "rightwards") {
        bars <- redirect(cbind(x0,y0,x1,y1),0);
        arcs <- redirect(ret$arcs,0);
        bars[,c(1,3)] <- tH * bars[,c(1,3)];
        arcs[,c(1,3)] <- tH * arcs[,c(1,3)];
        
        # xlim <- c(0, 1 + ofs);
        # ylim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
        
        ret$segs[-1, c(1,3)] <- tH * ret$segs[-1, c(1,3)]; 
        
      }
      else if (direction == "leftwards") {
        bars <- redirect(cbind(x0,y0,x1,y1),pi);
        bars[,c(2,4)] <- abs(bars[,c(2,4)]);
        arcs <- redirect(ret$arcs,pi);
        arcs[,c(2,4)] <- abs(arcs[,c(2,4)]);
        
        
        bars[,c(1,3)] <- tH * bars[,c(1,3)];
        arcs[,c(1,3)] <- tH * arcs[,c(1,3)];
        
        
        ret$segs[-1, c(1,3)] <- -tH * ret$segs[-1, c(1,3)];
        
        # xlim <- rev(-1*c(0, 1 + ofs));
        # ylim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
      }
      else if (direction == "downwards") {
        bars <- BAMMtools:::redirect(cbind(x0,y0,x1,y1),-pi/2);
        arcs <- BAMMtools:::redirect(ret$arcs,-pi/2);
        
        bars[,c(2,4)] <- tH * bars[,c(2,4)];
        arcs[,c(2,4)] <- tH * arcs[,c(2,4)];
        
        
        ret$segs <- BAMMtools:::redirect(ret$segs, -pi/2);
        ret$segs[,c(2,4)] <- tH * ret$segs[,c(2,4)];
        
        # xlim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
        # ylim <- rev(-1*c(0, 1 + ofs));	
      }
      else if (direction == "upwards") {
        bars <- BAMMtools:::redirect(cbind(x0,y0,x1,y1),pi/2);
        bars[,c(1,3)] <- abs(bars[,c(1,3)]);
        arcs <- BAMMtools:::redirect(ret$arcs,pi/2);
        arcs[,c(1,3)] <- abs(arcs[,c(1,3)]);
        
        bars[,c(2,4)] <- tH * bars[,c(2,4)];
        arcs[,c(2,4)] <- tH * arcs[,c(2,4)];
        
        ret$segs <- BAMMtools:::redirect(ret$segs, pi/2);
        ret$segs[,c(1,3)] <- abs(ret$segs[,c(1,3)]);
        ret$segs[,c(2,4)] <- tH * ret$segs[,c(2,4)];
        
        # xlim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
        # ylim <- c(0, 1 + ofs);
      }
      if (is.null(xlim) && direction == "rightwards") xlim <- c(0, tH + ofs);
      if (is.null(xlim) && direction == "leftwards") xlim <- c(-(tH + ofs), 0);
      if (is.null(ylim) && (direction == "rightwards" || direction == "leftwards")) ylim <- c(0, phy$Nnode);  
      
      if (is.null(xlim) && (direction == "upwards" || direction == "downwards")) xlim <- c(0, phy$Nnode);
      if (is.null(ylim) && direction == "upwards") ylim <- c(0, tH + ofs);
      if (is.null(ylim) && direction == "downwards") ylim <- c(-(tH + ofs), 0);  
      
      
      plot.window(xlim = xlim, ylim = ylim);
      segments(bars[-1,1], bars[-1,2], bars[-1,3], bars[-1,4], col = edge.color[-1], lwd = lwd, lend = 2);
      isTip <- phy$edge[, 2] <= phy$Nnode + 1;
      isTip <- c(FALSE, isTip);
      segments(arcs[!isTip, 1], arcs[!isTip, 2], arcs[!isTip, 3], arcs[!isTip, 4], col = arc.color[!isTip], lwd = lwd, lend = 2);  
      if (labels) {
        if (direction == "rightwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 4, offset = 0.25)
        else if (direction == "leftwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 2, offset = 0.25)
        else if (direction == "upwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 4, srt = 90, offset = 0)
        else if (direction == "downwards")
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 2, srt = 90, offset = 0);
      }
    }
    # if (legend) {
    # #rateLegend(colorobj$colsdensity, logcolor);
    # if (is.null(color.interval)) {
    # barLegend(pal, colorbreaks, fig=c(0.9,1,0.25,0.75), side=2);
    # } else {
    # barLegend(pal, colorbreaks, fig=c(0.9,1,0.25,0.75), side=2, colpalette=colorobj$colpalette);
    # }
    # }
  }
  index <- order(as.numeric(rownames(ret$segs)));
  if (show) {
    if (method == "phylogram") {
      assign("last_plot.phylo", list(type = "phylogram", direction = direction, Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
    } else if (method == "polar") {
      assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 3], yy = ret$segs[index, 4], theta = ret$segs[index, 5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
    }
    if (legend) {
      BAMMtools:::addBAMMlegend(x = list(coords = ret$segs[-1, ], colorbreaks = colorobj$breaks, palette = colorobj$colpalette, colordens = colorobj$colsdensity), location = 'right')
    }
  }
  if (par.reset) {
    par(op);
  }
  invisible(list(coords = ret$segs[-1, ], colorbreaks = colorobj$breaks, palette = colorobj$colpalette, colordens = colorobj$colsdensity, plotvals = keep));
}

addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}


old_colorMap = function(x, pal, NCOLORS)
{
  dpal = c('BrBG','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
  colset = numeric(length(x));
  if(length(pal) == 3)
  {
    colpalette = colorRampPalette(pal,space='Lab')(NCOLORS);	
  }
  else if(pal %in% dpal)
  {
    colpalette = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
  }
  else if(pal == 'temperature')
  {
    colpalette = gplots::rich.colors(NCOLORS);	
  }
  
  bks = quantile(x, seq(0,1,length.out=(NCOLORS+1)));
  for(i in 2:length(bks))
  {
    if(i == 2)
    {
      colset[x < bks[2]] = colpalette[1];	
    }
    else if(i == length(bks))
    {
      colset[x >= bks[length(bks)-1]] = colpalette[length(bks)-1];
    }
    else
    {
      colset[x >= bks[i-1] & x < bks[i]] = colpalette[i-1];
    }
  }
  return(colset);
}

addBAMMlegend_sonal <- function(colorbreaks, pal, direction, side, location = 'topleft', nTicks = 2, stretchInterval = FALSE, shortFrac = 0.02, longFrac = 0.3, axisOffset = 0.002, cex.axis = 0.8, labelDist = 0.7, ...) {
  #location xmin,xmax,ymin,ymax
  
  if (hasArg('corners')) {
    stop('Error: some options have been deprecated. Please consult the documentation.')
  }
  
  if(!hasArg('direction')) {
    direction <- 'auto'
  }
  
  if (!direction %in% c('auto', 'vertical', 'horizontal')) {
    stop("direction must be auto, vertical or horizontal.");
  }
  
  if (is.character(location)) {
    if (!location %in% c('bottomleft','bottomright','topleft','topright','bottom','top','left','right')) {
      stop('location is not recognized.');
    }
  }
  
  # If there are duplicate colors, then this color ramp is the result of
  # a specified color.interval. If stretchInterval is TRUE, then rather
  # than include many duplicate colors, we will only include the color.interval
  # range.
  # intervalSide = top means that the top range of the color palette has
  # duplicate colors
  intervalSide <- NULL
  if (length(unique(pal)) != length(pal) & stretchInterval) {
    uniquePal <- which(!duplicated(pal))
    if (uniquePal[2] != (uniquePal[1] + 1)) {
      uniquePal[1] <- uniquePal[2] - 1
    }
    colorbreaks <- colorbreaks[c(uniquePal, uniquePal[length(uniquePal)] + 1)]
    pal <- pal[uniquePal]
    if (identical(x$palette[1], x$palette[2]) & identical(tail(x$palette, 1), tail(x$palette, 2)[1])) {
      intervalSide <- 'both'
    } else if (identical(x$palette[1], x$palette[2]) & !identical(tail(x$palette, 1), tail(x$palette, 2)[1])) {
      intervalSide <- 'bottom'
    } else if (!identical(x$palette[1], x$palette[2]) & identical(tail(x$palette, 1), tail(x$palette, 2)[1])) {
      intervalSide <- 'top'
    }
  } 
  
  n <- length(colorbreaks);
  
  #return plot region extremes and define outer coordinates
  minX <- grconvertX(par('fig')[1], from = 'ndc', to = 'user') 
  maxX <- grconvertX(par('fig')[2], from = 'ndc', to = 'user')
  minY <- grconvertY(par('fig')[3], from = 'ndc', to = 'user')
  maxY <- grconvertY(par('fig')[4], from = 'ndc', to = 'user')
  
  xrange <- maxX - minX
  yrange <- maxY - minY
  minX <- minX + xrange * 0.05
  maxX <- maxX - xrange * 0.05
  minY <- minY + yrange * 0.05
  maxY <- maxY - yrange * 0.05
  
  if (is.character(location)) {
    
    if (location == 'topleft' & direction %in% c('auto', 'vertical')) {
      location <- vector('numeric', length = 4);
      location[1] <- minX
      location[2] <- minX + (maxX - minX) * shortFrac
      location[3] <- maxY - (maxY - minY) * longFrac
      location[4] <- maxY
    } else
      
      if (location == 'topleft' & direction == 'horizontal') {
        location <- vector('numeric', length = 4);
        location[1] <- minX
        location[2] <- minX + (maxX - minX) * longFrac
        location[3] <- maxY - (maxY - minY) * shortFrac
        location[4] <- maxY
      } else
        
        if (location == 'topright' & direction %in% c('auto', 'vertical')) {
          location <- vector('numeric', length = 4);
          location[1] <- maxX - (maxX - minX) * shortFrac
          location[2] <- maxX
          location[3] <- maxY - (maxY - minY) * longFrac
          location[4] <- maxY
        } else
          
          if (location == 'topright' & direction == 'horizontal') {
            location <- vector('numeric', length = 4);
            location[1] <- maxX - (maxX - minX) * longFrac
            location[2] <- maxX
            location[3] <- maxY - (maxY - minY) * shortFrac
            location[4] <- maxY
          } else
            
            if (location == 'bottomleft' & direction %in% c('auto', 'vertical')) {
              location <- vector('numeric', length = 4);
              location[1] <- minX
              location[2] <- minX + (maxX - minX) * shortFrac
              location[3] <- minY
              location[4] <- minY + (maxY - minY) * longFrac
            } else
              
              if (location == 'bottomleft' & direction == 'horizontal') {
                location <- vector('numeric', length = 4);
                location[1] <- minX
                location[2] <- minX + (maxX - minX) * longFrac
                location[3] <- minY
                location[4] <- minY + (maxY - minY) * shortFrac
              } else
                
                if (location == 'bottomright' & direction %in% c('auto', 'vertical')) {
                  location <- vector('numeric', length = 4);
                  location[1] <- maxX - (maxX - minX) * shortFrac
                  location[2] <- maxX
                  location[3] <- minY
                  location[4] <- minY + (maxY - minY) * longFrac
                } else
                  
                  if (location == 'bottomright' & direction == 'horizontal') {
                    location <- vector('numeric', length = 4);
                    location[1] <- maxX - (maxX - minX) * longFrac
                    location[2] <- maxX
                    location[3] <- minY
                    location[4] <- minY + (maxY - minY) * shortFrac 
                  } else
                    
                    if (location == 'left') {
                      location <- vector('numeric', length = 4);
                      location[1] <- minX
                      location[2] <- minX + (maxX - minX) * shortFrac
                      location[3] <- mean(par('usr')[3:4]) - ((maxY - minY) * longFrac)/2
                      location[4] <- mean(par('usr')[3:4]) + ((maxY - minY) * longFrac)/2
                      direction <- 'vertical'
                    } else
                      
                      if (location == 'right') {
                        location <- vector('numeric', length = 4);
                        location[1] <- maxX - (maxX - minX) * shortFrac
                        location[2] <- maxX
                        location[3] <- mean(par('usr')[3:4]) - ((maxY - minY) * longFrac)/2
                        location[4] <- mean(par('usr')[3:4]) + ((maxY - minY) * longFrac)/2
                        direction <- 'vertical'
                      } else
                        
                        if (location == 'top') {
                          location <- vector('numeric', length = 4);
                          location[1] <- mean(par('usr')[1:2]) - ((maxX - minX) * longFrac)/2
                          location[2] <- mean(par('usr')[1:2]) + ((maxX - minX) * longFrac)/2
                          location[3] <- maxY - (maxY - minY) * shortFrac
                          location[4] <- maxY
                          direction <- 'horizontal'
                        } else
                          
                          if (location == 'bottom') {
                            location <- vector('numeric', length = 4);
                            location[1] <- mean(par('usr')[1:2]) - ((maxX - minX) * longFrac)/2
                            location[2] <- mean(par('usr')[1:2]) + ((maxX - minX) * longFrac)/2
                            location[3] <- minY
                            location[4] <- minY + (maxY - minY) * shortFrac
                            direction <- 'horizontal'
                          }
  }
  
  
  # infer direction based on dimensions of legend box
  if (direction == 'auto') {
    if (((location[2] - location[1]) / (par('usr')[2] - par('usr')[1])) >= ((location[4] - location[3]) / (par('usr')[4] - par('usr')[3]))) {
      direction <- 'horizontal';
    } else {
      direction <- 'vertical';
    }
  }
  
  if (direction == 'horizontal') {
    axisOffset <- axisOffset * (par('usr')[4] - par('usr')[3]);
  } else if (direction == 'vertical') {
    axisOffset <- axisOffset * (par('usr')[2] - par('usr')[1]);
  }
  
  #determine side for labels based on location in plot and direction
  if (!hasArg('side')) {
    if (direction == 'vertical') { #side = 1 or 4
      if (mean(location[1:2]) <= mean(par('usr')[1:2])) {
        side <- 4;
      } else {
        side <- 2;
      }
    }
    if (direction == 'horizontal') { #side = 2 or 3
      if (mean(location[3:4]) > mean(par('usr')[3:4])) {
        side <- 1;
      } else {
        side <- 3;
      }
    }
  }
  
  if (direction == 'horizontal') {
    x <- seq(from = location[1], to = location[2], length.out = n);
    width <- location[3:4];
  } else {
    x <- seq(from = location[3], to = location[4], length.out = n);
    width <- location[1:2];
  }
  
  #get bin coordinates
  x <- rep(x,each = 2);
  x <- x[-c(1,length(x))];
  x <- matrix(x, ncol = 2, byrow = TRUE);
  
  #find tick locations
  #get equivalent rate bins
  z <- rep(colorbreaks,each = 2);
  z <- z[-c(1,length(z))];
  z <- matrix(z, ncol = 2, byrow = TRUE);
  
  tx <- trunc(seq(from = 1, to = nrow(x), length.out = nTicks + 2));
  tickLocs <- x[tx,1]
  tx <- z[tx,1]
  tickLocs[length(tickLocs)] <- max(x[,2])
  tx[length(tx)] <- max(z[,2])	
  
  #plot bar
  if (direction == 'horizontal') {
    rect(xleft = x[,1], ybottom = width[1], xright = x[,2], ytop = width[2], border = pal, col = pal, xpd = NA);
  } else {
    rect(xleft = width[1], ybottom = x[,1], xright = width[2], ytop = x[,2], border = pal, col = pal, xpd = NA);
  }
  
  #add tickmarks
  tickLabels <- as.character(signif(tx, 2));
  if (!is.null(intervalSide)) {
    if (intervalSide == 'top' | intervalSide == 'both') {
      tickLabels[length(tickLabels)] <- paste('\u2265', tickLabels[length(tickLabels)])
    }
    if (intervalSide == 'bottom' | intervalSide == 'both') {
      tickLabels[1] <- paste('\u2264', tickLabels[1])
    }
  }
  if (side == 1) { #bottom
    axis(side, at = tickLocs, pos = location[3] - axisOffset, labels = tickLabels, xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  } 
  if (side == 3) { #top
    axis(side, at = tickLocs, pos = location[4] + axisOffset, labels = tickLabels, xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  }
  if (side == 2) { #left
    axis(side, at = tickLocs, pos = location[1] - axisOffset, labels = tickLabels, xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  }
  if (side == 4) { #right
    axis(side, at = tickLocs, pos = location[2] + axisOffset, labels = tickLabels, xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
  }
  invisible(list(coords = x, width = width, pal = pal, tickLocs = tickLocs, labels = tx))
}