


# assignSpanningSet
# Find spanning taxon pairs for each node,
#    required to create event data object

assignSpanningSet <- function(phy, slow = FALSE){
	
	nodeset <- character(max(phy$edge))
	nodeset[1:length(phy$tip.label)] <- phy$tip.label
	ntips <- length(phy$tip.label)
	rootnode <- ntips + 1
	
	spanmat <- matrix(NA, nrow=length(nodeset), ncol=2)
	spanmat[1:ntips,1] <- phy$tip.label
	
	if (slow){
		for (ii in rootnode:nrow(phy$edge)){
			xc <- extract.clade(phy, node = ii)$tip.label
			spanmat[ii,] <- xc[c(1, length(xc))]
		}
		
	}else{
		for (ii in 1:length(phy$tip.label)){
			tipname <- phy$tip.label[ii]
			curr_node <- ii
			DONE <- FALSE
			while (curr_node != rootnode & !DONE){
				anc <- phy$edge[,1][phy$edge[ , 2] == curr_node]
				if (nodeset[anc] == ""){
					nodeset[anc] <- tipname
					curr_node <- anc
				}else{
					DONE <- TRUE
				} 
			}
		}		
		for (ii in rootnode:max(phy$edge)){
			dset <- nodeset[phy$edge[,2][phy$edge[,1] == ii]]
			spanmat[ii, ] <- dset
		}
	
	}
 
	return(spanmat)
	
}

# clads_to_eventdata
#	args:
#	map_clads: clads output object w stored MAP rates (as provided)
#   delta: must be smaller than the smallest branch length in tree
#   Tree must be fully resolved; no zero length branches.
 

clads_to_eventdata <- function(map_clads, delta = 0.00001){
	
	tree <- as.phylo(map_clads$tree)
	tree <- check_and_fix_ultrametric(tree)
  
	phy <- BAMMtools:::getStartStopTimes(tree)

	spanmat <- assignSpanningSet(phy)
	
	generation <- rep(1, nrow(spanmat))
	leftchild  <- spanmat[,1]
	rightchild <- spanmat[,2]
	
	abstime     <- numeric(nrow(spanmat))
	lambdainit  <- numeric(nrow(spanmat)) 
	lambdashift <- numeric(nrow(spanmat))
	muinit      <- numeric(nrow(spanmat))
	mushift     <- numeric(nrow(spanmat))
	
	# iterate over nodes, 
	#   fill in rates times etc
	for (ii in 1:nrow(phy$edge)){
 		fnode <- phy$edge[ii,2]
 		lambdainit[fnode] <- map_clads$lambdai_map[ii]
 		abstime[fnode]    <- phy$begin[ii] + delta
	}
	
	# The root process is still not filled in;
	#   just take mean rate for 2 desc branches:
	rootnode <- length(phy$tip.label) + 1
	root_rates <- map_clads$lambdai_map[phy$edge[,1] == rootnode]
	lambdainit[rootnode] <- mean(root_rates)
	if (length(map_clads$eps_map == 1)){
		muinit <- map_clads$eps_map * lambdainit
	}else{
		stop("branch-specific epsilon currently not supported\n")
	}
	
	dff <- data.frame(generation, leftchild, rightchild, abstime, lambdainit, lambdashift, muinit, mushift, stringsAsFactors=F)
	
	dff <- dff[order(dff$abstime), ]
	
	phy <- read.tree(text = write.tree(ladderize(tree)))
	
	ed <- getEventData(phy, eventdata = dff, burnin=0, nsamples=1)
	return(ed)
}


# check_and_fix_ultrametric
#     checks trees to see if they pass ape ultrametricity test.
# If not, it computes the differential root-to-tip distance across all tips.
# It adds the appropriate quantity to each terminal branch length to ensure that 
# tree passes ultrametric test.
# Note: this is only a valid method of making trees ultrametric when the 
# 	non-ultrametricity is due to small numerical discrepancies, e.g., 
#   rounding or other floating point issues during phylogeny construction.
# 

check_and_fix_ultrametric <- function(phy){
	
	if (!is.ultrametric(phy)){
		
		vv <- vcv.phylo(phy)
		dx <- diag(vv)
		mxx <- max(dx) - dx
		for (i in 1:length(mxx)){
			phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
		}
		if (!is.ultrametric(phy)){
			stop("Ultrametric fix failed\n")
		}	
	}
	
	return(phy)
}




