#Calculates species-specific diversification rate, as per Jetz et al. 2012
#returns a vector of diversification rates, each element named by tip label.
#author: Pascal Title
#date: 31 Jan 2013

# for each species i,
	# get branch lengths from root to tip
 		# sum(edge length j * (1/(2 ^ (j-1)))) ^ (-1)


jetzDivRates <- function(tree) {
	
	spRate <- function(sp, tree) {
		#get branch lengths from root to tip
		edges <- vector()
		daughterNode <- match(sp, tree$tip.label)
		while (daughterNode != (length(tree$tip.label) + 1)) {
			parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
			edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
			daughterNode <- parentNode
		}
		
		res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
		res <- res ^ (-1)
		
		return(res)
	}
	
	rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
	names(rates) <- tree$tip.label
	
	return(rates)
	
}







