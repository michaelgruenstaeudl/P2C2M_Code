gene.tree.prob <-
function(sptree, gtree, association, ploidy=1){
	
	species<-sptree$tip.label
	nspec<-length(species)
	nbranch<-(2*nspec)-1

	tiplist<-list()
	for(i in 1:nspec){
		tiplist[[species[i]]]<-association[which(association[,1]==species[i]),2]
		}

	####this subfunction returns the tip numbers for a given node
	node.tips<-function (phy, node){
 		n<- length(phy$tip.label)
 		if (node<= n){node}else{
 			l<- numeric()
 			d<- phy$edge[which(phy$edge[,1]==node),2]
 			for(j in d){if(j<= n){l<- c(l, j)}else{l<-c(l, node.tips(phy,j))}}
       		l}}
	#####end subfunction
	####subfunction that gets a tree containing descendants of all gene lineages that pass through a branch 
	######do I really need that node conditional??? or should I just deal with it later?
	get.lin<-function(sptree, gtree, association, node){	
		species<-sptree$tip.label
		if((length(species)+1)!=node){
			spec.desc<-species[node.tips(sptree, node)]
			gtree.desc<-c()
			for(i in 1:length(spec.desc)){
				gtree.desc<-c(gtree.desc, association[which(association[,1]==spec.desc[i]),2])
				}
			match(gtree.desc, gtree$tip.label)->index
			pruned<-drop.tip(gtree, gtree$tip.label[-index])
			}else{pruned<-gtree}		
		return(pruned)
		}
	#####end subfunction
	#####prepare branch demographic info...
	demo.maker<-function(sptree, nspec){
		demo<-cbind(sptree$edge[,2], sptree$edge.length)
		demo<-rbind(demo, c((nspec+1), Inf))
		demo<-cbind(demo, sptree$dmv)
		demo<-demo[order(demo[,1]),]
		sbt<-branching.times(sptree)
		sbt<-sbt[order(as.numeric(names(sbt)))]
		sbt<-c(rep(0, nspec), sbt)
		demo<-cbind(demo, sbt)
		colnames(demo)<-c("node", "length", "dmv", "sbt")
		rownames(demo)<-c(1:length(sbt))
		return(demo)
		}
	#####


	###this function calculates the probability of a the coalescent process within a branch of the species tree
	
	b.prob<-function(sptree, gtree, demo, node){
		#first, we get all descendants of all lineages that pass through the specified branch (node)
		pruned<-get.lin(sptree, gtree, association, node)
		#then we get the branching times for all descendents and sort them
		gbt<-sort(branching.times(pruned))
		#we add a time '0' to the starting point
		gbt<-c(0, gbt)
		#we make the vector a matrix, with col1 a series of times and col2 the number of lineages that exist at that time point
		gbt<-cbind(gbt, length(gbt):1)
		#get the starting point of the branch
		start<-demo[node,"sbt"]
		#get the end point of the branch
		end<-(demo[node, "sbt"]+demo[node, "length"])
		#get how many lineages enter the branch
		enter<-gbt[gbt[,1]==max(gbt[gbt[,1]<=start,1]),2]
		#get how many lineages exit the branch
		exit<-gbt[gbt[,1]==max(gbt[gbt[,1]<=end,1]),2]
		
		####if there's only one lineage in the branch, exit with probability 1. 
		if(start==1){
			return(1)
			}
		#### if more, calculate coalescent probability. 
		else{
			#first, toss out all coalescences occurring before or after the branch under consideration
			gbt<-gbt[gbt[,1]>=start & gbt[,1]<=end,]
			#ensure gbt matrix row 1 is (start,enter), and the last row is (end, exit)
			if(start!=0){	
				gbt<-rbind(c(start, enter), gbt, c(end, exit))
				}else{gbt<-rbind(gbt, c(end, exit))}

			#now we turn gbt into a matrix ("waits") containing waiting times between coalescent events (and the end of the branch)
				#and the number of lineages in each waiting time. 				
			waits<-gbt[2:length(gbt[,1]), 1]-gbt[1:length(gbt[,1])-1, 1]
			waits<-cbind(waits, gbt[1:length(gbt[,1])-1, 2])
			
			#get waits that correspond to m...n+1 and n in R&Y formulation
			subw<-(length(waits[,1])-1)
			finalwait<-subw+1
			subw<-1:subw
			
			#lambda is the coalescent branching rate for each interval
			lambda<-((waits[subw,2]*(waits[subw,2]-1))/(2*demo[node,"dmv"]*ploidy))
			#exponent is the exponential term of the exponential distribution
			exponent<-exp(-lambda*waits[subw,1])
			
			#now multiply 1/theta * exponent for all waiting times and then multiply by probability of no coalescence from the last coalescence to the end of the branch
			
			coalprob<-prod(((demo[node,"dmv"]*ploidy)^-1)*exponent)
			noncoalprob<-exp(-((2*demo[node,"dmv"]*ploidy)^-1)*waits[finalwait,2]*(waits[finalwait,2]-1)*waits[finalwait,1])
			
			if(waits[finalwait,1]==Inf){
				noncoalprob<-1
				}
			
			probability<-coalprob*noncoalprob
			
			}
		
		
		return(probability)
		}

	######now actually calculate the whole gene tree probability
		####remember to exclude species tip branches with only one allele sampled
	
	demo<-demo.maker(sptree, nspec)
	numtips<-lapply(tiplist, length)
	if(any(numtips==1)){
		demo2<-demo[-which(numtips==1),]}	
		else{demo2<-demo}
	
	lnP<-c()
	nodes<-demo2[,"node"]
	names(nodes)<-NULL
	for(i in 1:length(nodes)){
		lnP<-c(lnP, log(b.prob(sptree, gtree, demo, nodes[i])))
		}
	
	return((sum(lnP)))
	}

