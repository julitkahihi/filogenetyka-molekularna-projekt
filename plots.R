#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ape); library(phangorn) })
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop('Usage: Rscript src/plots_v2.R /path/to/package')
base <- normalizePath(args[1], mustWork=TRUE)
figdir <- file.path(base,'results','figures'); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
tabdir <- file.path(base,'tables'); dir.create(tabdir,recursive=TRUE,showWarnings=FALSE)
igf_path <- file.path(base,'data','igf1','IGF1_ML.treefile')
gen_path <- file.path(base,'data','reference','GENOME_ML.treefile')
if (!file.exists(igf_path) || !file.exists(gen_path)) stop('Need both IQ-TREE2 .treefile files')
clean <- function(x) gsub('_',' ',sub('\\|.*$','',trimws(as.character(x))),fixed=TRUE)
meta <- read.delim(file.path(base,'inputs','species.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
meta$species <- clean(meta$species)
body <- read.delim(file.path(base,'inputs','body_mass_kg.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c('species','mass_kg') %in% names(body)))
body$species <- clean(body$species)
body$mass_kg <- suppressWarnings(as.numeric(gsub(',','.',body$mass_kg,fixed=TRUE)))
if (!('mass_source' %in% names(body))) body$mass_source <- 'user_supplied'
body <- body[!duplicated(body$species,fromLast=TRUE),]

metadata <- function(t) {
  sp <- clean(t$tip.label)
  x <- meta[match(sp,meta$species),,drop=FALSE]
  if (anyNA(x$clade)) stop('Taxa missing clade: ', paste(sp[is.na(x$clade)],collapse=', '))
  x$species <- sp
  x$mass_kg <- body$mass_kg[match(sp,body$species)]
  x$mass_source <- body$mass_source[match(sp,body$species)]
  x$rodent <- x$clade %in% c('Castorimorpha','Sciuromorpha','Hystricomorpha','Myomorpha','Anomaluromorpha')
  x
}
colors <- c(Castorimorpha='#C44367',Sciuromorpha='#CF5830',Hystricomorpha='#B47705',Myomorpha='#807800',
            Lagomorpha='#0075AD',Marsupialia='#7164C9',Primates='#B044A3',Carnivora='#21844B',Cetacea='#008D73',
            Artiodactyla='#618F00',Perissodactyla='#855BBA',Chiroptera='#008C8C',Eulipotyphla='#00A3AF',
            Proboscidea='#C3418D',Cingulata='#008EA8',Scandentia='#8059A5')
all_groups <- sort(unique(meta$clade))
missing_clades <- setdiff(all_groups,names(colors))
if (length(missing_clades)) colors <- c(colors,setNames(grDevices::hcl.colors(length(missing_clades),'Dark 3'),missing_clades))

split_key <- function(t) {
  n <- length(t$tip.label)
  tiplabels <- t$tip.label
  edge <- t$edge
  desc <- function(v) {
    if (v <= n) return(tiplabels[v])
    kids <- edge[edge[,1]==v,2]
    unlist(lapply(kids,desc),use.names=FALSE)
  }
  vapply(n+seq_len(t$Nnode),function(v) {
    left <- sort(desc(v)); right <- sort(setdiff(tiplabels,left))
    if (min(length(left),length(right)) < 2) return(NA_character_)
    l <- paste(left,collapse='|'); r <- paste(right,collapse='|')
    if (length(left)<length(right)) l else if (length(right)<length(left)) r else min(l,r)
  },'')
}
boot_parse <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_real_)
  p <- strsplit(s,'/',fixed=TRUE)[[1]]
  z <- suppressWarnings(as.numeric(tail(p,1)))
  if (length(z)==1 && is.finite(z) && z>=0 && z<=100) z else NA_real_
}
root_supported <- function(t) {
  s <- if (is.null(t$node.label)) rep(NA_real_,t$Nnode) else vapply(t$node.label,boot_parse,0.0)
  spl <- split_key(t)
  valid <- !is.na(spl) & is.finite(s)
  dict <- split(s[valid],spl[valid])
  otu <- which(clean(t$tip.label)=='Monodelphis domestica')
  if (length(otu)==1) {
    tt <- ape::root(t,outgroup=t$tip.label[otu],resolve.root=TRUE)
    status <- 'Monodelphis domestica'
  } else {
    tt <- phangorn::midpoint(t)
    status <- 'MIDPOINT: outgroup unavailable'
  }
  tt <- ape::ladderize(tt)
  new <- split_key(tt)
  sup <- vapply(new,function(k) if (!is.na(k) && !is.null(dict[[k]])) dict[[k]][1] else NA_real_,0.0)
  list(tree=tt,support=sup,root=status)
}

# reorder a rooted tree to approximate a target tip order, preserving topology
reorder_to_match <- function(t, target_order) {
  target_rank <- setNames(seq_along(target_order), target_order)
  n <- length(t$tip.label)
  parents <- t$edge[,1]
  kids_all <- t$edge[,2]
  root <- setdiff(parents, kids_all)[1]
  children <- split(kids_all, parents)
  edge_index <- setNames(seq_len(nrow(t$edge)), paste(t$edge[,1], t$edge[,2], sep='>'))
  desc_cache <- vector('list', n + t$Nnode)
  desc_fun <- function(v) {
    if (!is.null(desc_cache[[v]])) return(desc_cache[[v]])
    out <- if (v <= n) t$tip.label[v] else unlist(lapply(children[[as.character(v)]], desc_fun), use.names=FALSE)
    desc_cache[[v]] <<- out
    out
  }
  reorder_node <- function(v) {
    if (v <= n) return()
    kids <- children[[as.character(v)]]
    for (k in kids) if (k > n) reorder_node(k)
    kids <- children[[as.character(v)]]
    scores <- vapply(kids, function(k) mean(target_rank[desc_fun(k)], na.rm=TRUE), 0.0)
    children[[as.character(v)]] <<- kids[order(scores)]
  }
  reorder_node(root)
  new_edge <- matrix(numeric(0), ncol=2)
  new_len <- if (!is.null(t$edge.length)) numeric(0) else NULL
  emit <- function(v) {
    if (v <= n) return()
    kids <- children[[as.character(v)]]
    for (k in kids) {
      idx <- edge_index[[paste(v,k,sep='>')]]
      new_edge <<- rbind(new_edge, t$edge[idx,,drop=FALSE])
      if (!is.null(new_len)) new_len <<- c(new_len, t$edge.length[idx])
      emit(k)
    }
  }
  emit(root)
  t$edge <- new_edge
  if (!is.null(new_len)) t$edge.length <- new_len
  t
}

input_igf <- ape::read.tree(igf_path); input_gen <- ape::read.tree(gen_path)
input_igf$tip.label <- clean(input_igf$tip.label)
input_gen$tip.label <- clean(input_gen$tip.label)
if (anyDuplicated(input_igf$tip.label) || anyDuplicated(input_gen$tip.label)) stop('Nonunique species tips after name normalization')
common <- intersect(input_igf$tip.label,input_gen$tip.label)
if (length(common) < 6) stop('Fewer than 6 common taxa between IGF1 and genome reference')
cat('Common taxa:',length(common),'\n')
prune <- function(t) { unused <- which(!t$tip.label %in% common); if (length(unused)) ape::drop.tip(t,unused) else t }
igf <- root_supported(prune(input_igf)); gen <- root_supported(prune(input_gen))
ref_tax <- metadata(gen$tree)
if (!setequal(igf$tree$tip.label,gen$tree$tip.label)) stop('Tree species sets differ after pruning')
if (!all(is.finite(igf$tree$edge.length)) || !all(is.finite(gen$tree$edge.length))) stop('Need finite branch lengths')

write.table(data.frame(species=ref_tax$species,clade=ref_tax$clade,mass_kg=ref_tax$mass_kg,mass_source=ref_tax$mass_source),
            file.path(tabdir,'plot_taxa.tsv'),quote=FALSE,sep='\t',row.names=FALSE,na='')

a <- igf$tree; b <- gen$tree
RF <- as.numeric(phangorn::RF.dist(a,b,normalize=FALSE,rooted=FALSE))
RF_norm <- as.numeric(phangorn::RF.dist(a,b,normalize=TRUE,rooted=FALSE))
d1 <- ape::cophenetic.phylo(a); d2 <- ape::cophenetic.phylo(b)
ix <- sort(common)
if (is.null(rownames(d1)) || is.null(rownames(d2))) stop('Missing species names in a patristic distance matrix')
d1 <- d1[ix,ix,drop=FALSE]; d2 <- d2[ix,ix,drop=FALSE]
mask <- upper.tri(d1)
x <- as.numeric(d1[mask]); y <- as.numeric(d2[mask])
valid <- is.finite(x)&is.finite(y)&x>0&y>0
spearman <- if (sum(valid)>3) suppressWarnings(cor(x[valid],y[valid],method='spearman')) else NA_real_
pearson <- if (sum(valid)>3) suppressWarnings(cor(x[valid],y[valid],method='pearson')) else NA_real_
pair_i <- which(mask,arr.ind=TRUE)
pairdist <- data.frame(species_1=ix[pair_i[,1]],species_2=ix[pair_i[,2]],igf1_distance=x,genomic_marker_distance=y)
med1 <- median(x[x>0 & is.finite(x)]); med2 <- median(y[y>0 & is.finite(y)])
pairdist$igf1_normalized <- pairdist$igf1_distance/med1
pairdist$genome_normalized <- pairdist$genomic_marker_distance/med2
pairdist$log2_IGF1_vs_genome <- log2(pairdist$igf1_normalized/pairdist$genome_normalized)
write.table(pairdist,file.path(tabdir,'pairwise_distances.tsv'),row.names=FALSE,quote=FALSE,sep='\t')

def_m <- function(mat,sp) median(mat[sp,setdiff(ix,sp)],na.rm=TRUE)
rel <- data.frame(species=ix,clade=meta$clade[match(ix,meta$species)],
                  IGF1_median_distance=vapply(ix,function(sp) def_m(d1,sp),0.0),
                  genome_median_distance=vapply(ix,function(sp) def_m(d2,sp),0.0))
rel$IGF1_normalized <- rel$IGF1_median_distance/med1
rel$genome_normalized <- rel$genome_median_distance/med2
rel$relative_IGF1_divergence <- rel$IGF1_normalized/rel$genome_normalized
rel$log2_relative_IGF1 <- log2(rel$relative_IGF1_divergence)
rel$mass_kg <- body$mass_kg[match(rel$species,body$species)]
write.table(rel,file.path(tabdir,'IGF1_relative_divergence.tsv'),row.names=FALSE,quote=FALSE,sep='\t',na='')
metrics <- data.frame(n_common=length(common),n_IGF1=length(input_igf$tip.label),n_genomic=length(input_gen$tip.label),
                      RF=RF,RF_normalized=RF_norm,spearman_patristic=spearman,pearson_patristic=pearson,
                      rooting_igf1=igf$root,rooting_genome=gen$root)
write.table(metrics,file.path(tabdir,'tree_comparison_metrics.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')

point_cex <- function(m) {
  x <- (log10(m)-log10(0.015))/(log10(25)-log10(0.015))
  0.7+2.1*sqrt(pmin(1,pmax(0,x)))
}
edge_colours <- function(t,meta_t) {
  n <- length(t$tip.label); edge <- t$edge
  groups <- function(node) {
    if (node<=n) return(meta_t$clade[node])
    unique(unlist(lapply(edge[edge[,1]==node,2],groups),use.names=FALSE))
  }
  vapply(edge[,2],function(k) { cl <- groups(k); if (length(cl)==1) unname(colors[cl]) else '#585858' },'')
}
clades_to_show <- intersect(c('Castorimorpha','Sciuromorpha','Hystricomorpha','Myomorpha', sort(setdiff(all_groups,c('Castorimorpha','Sciuromorpha','Hystricomorpha','Myomorpha')))), unique(ref_tax$clade))


# ---- Minimal publication-style plotting ----
draw_tree <- function(obj, label, bootstrap_min=70, cex_labels=.86) {
  t <- obj$tree
  n <- length(t$tip.label)
  mm <- metadata(t)
  tip_cols <- unname(colors[mm$clade])
  tmp <- t
  tmp$edge.length <- rep(1, nrow(tmp$edge))
  maxd <- max(ape::node.depth.edgelength(tmp)[seq_len(n)])
  ape::plot.phylo(tmp, show.tip.label=FALSE, use.edge.length=TRUE,
                  edge.color=edge_colours(t,mm), edge.width=1.6,
                  direction='rightwards', type='phylogram', no.margin=TRUE,
                  x.lim=c(-.03*maxd, 1.53*maxd))
  xy <- get('last_plot.phylo', envir=ape:::.PlotPhyloEnv)
  known <- mm$rodent & is.finite(mm$mass_kg) & mm$mass_kg>0
  sizes <- rep(.38,n)
  sizes[known] <- point_cex(mm$mass_kg[known])
  sizes[mm$rodent & !known] <- .65
  pch <- ifelse(mm$rodent & !known, 1, 21)
  points(xy$xx[seq_len(n)],xy$yy[seq_len(n)],pch=pch,
         cex=sizes,col=tip_cols,bg=tip_cols,lwd=1)
  text(xy$xx[seq_len(n)]+maxd*.018,xy$yy[seq_len(n)],
       labels=clean(t$tip.label),adj=0,col=tip_cols,font=3,cex=cex_labels)
  nodes <- n+seq_len(t$Nnode)
  show <- which(is.finite(obj$support) & obj$support>=bootstrap_min)
  if (length(show))
    text(xy$xx[nodes[show]]+maxd*.008,xy$yy[nodes[show]]+.28,
         labels=round(obj$support[show]),col='#323232',cex=.53,font=2)
  title(main=label,adj=0,line=.2,cex.main=1.12)
}

legend_page <- function(ncolumns=8) {
  par(mar=c(0, .5, 0, .5), xpd=NA)
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
  nc <- min(ncolumns,length(clades_to_show))
  for (i in seq_along(clades_to_show)) {
    x <- .008 + ((i-1) %% nc)*(.99/nc)
    y <- .76 - floor((i-1)/nc)*.30
    points(x,y,pch=15,col=colors[clades_to_show[i]],cex=.95)
    text(x+.011,y,clades_to_show[i],adj=0,cex=.66)
  }
  weights <- c(.02,.1,1,10)
  for (i in seq_along(weights)) {
    x <- .10+(i-1)*.12
    points(x,.12,pch=21,col='#777777',bg='#777777',cex=point_cex(weights[i]))
    text(x+.015,.12,paste0(weights[i],' kg'),adj=0,cex=.69)
  }
  points(.60,.12,pch=1,cex=.70,col='#777777')
  text(.617,.12,'Unknown mass',adj=0,cex=.69)
}

make_plot <- function(filename,w,h,fun) {
  target <- file.path(figdir,filename)
  if (capabilities('cairo')) png(target,width=w,height=h,res=240,type='cairo-png',bg='white')
  else png(target,width=w,height=h,res=240,bg='white')
  tryCatch(fun(),finally=dev.off())
  pdf_path <- sub('\\.png$','.pdf',target)
  pdf(pdf_path,width=w/240,height=h/240,onefile=TRUE,useDingbats=FALSE)
  tryCatch(fun(),finally=dev.off())
  cat('Figure:',target,'\n')
}

make_plot('01_IGF1_tree.png',4200,3050,function() {
  layout(matrix(c(1,2),ncol=1),heights=c(9,1.15))
  par(mar=c(.25,.6,1.2,.4),xpd=FALSE)
  draw_tree(igf,'IGF1')
  legend_page()
})
make_plot('02_GENOME_reference_tree.png',4200,3050,function() {
  layout(matrix(c(1,2),ncol=1),heights=c(9,1.15))
  par(mar=c(.25,.6,1.2,.4),xpd=FALSE)
  draw_tree(gen,'BUSCO reference')
  legend_page()
})
make_plot('03_IGF1_vs_GENOME_side_by_side.png',6300,3050,function() {
  layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(9,1.15))
  par(mar=c(.25,.45,1.2,.3),xpd=FALSE)
  draw_tree(igf,'IGF1',cex_labels=.77)
  par(mar=c(.25,.45,1.2,.3),xpd=FALSE)
  draw_tree(gen,'BUSCO reference',cex_labels=.77)
  legend_page()
})

# Coordinates from directed rooted tree edges.
# IMPORTANT: every leaf must receive its own y from a root-to-tip traversal.
coords <- function(t) {
  n <- length(t$tip.label)
  edge <- t$edge
  roots <- setdiff(unique(edge[,1]), unique(edge[,2]))
  if (length(roots)!=1L) stop('Tanglegram: expected exactly one root')
  root <- roots[1]
  children <- split(edge[,2], as.character(edge[,1]))
  yy <- rep(NA_real_,n+t$Nnode)
  counter <- 0L
  assign_y <- function(v) {
    if (v<=n) {
      counter <<- counter+1L
      yy[v] <<- counter
    } else {
      cc <- children[[as.character(v)]]
      if (!length(cc)) stop('Tanglegram: internal node has no descendants')
      for (k in cc) assign_y(k)
      yy[v] <<- mean(range(yy[cc]))
    }
  }
  assign_y(root)  # FIX: absent from v2; without this every leaf is drawn at y=0
  if (counter!=n || any(!is.finite(yy)) || length(unique(yy[seq_len(n)]))!=n)
    stop('Tanglegram: failed to calculate distinct y positions for all tips')
  depth <- numeric(n+t$Nnode)
  assign_depth <- function(v) {
    cc <- children[[as.character(v)]]
    if (length(cc)) for (k in cc) {
      depth[k] <<- depth[v]+1
      if (k>n) assign_depth(k)
    }
  }
  assign_depth(root)
  list(x=depth,y=(yy-.5)/n,edge=edge,root=root,
       order=t$tip.label[order(yy[seq_len(n)])])
}

# Tree rotations change tip display order only, not branching relationships.
# Uses the existing reorder_to_match() defined above.
optimize_tanglegram <- function(left,right,rounds=6L) {
  best_left <- left; best_right <- right
  c0 <- function(l,r) {
    p <- match(coords(l)$order,coords(r)$order)
    n <- length(p)
    sum(outer(p,p,'>')[upper.tri(matrix(FALSE,n,n))])
  }
  best <- c0(left,right)
  for (iter in seq_len(rounds)) {
    right <- reorder_to_match(right,coords(left)$order)
    left  <- reorder_to_match(left,coords(right)$order)
    new <- c0(left,right)
    if (new < best) {
      best <- new; best_left <- left; best_right <- right
    }
  }
  cat('Tanglegram crossings after rotations:',best,'\n')
  list(left=best_left,right=best_right,crossings=best)
}

make_plot('04_IGF1_vs_GENOME_tanglegram.png',5700,3100,function() {
  par(mar=c(.25,.35,1,.35),xpd=FALSE)
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
  opt <- optimize_tanglegram(a,b)
  A <- coords(opt$left); B <- coords(opt$right)
  n <- length(opt$left$tip.label)
  left_names <- opt$left$tip.label
  right_names <- opt$right$tip.label
  idx <- match(left_names,right_names)
  if (anyNA(idx)) stop('Tanglegram: some species could not be matched')
  tax <- setNames(ref_tax$clade,ref_tax$species)

  # Left tree points right and right tree points left.
  X1 <- .04 + .395*A$x/max(A$x)
  X2 <- .96 - .395*B$x/max(B$x)
  Y1 <- .04 + .92*A$y
  Y2 <- .04 + .92*B$y
  draw_edges <- function(P,XX,YY) {
    e <- P$edge
    for (j in seq_len(nrow(e))) {
      parnode <- e[j,1]; child <- e[j,2]
      segments(XX[parnode],YY[child],XX[child],YY[child],col='#656565',lwd=.90)
    }
    for (p in unique(e[,1])) {
      cc <- e[e[,1]==p,2]
      segments(XX[p],min(YY[cc]),XX[p],max(YY[cc]),col='#656565',lwd=.90)
    }
  }
  draw_edges(A,X1,Y1)
  draw_edges(B,X2,Y2)

  # Every connector is the same species in the two trees.
  link_col <- grDevices::adjustcolor(unname(colors[tax[left_names]]),alpha.f=.49)
  segments(X1[seq_len(n)],Y1[seq_len(n)],X2[idx],Y2[idx],col=link_col,lwd=1.10)
  points(X1[seq_len(n)],Y1[seq_len(n)],pch=16,cex=.59,col=colors[tax[left_names]])
  points(X2[seq_len(n)],Y2[seq_len(n)],pch=16,cex=.59,col=colors[tax[right_names]])

  # Only two panel labels and the topology metric; no redundant footer.
  text(.04,.987,'IGF1',adj=0,cex=1.17,font=2)
  text(.96,.987,'BUSCO reference',adj=1,cex=1.17,font=2)
  text(.50,.987,sprintf('RF = %.3f',RF_norm),adj=.5,cex=.95,col='#555555')
})

make_plot('05_IGF1_relative_divergence.png',3500,2250,function() {
  rr <- rel[is.finite(rel$log2_relative_IGF1),]
  rr <- rr[order(rr$log2_relative_IGF1),]
  if (nrow(rr)>24) rr <- rr[c(seq_len(12),seq.int(nrow(rr)-11,nrow(rr))),]
  xpos <- rr$log2_relative_IGF1
  yy <- seq_len(nrow(rr))
  rx <- range(c(xpos,0),finite=TRUE)
  span <- diff(rx); if (span<=0) span <- 1
  par(mar=c(3.7,10.3,1.5,.8),xpd=FALSE)
  plot(xpos,yy,type='n',yaxt='n',ylab='',xlim=rx+c(-.08,.08)*span,
       ylim=c(.4,nrow(rr)+.6),xlab='log2 relative IGF1 divergence',bty='l')
  abline(v=0,col='#888888',lty=2,lwd=1)
  segments(0,yy,xpos,yy,col='#B6B6B6',lwd=.95)
  points(xpos,yy,pch=16,cex=.93,col=colors[rr$clade])
  axis(2,at=yy,labels=rr$species,las=2,tick=FALSE,cex.axis=.79,line=-.5)
  title(main='IGF1 relative divergence',adj=0,cex.main=1.10)
})

make_plot('06_pairwise_distances.png',3400,2500,function() {
  par(mar=c(4.3,4.5,1.6,.8))
  u <- pairdist$genome_normalized; v <- pairdist$igf1_normalized
  good <- is.finite(u)&is.finite(v)&u>=0&v>=0
  lim <- range(c(u[good],v[good]))
  smoothScatter(u[good],v[good],nbin=250,nrpoints=0,
                colramp=colorRampPalette(c('#FFFFFF','#E4EDF4','#B2C9DD','#648CAB','#264864')),
                xlab='BUSCO normalized distance',
                ylab='IGF1 normalized distance',
                xlim=lim,ylim=lim)
  abline(0,1,lty=2,col='#BD6C61',lwd=1.45)
  title(main=sprintf('Pairwise divergence  |  Spearman rho = %.2f',spearman),cex.main=1.10)
})
cat('DONE. FIGURES:',figdir,'\n')
