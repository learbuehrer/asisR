##############################################################################
# Functions for the MCMC ABN analysis
##############################################################################

##############################################################################
# 1. Find optimal no. of parent nodes
##############################################################################
findMaxParents <- function(METHOD= "mle",  SCORE ="bic",              ## method. cache parameters
                           dat, dists, banmat=NULL, retainmat=NULL,   ## data cache parameters
                           n.cores=11, FILENAMEbase, PRIOR=1,         ## computational parameters
                           group.var=NULL)         
{# initialize function, save starting time and print start:
  cat(paste("\nRun the exact search across incremental parent limits with method:",METHOD))
  starttime <- Sys.time()
  
  # initialize empty data frame for net scores:
  if(is.null(group.var)){
    novars <- ncol(dat)}  else {
      novars <- ncol(dat)-1}
  
  # find the minimal number of max parents, according to the retain matrix:
  if(!is.null(retainmat)){startvalue <- max(1,unname(rowSums(retainmat)[rowSums(retainmat) > 1]))}  else
  {startvalue = 1}
  
  tmpscores <- vector(length = novars)
  net.scores <-    data.frame(npar = NULL,
                              scoretype = NULL,
                              scorevalue = NULL)
  
  # run abn per number of parents possible:
  clust <-
    parallel::makeCluster(n.cores, outfile = paste0(FILENAMEbase, "maxparents.log"))
  doParallel::registerDoParallel(cl = clust)
  net.scores <- foreach(
    i = startvalue:novars,
    .combine = 'rbind',
    .packages = c("abn", "mcmcabn"),
    .inorder = TRUE) %dopar% {
      max.par <- i      
      mycache <- buildScoreCache(
        data.df = as.data.frame(dat),
        data.dists = dists,
        dag.banned = banmat,
        dag.retained = retainmat,
        max.parents = max.par,
        method = METHOD,
        group.var = group.var)
      
      dag.mP <- mostProbable(score.cache = mycache,
                             score = SCORE, prior.choice=PRIOR)
      fabn.mP <- fitAbn(object = dag.mP,
                        method = METHOD, prior.choice=PRIOR)
      
      return(list(i=list(dag.mP, c(i, SCORE, fabn.mP[[SCORE]]))))}
  stopCluster(clust)
  endtime <- Sys.time()
  cat(paste("\n****************************\nEnd MCMC ABN. maxparent searching, Time used [h]:", 
            round(difftime(endtime, starttime, units = "hours"), 4)))
  
  
  ## format net.scores, such that all information is in a nice structure:
  net.scores.dags.prelim <- unlist(net.scores, recursive = F)
  net.scores.dags <- net.scores.dags.prelim[seq(1,length(net.scores.dags.prelim), 2)]
  net.scores <- net.scores.dags.prelim[seq(2,length(net.scores.dags.prelim), 2)]
  
  net.scores <- data.frame(npar = as.integer(unlist(net.scores)[seq(1,length(unlist(net.scores)), 3)]),
                           scoretype = unlist(net.scores)[seq(2,length(unlist(net.scores)), 3)],
                           scorevalue = as.numeric(unlist(net.scores)[seq(3,length(unlist(net.scores)), 3)]))
  rownames(net.scores) <- NULL
  
  
  # find the first number of parents with the lowest aic/bic score
  x <- net.scores %>%
    filter(scoretype == SCORE)
  max.par <- x$npar[which(x$scorevalue == min(x$scorevalue))][1]
  cat(paste("\n****************************\nmax. parents:", max.par))
  
  return(list("max.par"= max.par, 
              "net.scores" = net.scores,
              "net.scores.dags" = net.scores.dags))}


##############################################################################
# 2. Looking at the maxparents:
##############################################################################
# Plof of Score over the number of max parents:
ScorePlotMaxParent <- function(net.scores,FILENAMEbase){
  df <- net.scores %>%
    mutate(
      scoretype = as.factor(scoretype),
      npar = as.factor(npar),
      scorevalue = as.numeric(scorevalue)) %>%
    
    group_by(scoretype) %>%
    reframe(
      scorevalue.norm = (1 - abs(1 - scorevalue / max(scorevalue))) * 100,
      .groups = "keep",
      npar = npar,
      scorevalue = scorevalue)
  
  plt.rel <- ggplot(df, aes(x=npar, color = scoretype, y = scorevalue.norm)) +
    geom_point(aes(shape = scoretype, size = 1), show.legend = F)+
    labs(title =  "Relative network score per no. parent nodes",
         x="number of parent nodes",
         y="network score [%]") +
    scale_x_discrete(limits = unique(df$npar))+
    theme_bw()+
    labs(caption = paste0("Score type = ", net.scores$scoretype[1]))
  
  plt.abs <- ggplot(df, aes(color = scoretype, y = scorevalue, group = scoretype)) +
    geom_point(aes(x=npar, shape = scoretype, size = 1), show.legend = F)+
    labs(title = "Absolute network score per no. parent nodes",
         x="number of parent nodes",
         y="network score") +
    scale_x_discrete(limits = unique(df$npar))+
    theme_bw()+
    labs(caption = paste0("Score type = ", net.scores$scoretype[1]))
  
  plt.comb <- cowplot::plot_grid(plt.rel, plt.abs, labels = "AUTO", ncol = 1)
  plt.comb
  
  ggsave(plt.comb, filename = paste0(FILENAMEbase, "ScorePlotMaxParent.png"),
         width = 7,height = 9, dpi=600)}


# DAGs for the max parents:
DAGMaxParent <- function(net.scores.dags, net.scores, par1, par2, retainmat = NULL,
                         FILENAMEbase, filename, width = 12, height = 9, dpi = 600) {
  width_px <- width * dpi
  height_px <- height * dpi
  
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  par(mfrow = c(par1, par2))
  
  if (!is.null(retainmat)) {
    parents <- max(1, unname(max(rowSums(retainmat))))
  } else {
    parents <- 1
  }
  
  for (idx in seq_along(net.scores.dags)) {
    plot(net.scores.dags[[idx]])
    box("figure", col = "gray", lwd = 1)
    title(sub = paste0("DAG with max.par = ", parents, 
                       ", netscore = ", round(net.scores[idx, 3], 3)))
    parents <- parents + 1
  }
  
  dev.off()
}


##############################################################################
# 3. Fit MCMCM abn() with max.par
##############################################################################
mcmcmABN <- function(METHOD= "mle",  SCORE ="bic",              ## method. cache parameters
                     dat, dists, banmat=NULL, retainmat=NULL,   ## data cache parameters
                     n.cores=11, FILENAMEbase,  FILENAME,       ## computational parameters
                     max.par = length(dat)-1,                   ## optimal number of max.par
                     MCMC.SEEDS = c("11", "12", "13", "14"),    ## set seeds for evaluation
                     MCMC.SCHEME = c(500,0,0),                  ## Meaning: c(RETURN.DAGS, THINNING, BURNIN.LEN)
                     MCMC.PRIOR = 2,                            ## 2 = Koivisto Prior, 1 = uniformative prior
                     PROB.REV = 0.05,                           ## prob of selecting a new edge reversal
                     PROB.MBR = 0.05,                           ## prob of selecting a Markov blanket resampling move
                     fileaddendum="",                           ## file name addendum for cross-validation
                     group.var=NULL                             ## variable for moxed effects
)                          

# REV and MBR are efficient in producing high scoring structure but 
# computationally costly compared to classical MCMC jumps.

{## fit abn with the given number of max.par:
  cat(paste("\nfit ABN with max.par:"))
  starttime <- Sys.time()
  mycache.maxpar <- buildScoreCache(
    data.df = as.data.frame(dat),
    data.dists = dists,
    dag.banned = banmat,
    dag.retained = retainmat,
    max.parents = max.par,
    group.var = group.var,
    method = METHOD)
  
  dag.maxpar <- mostProbable(score.cache = mycache.maxpar,
                             score = SCORE)
  fabn.maxpar <- fitAbn(object = dag.maxpar,
                        method = METHOD)
  endtime <- Sys.time()
  cat(paste("\n***********************\nEnd ABN with max.par. Time used [h]:", 
            round(difftime(endtime, starttime, units = "hours"), 4)))
  
  ### Bootstrapping:
  starttime <- Sys.time()
  cat(paste("\n***********************\nnumber of seeds:", length(MCMC.SEEDS)))
  cat(paste("\nStart MCMC ABN with max.par:", max.par))
  cat(paste("\nMCMC scheme (return dags, thinning, burnin length):", paste(MCMC.SCHEME, collapse = ", ")))
  cat(paste("\nStarting Time:", starttime))
  clust <-
    parallel::makeCluster(length(MCMC.SEEDS),
                          outfile = paste0(FILENAMEbase, "mcmcabnmaxpar_",fileaddendum, ".log"))
  doParallel::registerDoParallel(cl = clust)
  mcmc.out.list <- foreach(
    SEED = MCMC.SEEDS,
    .packages = c("abn", "mcmcabn"),
    .inorder = TRUE
  ) %dopar% {
    mcmcabn(
      score.cache = mycache.maxpar,
      score = SCORE,
      data.dists = dists,
      max.parents = max.par,
      mcmc.scheme = MCMC.SCHEME,
      seed = SEED,
      verbose = FALSE,
      start.dag = dag.maxpar$dag,
      prob.rev = PROB.REV,
      prob.mbr = PROB.MBR,
      prior.choice = MCMC.PRIOR) }
  stopCluster(clust)
  endtime <- Sys.time()
  cat(paste("\n***********************\n End MCMC ABN. Time used [h]:", 
            round(difftime(endtime, starttime, units = "hours"), 4)))
  
  ### Save raw data
  save(list = ls(),
       file = paste0(FILENAMEbase, FILENAME, "_final.RData"))
  cat("\n***********************\nFinal data written to:",
      paste0(FILENAMEbase, FILENAME, "_final.RData"))
  
  return(list("mcmc.out.list" = mcmc.out.list,
              "dag.maxpar" = dag.maxpar,
              "fabn.maxpar" = fabn.maxpar))}




##############################################################################
# 4. Postpreprocessing
##############################################################################
postBURNin <- function(mcmc.out.list, burnin.length){
  for (d in 1:length(mcmc.out.list)){
    mcmc.out.list[[d]]$burnin <- burnin.length
    mcmc.out.list[[d]]$dags <- mcmc.out.list[[d]]$dags[,,-(1:burnin.length)]
    mcmc.out.list[[d]]$scores <- mcmc.out.list[[d]]$scores[-(1:burnin.length)]
    mcmc.out.list[[d]]$alpha <- mcmc.out.list[[d]]$alpha[-(1:burnin.length)]
    mcmc.out.list[[d]]$method <- mcmc.out.list[[d]]$method[-(1:burnin.length)]
    mcmc.out.list[[d]]$rejection <- mcmc.out.list[[d]]$rejection[-(1:burnin.length)]
    mcmc.out.list[[d]]$heating <- mcmc.out.list[[d]]$heating[-(1:burnin.length)]
  }
  return(mcmc.out.list)
}


postTHINN <- function(mcmc.out.list, thinningsteps){
  thin <- c(TRUE, rep(FALSE, thinningsteps-1)) # see recycle rules
  mcmc.out.list.thinned <- mcmc.out.list
  
  for (d in 1:length(mcmc.out.list.thinned)){
    if(thinningsteps > mcmc.out.list.thinned[[d]]$iterations){
      stop(paste0("More thinning steps (", thinningsteps, ") than mcmc iterations (", mcmc.out.list.thinned[[d]]$iterations, ")."))
    }
    mcmc.out.list.thinned[[d]]$thinning <- thinningsteps
    mcmc.out.list.thinned[[d]]$dags <- mcmc.out.list.thinned[[d]]$dags[,,thin]
    mcmc.out.list.thinned[[d]]$scores <- mcmc.out.list.thinned[[d]]$scores[thin]
    mcmc.out.list.thinned[[d]]$alpha <- mcmc.out.list.thinned[[d]]$alpha[thin]
    mcmc.out.list.thinned[[d]]$method <- mcmc.out.list.thinned[[d]]$method[thin]
    mcmc.out.list.thinned[[d]]$rejection <- mcmc.out.list.thinned[[d]]$rejection[thin]
  }
  return(mcmc.out.list.thinned)}



##############################################################################
# 5. MCMC Quality Assessment: 
##############################################################################
## Gelman Diagnostics:
MCMCQualityAssessment <- function(mcmcoutputburnedthinned,mcmcoutput,
                                  FILENAMEbase,filename,  width = 6, height = 6, dpi = 600,
                                  ymax= 1.5){
  cat("QUALITY ASSESSMENT:\n*********************\nGelman: of higher than 1.1 or 1.2, run chain longer to improve convergence\n")
  
  ## Gelman:
  print(gelman.diag(x = mcmcoutputburnedthinned,autoburnin = F))
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  gelman.plot(mcmcoutput, autoburnin = F, main="Gelman Plot", ylim=c(1,ymax))
  dev.off()
  
  ## Raftery:
  cat("\n*********************\nRaftery: calculate no. of iterations and no. of burn-ins to satisfy specified condition\n*********************")
  print(raftery.diag((mcmcoutputburnedthinned)))
  
  ##Heidelberg and Welch:
  cat("\n*********************\nHeidelberg and Welch:  test H0: The Markov Chain is from a stationary distribution. If not passed, chain must run longer
      \n*********************")
  for (chain in 1:length(mcmcoutputburnedthinned)){
    cat("\n------------------------\n")
    print(paste("Chain no: ", chain))
    print(heidel.diag(mcmcoutputburnedthinned[[chain]]))}
}


## Trace Plot
TracePlot <- function(modeloutput, mcouthin1, mcoutthin2, mcoutthin3, mcoutthin4,
                      MCMC.SEEDS, dat, dists, METHOD, filename, FILENAMEbase){
  fabn <-fitAbn(dag = modeloutput$dag.maxpar$dag,
                data.df = dat,
                data.dists = dists, 
                method = METHOD) 
  
  max.score <- -fabn$bic
  
  
  dta <- data.frame(mcouthin1[2:4],
                    mcoutthin2[2:4],
                    mcoutthin3[2:4],
                    mcoutthin4[2:4]) 
  dta <- dta[,seq(1,3*length(MCMC.SEEDS),3)]
  names(dta) <- c("Run1","Run2", "Run3", "Run4")
  dta$X <- (1:length(dta$Run1))
  dta <- reshape2::melt(dta, "X", value.name = "scores")
  dta$cummax[1] <- dta$scores[1]
  for (i in 2:length(dta$scores)) {
    if (dta$scores[i] > dta$cummax[i - 1]) {
      dta$cummax[i] <- dta$scores[i]
    } else {dta$cummax[i] <- dta$cummax[i - 1]}}
  
  original_plot <- ggplot(data = dta, aes_string(x = "X", y="scores", color = "variable")) +
    geom_line(alpha = 0.8,lwd=1.1) +
    geom_hline(yintercept = max.score,linetype = "dashed", color = "red", alpha = 1) +
    geom_text(aes(25, max.score, label = round(max.score,digits = 2), vjust = -0.5), color = "red", check_overlap = TRUE) +
    labs(x = "DAG index", y = "DAG scores", colour = "MCMC:") +
    ggpubr::theme_pubr()+
    ylim(min(dta$scores),max(dta$scores)) 
  print(original_plot)
  
  # Plot
  y_density <- cowplot::axis_canvas(original_plot, axis = "y", coord_flip = TRUE) +
    geom_density(data = dta, aes_string(x = "scores",fill = "factor(variable)"), alpha = 0.5) +
    coord_flip()
  y_density
  
  cummax_plt <- ggplot(data = dta, aes_string(x = "X", y="cummax", color = "variable")) +
    geom_line(alpha = 0.8,lwd=1.1, inherit.aes = T) +
    geom_point(aes_string(color = "variable"), inherit.aes = T)+
    geom_hline(yintercept = max.score,linetype = "dashed", color = "red", alpha = 1) +
    geom_text(aes(25, max.score, label = round(max.score,digits = 2), vjust = -0.5), color = "red", check_overlap = TRUE) +
    labs(x = "DAG index", y = "DAG scores", colour = "MCMC:") +
    ggpubr::theme_pubr()+
    ylim(min(dta$scores),max(dta$scores)) 
  cummax_plt
  
  # create the combined plot
  combined_plot <- cowplot::ggdraw(cowplot::insert_yaxis_grob(plot = original_plot, grob = y_density, position = "right"))
  combined_cummax_plot <- cowplot::ggdraw(cowplot::insert_yaxis_grob(plot = original_plot, grob = cummax_plt, position = "right"))
  ggsave(paste0(FILENAMEbase, filename, ".png"),
         plot = cowplot::ggdraw(cowplot::insert_yaxis_grob(plot = original_plot, grob = y_density, position = "right")),
         width = 9,height = 7, dpi = 600)
  dev.off()}

##############################################################################
# 6. Assessment of the DAG:
##############################################################################
plotDAGthresholds <- function(threshold, dat, dists, dag,  
                              FILENAMEbase,filename, width = 12, height = 9, dpi = 600){
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  
  if((length(threshold))==1){
    par(mfrow=c(1,2))
    dag_threshold <- dag
    dag_threshold[dag_threshold > threshold] <- 1
    dag_threshold[dag_threshold <= threshold] <- 0
    
    edgestren <- dag
    edgestren[edgestren <= threshold] <- 0
    
    cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                         data.df = dat,
                                         data.dists = dists,
                                         digits = 2,
                                         edge.strength = edgestren,
                                         plot = T, 
                                         main=threshold)
    title(paste0("Threshold = ",threshold))
    
    # control with bnlearn
    dagBN = bnlearn::empty.graph(names(dat))
    bnlearn::amat(dagBN) = t(dag_threshold)
    bnlearn::graphviz.plot(dagBN,
                           shape = "rectangle",
                           main = "MCMC ABN DAG")
  }
  if((length(threshold)>1)){
    par(mfrow= c(2, ceiling(length(threshold)/2)))
    for(i in 1:length(threshold)){
      dag_threshold <- dag
      dag_threshold[dag_threshold > threshold[i]] <- 1
      dag_threshold[dag_threshold <= threshold[i]] <- 0
      
      edgestren <- round(dag, 2)
      edgestren[edgestren <= threshold[i]] <- 0
      
      cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                           data.df = dat,
                                           data.dists = dists,
                                           digits = 2,
                                           edge.strength = edgestren,
                                           plot = T, 
                                           main=threshold[i])
      
      title(sub = paste0("Threshold = ",threshold[i], ", n = ", nrow(dat)))}}
  dev.off()}

# Arc Strength:
arc.stren.threshold <- function(strengthdag = dag, method = "l1",
                                width = 12, height = 8, dpi = 600) {
  # Handle empty input
  if (length(strengthdag) == 0)
    return(0)
  
  # Compute ECDF and knots
  e <- ecdf(strengthdag)
  u <- knots(e)
  
  # Define norm function based on method
  norm <- switch(method,
                 "l1" = function(p) {
                   sum(diff(unique(c(0, u, 1))) * abs(e(unique(c(0, u[u < 1]))) - p))
                 },
                 "l2" = function(p) {
                   sum(diff(unique(c(0, u, 1))) * (e(unique(c(0, u[u < 1]))) - p)^2)
                 },
                 stop("Unsupported method. Use 'l1' or 'l2'."))
  
  # Optimize norm
  p0 <- optimize(f = norm, interval = c(0, 1))$minimum
  
  # Check boundaries
  if (norm(1) < norm(p0)) p0 <- 1
  if (norm(0) < norm(p0)) p0 <- 0
  
  # Compute threshold
  arc.stren.sign.threshold <- round(quantile(strengthdag, p0, type = 1, names = FALSE), 3)
  
  return(arc.stren.sign.threshold)}


# EDCF-plot
ecdf_plot_ggplot <- function(strengthdag = dag,
                             arc.stren.sign.threshold = arc.stren.sign.threshold,
                             filename, FILENAMEbase,
                             width = 12, height = 8, dpi = 600) {
  # Handle empty input
  if (length(strengthdag) == 0)
    return(0)
  
  # Compute ECDF
  e <- ecdf(strengthdag)
  x_vals <- sort(unique(strengthdag))
  y_vals <- e(x_vals)
  df <- data.frame(x = x_vals, y = y_vals)
  
  # Add (0, 0) manually
  df <- data.frame(x = c(0, x_vals), y = c(0, y_vals))
  
  
  # Create ggplot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_hline(yintercept = 1, color = "black", linewidth = 0.5, lty = 2) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, lty = 2) +
    geom_step(color = "black") +
    geom_point()+
    geom_vline(xintercept = arc.stren.sign.threshold, linetype = "dashed", color = "darkred") +
    geom_vline(xintercept = 0.5, color = "gray", linewidth = 0.5) +
    labs(title = paste0("ECDF(DAG), threshold = ", arc.stren.sign.threshold),
         x = "Arc Strength", y = "Fn(x)") +
    theme_classic() +
    xlim(c(0, 1)) +
    ylim(c(0, 1))
  
  
  # Save plot
  ggsave(filename = paste0(FILENAMEbase, filename, ".png"),
         plot = p, width = width, height = height, dpi = dpi)
  
  return(invisible(p))}



# Threshold map:
ThresholdHeatmap <- function(dag, threshold,banmat,
                             filename, FILENAMEbase, width=12, height=8, dpi=600){
  ## banmat
  rownames(dag) <- colnames(dag) <- rownames(banmat)
  df <- melt(dag)
  
  a <- ggplot(df, aes(x = Var2, y = Var1,  fill=value)) +
    geom_tile(show.legend = T) +
    scale_fill_gradient(low="white", high="darkblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(caption= paste0("arcstrength = ", threshold))+
    ylim(rev(levels(df$Var2)))+
    scale_x_discrete(
      expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
      position = "top")+
    xlab("parent")+
    ylab("child")
  
  df$value[df$value<threshold] <- 0
  b <- ggplot(df, aes(x = Var2, y = Var1,  fill=value)) +
    geom_tile(show.legend = T) +
    scale_fill_gradient(low="white", high="darkblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(caption= paste0("arcstrength for threshold = ", threshold))+
    ylim(rev(levels(df$Var2)))+
    scale_x_discrete(
      expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
      position = "top")+
    xlab("parent")+
    ylab("child")
  
  ggarrange(a,b, nrow=1, ncol=2, labels = c("A", "B")) 
  
  ggsave(plot = ggarrange(a,b, nrow=1, ncol=2, labels = c("A", "B")),
         filename = paste0(FILENAMEbase, filename, ".png"),
         width = width,height = height, dpi = dpi)}


# Adjacency matrix
AdjacencyMatrix <- function(banmat, dag, threshold,
                            filename, FILENAMEbase, width=12, height=8, dpi=600){
  ## banmat
  df <- melt(banmat)
  
  a <- ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "palegreen3", high = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(caption= "Banmatrix")+
    ylim(rev(levels(df$Var2)))+
    scale_x_discrete(
      expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
      position = "top")+
    xlab("parent")+
    ylab("child")
  
  dag_threshold <- dag
  dag_threshold[dag_threshold > threshold] <- 1
  dag_threshold[dag_threshold <= threshold] <- 0
  rownames(dag_threshold) <- colnames(dag_threshold) <- rownames(banmat)
  
  dfdag <- melt(dag_threshold)
  b <- ggplot(dfdag, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(caption= "DAG fitted")+
    ylim(rev(levels(df$Var2)))+
    scale_x_discrete(
      expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
      position = "top")+
    xlab("parent")+
    ylab("child")
  
  ggarrange(a,b, nrow=1, ncol=2, labels = c("A", "B")) 
  
  ggsave(plot = ggarrange(a,b, nrow=1, ncol=2, labels = c("A", "B")),
         filename = paste0(FILENAMEbase, filename, ".png"),
         width = width,height = height, dpi = dpi)}



AdjacencyMatrixComparision <- function(banmat, dag_train, dag_test, 
                                       filename, FILENAMEbase, width=12, height=8, dpi=600){
  
  # get the variable names:
  df <- melt(banmat)
  
  # matrix for dag train: 
  dfdag <- melt(dag_train)
  a <- ggplot(dfdag, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(caption= "DAG train")+
    ylim(rev(levels(df$Var2)))+
    scale_x_discrete(
      expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
      position = "top")+
    xlab("parent")+
    ylab("child")
  
  
  # matrix for the test data:
  dfdag <- melt(dag_test)
  b <- ggplot(dfdag, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(caption= "DAG test")+
    ylim(rev(levels(df$Var2)))+
    scale_x_discrete(
      expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
      position = "top")+
    xlab("parent")+
    ylab("child")
  
  ggarrange(a,b, nrow=1, ncol=2, labels = c("A", "B")) 
  
  ggsave(plot = ggarrange(a,b, nrow=1, ncol=2, labels = c("A", "B")),
         filename = paste0(FILENAMEbase, filename, ".png"),
         width = width,height = height, dpi = dpi)}


DAGsabn <- function(seeds = c("0815", "1234"),dat, dists, banmat= NULL, retainmat= NULL,
                    max.parents=length(dat)-1, group.var = NULL,
                    filename, FILENAMEbase, width=12, height=8, dpi=600, prior=1, par1, par2, SCORE="bic"){
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  par(mfrow= c(par1,par2))
  
  
  for (i in 1:length(seeds)){
    chacheMLE <- buildScoreCache(data.df = dat,
                                 data.dists = dists,
                                 method = "mle",
                                 dag.banned = banmat,
                                 dag.retained = retainmat,
                                 centre = TRUE,
                                 max.parents = max.parents,
                                 group.var = group.var)
    
    dag <- mostProbable(score.cache = chacheMLE, prior.choice = prior,
                        score=SCORE)
    plot(dag)
    box("figure", col="gray", lwd = 1)
    title(sub = paste0("seed = ", seeds[i]))}
  
  
  dev.off()}


abnMCMCabnComparision <- function(dat, dists, dag,  threshold,prior=1, banmat= NULL,
                                  retainmat = NULL,
                                  FILENAMEbase,filename, width = 12, height = 9, dpi = 600,
                                  max.parents, seed= "1234",
                                  group.var=NULL){
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  par(mfrow=c(1,2))
  dag_threshold <- dag
  dag_threshold[dag_threshold > threshold] <- 1
  dag_threshold[dag_threshold <= threshold] <- 0
  edgestren <- dag
  edgestren[edgestren <= threshold] <- 0
  cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                       data.df = dat,
                                       data.dists = dists,
                                       digits = 2,
                                       edge.strength = edgestren,
                                       plot = T, 
                                       main=threshold)
  title(sub = paste0("mcmcabn(), threshold = ",threshold, ", n = ", nrow(dat)))
  set.seed(seed)
  chacheMLE1 <- buildScoreCache(data.df = dat,
                                data.dists = dists,
                                method = "mle",
                                dag.banned = banmat,
                                dag.retained = retainmat,
                                centre = TRUE,
                                max.parents = max.parents,
                                group.var = group.var)
  
  dagabn <- mostProbable(score.cache = chacheMLE1, prior.choice = prior)
  plot(dagabn)
  title(sub = paste0("abn(), seed = ", seed, ", n = ", nrow(dat)))
  
  # Close the device
  dev.off()}



###############################################################################
# 7. Functions for imputation of mcmcabn:
###############################################################################
run_MCMC_ABN_all_imputations <- function(dat_i_list, runname, dists, banmatASIS, retainmatASIS,
                                         modelmaxparents, MCMC.SEEDS, MCMC.SCHEME,
                                         burnin.length = 1000, thinningsteps = 2,
                                         keep, baseline = 0) {
  # Create results directory
  outdir <- file.path("results", runname)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  model_list <- list()
  dag_list <- list()
  
  for (i in seq_along(dat_i_list)) {
    imp_index <- i + baseline
    cat("\n▶ Running MCMC for imputation", imp_index, "...\n")
    
    # --- Data Preparation per dataset ---
    dat_i <- dat_i_list[[i]] %>%
      mutate(
        age_Consultation = scale(age_Consultation),
        maxdiam = scale(maxdiam),
        mR_Baseline = as.numeric(mR_BaselineMICE) - 1,
        mR_Discharge = as.numeric(mR_Discharge) - 1,
        mR_1 = as.numeric(mR_1) - 1
      ) %>%
      select(all_of(keep)) %>%
      na.omit()
    
    # --- Run MCMC model ---
    modelMCMC <- mcmcmABN(
      FILENAMEbase = paste0(outdir, "/"),
      dat = dat_i,
      dists = dists,
      banmat = banmatASIS,
      retainmat = retainmatASIS,
      FILENAME = paste0(runname, "_imp", imp_index),
      max.par = modelmaxparents$max.par,
      MCMC.SEEDS = MCMC.SEEDS,
      MCMC.SCHEME = MCMC.SCHEME,
      PROB.REV = 0.05,
      PROB.MBR = 0,
      MCMC.PRIOR = 1,
      SCORE = "bic"
    )
    
    # --- Postprocessing (burn-in + thinning) ---
    mcmc.out.list.burn <- postBURNin(modelMCMC$mcmc.out.list, burnin.length)
    mcmc.out.list.thin <- postTHINN(mcmc.out.list.burn, thinningsteps)
    
    # --- Extract DAGs from all 4 chains ---
    dag_arrays <- lapply(mcmc.out.list.thin, function(x) x$dags)
    list.mc.out.burn.thin.dag <- abind::abind(dag_arrays, along = 3)
    
    # --- Compute per-imputation consensus DAG ---
    dag_mean <- apply(list.mc.out.burn.thin.dag, 1:2, mean)
    
    # --- Save results ---
    saveRDS(modelMCMC, file.path(outdir, paste0("model_imp", imp_index, ".rds")))
    saveRDS(dag_mean, file.path(outdir, paste0("dag_imp", imp_index, ".rds")))
    
    # Store in memory too (optional)
    model_list[[paste0("imp", imp_index)]] <- modelMCMC
    dag_list[[paste0("imp", imp_index)]] <- dag_mean
  }
  
  return(list(models = model_list, dags = dag_list))
}

compute_consensus_dag_from_saved <- function(runname) {
  outdir <- file.path("results", runname)
  
  # List all DAG files
  dag_files <- list.files(outdir, pattern = "^dag_imp[0-9]+\\.rds$", full.names = TRUE)
  if (length(dag_files) == 0) stop("No dag_imp*.rds files found in ", outdir)
  
  # Read them all
  dag_list <- lapply(dag_files, readRDS)
  names(dag_list) <- gsub("^dag_|\\.rds$", "", basename(dag_files))
  
  # Check dimension consistency
  dims <- sapply(dag_list, function(x) paste(dim(x), collapse = "x"))
  if (length(unique(dims)) > 1)
    stop("DAGs do not have identical dimensions. Check input data consistency.")
  
  # Compute global consensus DAG
  cat("\n▶ Computing global consensus DAG from", length(dag_list), "imputations...\n")
  consensus_dag <- Reduce("+", dag_list) / length(dag_list)
  
  # Save consensus DAG
  saveRDS(consensus_dag, file.path(outdir, "consensus_dag.rds"))
  
  return(list(dags = dag_list, consensus = consensus_dag))
}



compute_consensus_fit <- function(results, consensus_dag, threshold = 0.5) {
  # Convert consensus DAG to binary adjacency matrix
  consensus_dag_bin <- ifelse(consensus_dag > threshold, 1, 0)
  
  # Fit consensus DAG on each imputed dataset
  fit_list <- lapply(results$models, function(model_i) {
    dat_i <- dat_i_list[[i]] %>%
      mutate(
        age_Consultation = scale(age_Consultation),
        maxdiam = scale(maxdiam),
        mR_Baseline = as.numeric(mR_BaselineMICE) - 1,
        mR_Discharge = as.numeric(mR_Discharge) - 1,
        mR_1 = as.numeric(mR_1) - 1
      ) %>%
      select(all_of(keep)) %>%
      na.omit()
    fabn <- fitAbn(
      dag = consensus_dag_bin,
      data.df = dat_i,
      data.dists = dists,
      method = "mle"
    )
    return(fabn)
  })
  
  # Pool parameter estimates across imputations
  param_names <- names(fit_list[[1]]$coef)
  pooled_estimates <- sapply(param_names, function(param) {
    mean(sapply(fit_list, function(f) f$coef[[param]]))
  })
  
  # Return a list with all fits and pooled consensus
  return(list(
    fit_per_imputation = fit_list,
    pooled_consensus_estimates = pooled_estimates,
    consensus_dag_bin = consensus_dag_bin
  ))
}

# Analysis:
analyze_MCMC_ABN_results <- function(modelMCMC, dag, dat, dists, banmatASIS, retainmatASIS,
                                     modelmaxparents, runname, impnumber,burnin.length, thinningsteps) {
  cat("\n▶ Analysis started for imputation", impnumber, "\n")
  # Extract mcmc outputs again for quality plots
  mcmc.out.list.burn <- postBURNin(modelMCMC$mcmc.out.list, burnin.length)
  mcmc.out.list.thin <- postTHINN(mcmc.out.list.burn, thinningsteps)
  mcmc.out.list.burn.thin <- mcmc.out.list.thin
  
  # Create score chains for diagnostics
  mc.out.burn.thin.score.list <- lapply(mcmc.out.list.burn.thin, function(x) mcmc(x$scores))
  list.mc.out.burn.thin.score <- mcmc.list(mc.out.burn.thin.score.list)
  
  mc.out.score.list <- lapply(modelMCMC$mcmc.out.list, function(x) mcmc(x$scores))
  list.mc.out.score <- mcmc.list(mc.out.score.list)
  
  ## prepare the data again:
  dat <- dat %>%
    mutate(age_Consultation = scale(age_Consultation),
           maxdiam = scale(maxdiam),
           mR_Baseline = as.numeric(mR_BaselineMICE) - 1,
           mR_Discharge = as.numeric(mR_Discharge) - 1,
           mR_1 = as.numeric(mR_1) - 1) %>%
    select(all_of(keep)) %>%
    na.omit()
  
  # create the repo for the plots:
  dir.create(paste0("./results/", runname, "/", impnumber, "/"), recursive = TRUE, showWarnings = FALSE)
  
  # Run quality assessment
  MCMCQualityAssessment(mcmcoutput = list.mc.out.score,
                        mcmcoutputburnedthinned = list.mc.out.burn.thin.score,
                        FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"),
                        filename = paste0("GelmanPlot_",impnumber ))
  
  # Trace plots
  TracePlot(mcouthin1 = mcmc.out.list.burn.thin[[1]],
            mcoutthin2 = mcmc.out.list.burn.thin[[2]],
            mcoutthin3 = mcmc.out.list.burn.thin[[3]],
            mcoutthin4 = mcmc.out.list.burn.thin[[4]],
            MCMC.SEEDS = c("0815", "0510", "1107", "2611"),
            dat = dat, dists = dists,
            METHOD = "mle", modeloutput = modelMCMC,
            filename = paste0("TracePlot_thinned", thinningsteps, "_burned", burnin.length, "_", impnumber),
            FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"))

  # Arc strength significance
  arc.stren.sign.threshold <- arc.stren.threshold(dag, method = "l1")
  
  # ECDF plots
  ecdf_plot_ggplot(strengthdag = dag,
                   arc.stren.sign.threshold = arc.stren.sign.threshold,
                   filename = paste0("ECDF_arc_sten_threshold2_", impnumber),
                   FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"))
  
  # Threshold-specific plots
  plotDAGthresholds(threshold = arc.stren.sign.threshold,
                    dag = dag, dists = dists, dat = dat,
                    FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"),
                    filename = paste0("DAG_arcstreng_threshold_", impnumber))
  
  ThresholdHeatmap(dag = dag, banmat = banmatASIS,
                   threshold = arc.stren.sign.threshold,
                   FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"),
                   filename = paste0("ThresholdHeatmap", impnumber))
  
  AdjacencyMatrix(banmat = banmatASIS,
                  dag = dag,
                  threshold = arc.stren.sign.threshold,
                  filename = paste0("AdjacencyMatrix",impnumber),
                  FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"))
  
  abnMCMCabnComparision(threshold = arc.stren.sign.threshold,
                        dag = dag, dists = dists, dat = dat,
                        prior = 1, banmat = banmatASIS, retainmat = retainmatASIS,
                        max.parents = modelmaxparents$max.par,
                        FILENAMEbase = paste0("./results/", runname, "/", impnumber, "/"),
                        filename = paste0("DAG_Comparison", impnumber))
  
  # CPDAG plots
  dag_bin <- dag
  dag_bin[dag_bin > arc.stren.sign.threshold] <- 1
  dag_bin[dag_bin <= arc.stren.sign.threshold] <- 0
  
  colnames(dag_bin) <- rownames(dag_bin) <- names(dat)
  
  png(paste0("results/", runname,"/",  impnumber ,"/CPDAG_", impnumber, ".png"), width = 14*600, height = 9*600, res = 600)
  par(mfrow = c(1, 2))
  
  dagBN <- bnlearn::empty.graph(names(dat))
  bnlearn::amat(dagBN) <- t(dag_bin)
  plot(dagBN, main = "DAG")
  
  cons.cpdag <- bnlearn::cpdag(dagBN)
  cons.amat <- bnlearn::amat(cons.cpdag)
  cons.amat[cons.amat == 1 & t(banmatASIS) == 1] <- 0
  
  cpdag.cons.amat <- bnlearn::empty.graph(names(dat))
  bnlearn::amat(cpdag.cons.amat) <- cons.amat
  plot(cpdag.cons.amat, main = "CPDAG II")
  dev.off()
  
  # SVG for presentation
  png(paste0("results/", runname,"/",  impnumber ,"/MCMC_", impnumber, ".png"), width = 14*600, height = 9*600, res = 600)
  par(mfrow = c(1, 1))
  edgestren <- dag
  edgestren[edgestren <= arc.stren.sign.threshold] <- 0
  dag_threshold <- dag
  dag_threshold[dag_threshold > arc.stren.sign.threshold] <- 1
  dag_threshold[dag_threshold <= arc.stren.sign.threshold] <- 0
  
  plotAbn(dag = dag_threshold, data.df = dat, data.dists = dists,
          digits = 2, edge.strength = edgestren, plot = TRUE,
          main = arc.stren.sign.threshold)
  dev.off()
  
  cat("\n▶ Analysis completed for imputation", impnumber, "!\n")}

load_MCMC_ABN_results <- function(runname) {
  outdir <- file.path("results", runname)
  
  # --- Find saved model and DAG files ---
  model_files <- list.files(outdir, pattern = "^model_imp[0-9]+\\.rds$", full.names = TRUE)
  dag_files   <- list.files(outdir, pattern = "^dag_imp[0-9]+\\.rds$", full.names = TRUE)
  
  if (length(model_files) == 0 && length(dag_files) == 0)
    stop("No saved model or DAG .rds files found in ", outdir)
  
  # --- Read all models and DAGs ---
  model_list <- lapply(model_files, readRDS)
  dag_list   <- lapply(dag_files, readRDS)
  
  # Ensure consistent ordering by imputation number
  get_imp_number <- function(fname) as.numeric(gsub("\\D", "", basename(fname)))
  order_idx_models <- order(sapply(model_files, get_imp_number))
  order_idx_dags   <- order(sapply(dag_files, get_imp_number))
  
  model_list <- model_list[order_idx_models]
  dag_list   <- dag_list[order_idx_dags]
  
  names(model_list) <- paste0("imp", sapply(model_files[order_idx_models], get_imp_number))
  names(dag_list)   <- paste0("imp", sapply(dag_files[order_idx_dags], get_imp_number))
  
  # --- Try to read global consensus DAG if available ---
  consensus_path <- file.path(outdir, "consensus_dag.rds")
  consensus_dag <- if (file.exists(consensus_path)) readRDS(consensus_path) else NULL
  
  # --- Return full results structure identical to run_MCMC_ABN_all_imputations ---
  result_list <- list(
    models = model_list,
    dags = dag_list,
    consensus = consensus_dag
  )
  
  cat("✅ Successfully loaded", length(model_list), "models and",
      length(dag_list), "DAGs from", outdir, "\n")
  
  if (!is.null(consensus_dag)) {
    cat("✅ Consensus DAG found and loaded.\n")
  } else {
    cat("⚠️ No consensus DAG found. You can compute it using compute_consensus_dag_from_saved().\n")
  }
  
  return(result_list)
}

