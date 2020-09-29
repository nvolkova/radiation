library(VariantAnnotation)
source('useful_functions.R')
source('plotting_functions.R')

data <- read.csv("sample_annotations.csv")
data$Sample <- as.character(data$Sample)
data$Genotype.new <- as.character(data$Genotype.new)
data$Code <- as.character(data$Code)
CD2Mutant <- sapply(data$Code, function(x) {
  t <- unlist(strsplit(x,split="[:]"))
  t[t=="NA"] <- ""
  if (t[4]!="") # mutagen exposure
    return(paste(t[3],substr(t[4],1,3),t[5],t[7],sep=":")) # genotype, mutagen, dose, experiment type, generation
  if (t[4]=="") # mutation accumulation
    return(paste(t[3],t[7],sep=":")) # genotype, experiment type, generation
})
names(CD2Mutant) <- data$Sample
CD2Mutant <- CD2Mutant[-c(grep("INCORRECT",data$Genotype.check))]
worms <- names(CD2Mutant)
worms <- sort(worms)
CD2Mutant <- CD2Mutant[worms]
data <- data[match(worms,data$Sample),]

mut.acc <- names(CD2Mutant)[is.na(data$Mutagen)]
normal_panel = names(CD2Mutant)[grep("N2:1",CD2Mutant[1:1000])][1:6]
mut.acc <- setdiff(mut.acc, normal_panel)

# Upload indels
indels_dedup <- sapply(worms, function(x) readVcf(paste0('PATH/TO/INDEL/VCFS/',x,'vcf.gz')))
# Upload subs
indels_dedup <- sapply(worms, function(x) readVcf(paste0('PATH/TO/SUBS/VCFS/',x,'vcf.gz')))
# Upload SVs
SVclust.new <- sapply(worms, function(x) readVcf(paste0('PATH/TO/SV/TABLES/',x,'.tsv')))

svmat <- do.call('rbind',SVclust.new)
svmat <- svmat[as.character(svmat$CHR1) == as.character(svmat$CHR2),]

svmat$CHR1 <- as.character(svmat$CHR1)
svmat$CHR2 <- as.character(svmat$CHR2)
svmat$Sample <- as.character(svmat$Sample)
svmat.new <- svmat[,c(1:4,12:13)]
svmat <- svmat[svmat$clust.type != 'some',]

for (w in unique(svmat$Sample)) {
  
  tmp <- svmat[svmat$Sample == w,]
  tmp1 <- tmp[,c(1:4,12:13)]
  
  for (j in unique(tmp$clust)) {
    
    clust.tmp <- tmp[tmp$clust == j,]
    
    tmp1 <- rbind(tmp1,
                  c(as.character(clust.tmp$CHR1[1]), min(as.numeric(clust.tmp$POS1)), 
                    as.character(clust.tmp$CHR2[1]), max(as.numeric(clust.tmp$POS2)),
                    clust.tmp$Sample[1], clust.tmp$clust.type[1]))
    
  }
  
  tmp1 <- tmp1[-c(1:nrow(tmp)),,drop = F]
  svmat.new <- svmat.new[-which(svmat.new$Sample == w),]
  svmat.new <- rbind(svmat.new, tmp1)
  print(w)
}
svmat.new$POS1 <- as.numeric(svmat.new$POS1)
svmat.new$POS2 <- as.numeric(svmat.new$POS2)
svmat <- svmat.new

chr_lens <- c(15072434,15279421,13783801,17493829,20924180,17718942)
names(chr_lens) <- c('I','II','III','IV','V','X')
df <- data.frame(name = names(chr_lens), length = chr_lens)

GG_save_pdf = function(list, filename, ...) {
  #start pdf
  pdf(filename, ...)
  #loop
  count <- 1
  for (p in list) {
    print(p)
    print(count)
    count = count +1
  }
  #end pdf
  dev.off()
  invisible(NULL)
}

# Clustering functions
isClustered <- function(vcf, p=1e-3, q=0.1, r=100){
  d <-  diff(vcf$start)
  w <- d > 1 & diff(as.numeric(vcf$seqnames)) == 0
  #	p <- 1e-3 # P N>Kat
  #	q <- 0.05 # P Kat>N
  P <- matrix(c(1-p,p,q, 1-q), ncol=2, byrow=TRUE) # Transition probabilities matrix
  p0 <- c(1,0)
  s <- c(mean(d[w]), r)
  dw <- d[w]
  l <- length(dw)
  T1 <- T2 <- matrix(0,ncol=l, nrow=2)
  T1[,1] <- log(c(q/(q+p), p/(q+p))) # log of Vitterbi path probabilities
  lP <- log(P)
  dg <- rbind(dgeom(dw, prob=1/s[1], log=TRUE), dgeom(dw, prob=1/s[2], log=TRUE)) 
  # state observation loglikelihood given the current state (none, kataegis)
  
  # Viterbi algorithm
  
  for(i in 2:l){
    x <- ((T1[,i-1] + lP) + dg[,i])
    T2[1,i] <- (x[1,1] < x[2,1])+1
    T2[2,i] <- (x[1,2] < x[2,2])+1
    T1[1,i] <- x[T2[1,i],1]
    T1[2,i] <- x[T2[2,i],2]
    #x <- T1[,i-1] + lP + rep(dg[,i],each=2) # previous Vitterbi path probability * transition probability * observation likelihood
    #T1[,i] <- sapply(1:2, function(state) max(x[,state]))
    #T2[,i] <- sapply(1:2,function(state) which.max(x[,state])) # most probable states - backpointer
  }
  finalT1 <- max(T1[,l] + log(c(q/(q+p), p/(q+p)))) # + probability of transition to final state, let's say 1 and then log(1)=0
  finalT2 <- which.max(T1[,l] + log(c(q/(q+p), p/(q+p))))# + 
  z <- numeric(l)
  z[l] <- finalT2
  z[l] <- 1 # this means that the backtrace starts from 1 (not necessarily)
  for(i in l:2){
    z[i-1] <- T2[z[i],i]
  }
  k <- numeric(nrow(vcf))
  k[-1][w][-1] <- z[-l]-1
  k[-nrow(vcf)][w][-1] <- (z[-l]-1) | k[-nrow(vcf)][w][-1]
  
  # Other clustered
  pc <- pgeom(dw, prob=1/s[1], lower.tail=TRUE)
  qc <- p.adjust(pc, "BH") < 0.05
  
  cl <- numeric(nrow(vcf))
  cl[-1][w] <- qc
  cl[-nrow(vcf)][w] <- qc | cl[-nrow(vcf)][w]
  
  clk <- factor(pmin(cl + 2*k,2), levels=0:2, labels=c("None","Clustered","Kataegis"))
  return(clk)
}

# Select samples to show
experiments <- unique(CD2Mutant[which(data$Mutagen == 'Radiation' & data$Type == 'mutagen' & data$Drug.concentration > 0)])

huge_tot <- list()
totcount <- 1
q <- list()
j <- 1
clrs <- c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta")
for (ex in experiments) {
  
  samset <- names(CD2Mutant)[CD2Mutant == ex & data$Type == 'mutagen']
    
  tot <- list()
    
  for (n in samset) {
      
    tmp_vcf <- vcfs_dedup[[n]]
    if (length(tmp_vcf) > 0) {
      #mcols(tmp_vcf)$Sample <- n
      tmp <- isMNV(tmp_vcf)
      tmp_vcf_df <- as.data.frame(granges(tmp_vcf))
      tmp_vcf_df$Sample <- n
      tmp_vcf_df$VAF <- geno(tmp_vcf)[['PM']][,2]
      tmp_vcf_df$Type <- NA
      if (sum(tmp)>0) {
        dnvs <- tmp_vcf_df[tmp,][seq(1,sum(tmp),2),]
        dnvs$Type <- 'DNV'
      } else {
        dnvs <- tmp_vcf_df[tmp,]
      }
      subs <- tmp_vcf_df[!tmp,]
      if (nrow(subs) >0) {
        subs <- subs[subs$seqnames!='MtDNA',]
        subs$Type <- paste0(subs$REF,'>',unlist(sapply(subs$ALT,as.character)))
        subs$Type[subs$Type == 'A>C'] <- 'T>G'
        subs$Type[subs$Type == 'A>G'] <- 'T>C'
        subs$Type[subs$Type == 'A>T'] <- 'T>A'
        subs$Type[subs$Type == 'G>A'] <- 'C>T'
        subs$Type[subs$Type == 'G>C'] <- 'C>G'
        subs$Type[subs$Type == 'G>T'] <- 'C>A'
      }
    } else {
      subs <- tmp_vcf_df
      dnvs <- tmp_vcf_df
    }
      
    vcf <- indels_dedup[[n]]
    indels <- as.data.frame(granges(vcf))
    if (nrow(indels) > 0) {
      indels <- indels[indels$seqnames!='MtDNA',]
      indels$Sample <- n
      indels$VAF <- (geno(vcf)$PU[,"TUMOUR"] + geno(vcf)$NU[,"TUMOUR"]) / 
        (geno(vcf)$PR[,"TUMOUR"] + geno(vcf)$NR[,"TUMOUR"]) 
      indels$Type <- ifelse(nchar(indels$REF) > nchar(unlist(sapply(indels$ALT,as.character))), yes = 'D', no = 'I')
      indels$Type[nchar(indels$REF) > 1 & nchar(unlist(sapply(indels$ALT,as.character))) > 1] <- 'DI'
    }
    
    if (nrow(subs) > 0) {
      if (nrow(indels) > 0) {
        if (nrow(dnvs) >0) tot[[n]] <- rbind(subs,indels,dnvs)
        else tot[[n]] <- rbind(subs,indels)
      } else {
        if (nrow(dnvs) >0) tot[[n]] <- rbind(subs,dnvs)
        else tot[[n]] <- subs
      }
    } else {
      if (nrow(indels) > 0) {
        if (nrow(dnvs) >0) tot[[n]] <- rbind(indels,dnvs)
        else tot[[n]] <- indels
      } else {
        if (nrow(dnvs) >0) tot[[n]] <- rbind(dnvs)
        else tot[[n]] <- NULL
      }
    }
      
    if (length(tot[[n]])>0) {
      tot[[n]]$Mode <- NA
      tot[[n]]$experiment <- CD2Mutant[n]
      if (nrow(subs)>0) 
        tot[[n]]$Mode[1:nrow(subs)] <- 'A'
      if (nrow(indels) >0)
        tot[[n]]$Mode[(nrow(subs)+1):(nrow(subs) + nrow(indels))] <- 'B'
      if (nrow(dnvs)>0)
        tot[[n]]$Mode[(nrow(subs)+nrow(indels)+1):nrow(tot[[n]])] <- 'C'
      rownames(tot[[n]]) <- NULL
      tot[[n]] <- tot[[n]][order(tot[[n]]$seqnames),]
      for (ch in levels(tot[[n]]$seqnames)) {
        tot[[n]][tot[[n]]$seqnames == ch,] <- tot[[n]][tot[[n]]$seqnames == ch,][order(tot[[n]]$start[tot[[n]]$seqnames == ch]),]
      }
      if (max(table(tot[[n]]$seqnames)) < 3) {
        tot[[n]]$clust <- 1
      } else {
        k <- isClustered(tot[[n]])
        tot[[n]]$clust <- as.character(k)
        tot[[n]]$clust[tot[[n]]$clust == 'Kataegis'] <- 5
        tot[[n]]$clust[tot[[n]]$clust != 5] <- 1
      }
    }
  }
  
  tot <- do.call('rbind', tot)
  huge_tot[[totcount]] <- tot
  totcount <- totcount + 1
  
  for (lll in unique(tot$start)) {
    if (sum(tot$start == lll) == 1) next
    if (var(tot$VAF[tot$start == lll]) < 0.01 & length(unique(tot$Type[tot$start == lll])) == 1)
      tot <- tot[-which(tot$start == lll)[-1],]
  }  
    
  if (sum(tot$clust>1)>0)
    q[[j]] <- ggplot() + geom_bar(data = df, aes(x = name, y = length), stat = 'identity',fill = 'white') +
      scale_y_continuous(labels = c('0 Mb','5 Mb', '10 Mb','15 Mb','20 Mb')) +
      geom_jitter(data = tot, aes(x = seqnames, y = start, col = Type, shape =Mode, size = clust)) +
      labs(title = paste0('Mutations across all samples from ',CD2Mutant[samset[1]], ' experiment'),
           x = 'Chromosome',y='Position') +
      scale_shape_discrete(labels = c('Substitutions','Indels','DNVs')) +
      scale_size_manual(values = c(1,5), labels = c('single','clustered')) +
      scale_color_manual(values=clrs[c(1:3,8:9,7,10,4:6)]) +
      guides(size=guide_legend(title="Clustering"), shape = guide_legend(title='Class')) + 
    theme(title = element_text(size=10))
  else
    q[[j]] <- ggplot() + geom_bar(data = df, aes(x = name, y = length), stat = 'identity',fill = 'white') +
      scale_y_continuous(labels = c('0 Mb','5 Mb', '10 Mb','15 Mb','20 Mb')) +
      geom_jitter(data = tot, aes(x = seqnames, y = start, col = Type, shape = Mode)) +
      labs(title = paste0('Mutations across all samples from ',CD2Mutant[samset[1]], ' experiment'),
           x = 'Chromosome',y='Position') +
      scale_shape_discrete(labels = c('Substitutions','Indels','DNVs')) +
      scale_color_manual(values=clrs[c(1:3,8:9,7,10,4:6)]) +
      guides(shape = guide_legend(title='Class')) + 
      theme(title = element_text(size=10))
  
  
  COLOR=c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD")
  names(COLOR) = c("TD", "DEL", "INV", "COMPLEX", "FOLDBACK", "MOVE", "TRSL")
  svs <- svmat[as.character(svmat$Sample) %in% samset,]
  svs$POS1 <- as.numeric(svs$POS1)
  svs$POS2 <- as.numeric(svs$POS2)
  svs$CHR1 <- match(svs$CHR1, df$name)
  svs$clust.type <- factor(svs$clust.type)
  
  if (nrow(svs)>0) {
    q[[j]] <- q[[j]] + geom_rect(data = svs,
                                 mapping = aes(xmin = CHR1-0.5, 
                                               xmax = CHR1+0.5, 
                                               ymin = POS1, 
                                               ymax = POS2+10000,
                                               fill = clust.type),
                                 alpha = 0.4) +
      scale_fill_manual(values = COLOR) + 
      guides(fill = guide_legend(title='SV'))
  }
  
  j = j+1
  print(j)
}

new_order <- order(experiments)
q <- q[new_order[c(51:57,1:6,28:29,58:62,103:105,
                   108:109, 22:23, 110:112, 30:35,
                   65:66,87:88,10:13,48,97:100,76:78,82:84,79:81,
                   89:90, 93:94,43:45,75,24:25,36,39,85:86,106:107,14:15,95:96,
                   18:21,26:27,46:47,91:92)]]
GG_save_pdf(q, filename = '~/Irradiation_location_plots.pdf', 8,7)

allclust <- do.call('rbind',huge_tot)
allclust <- allclust[allclust$clust>1,]

allclust$genotype <- data$Genotype.new[match(allclust$Sample, data$Sample)]


number_per_cluster <- list()
cluster_sizes <- list()
for (gene in unique(allclust$genotype)) {
  
  distances <- diff(allclust$start[allclust$genotype == gene])
  rle(as.numeric(abs(distances) < 100)) -> dist_rle
  number_per_cluster[[gene]] <- dist_rle$lengths[dist_rle$values == 1] + 1
  cluster_sizes[[gene]] <- NULL
  cur_clust <- 0
  for (j in 1:length(distances)) {
    if (abs(distances[j]) < 100)
      cur_clust <- cur_clust + distances[j]
    else {
      if (cur_clust > 0)
        cluster_sizes[[gene]] <- c(cluster_sizes[[gene]], cur_clust)
      cur_clust <- 0
    }
  }
  if (length(distances) == 1 & cur_clust > 0) cluster_sizes[[gene]] <- c(cluster_sizes[[gene]], cur_clust)
  
}

number_per_cluster$`rad-54B` <- 0
number_per_cluster$`ndx-4` <- 0
number_per_cluster$`him-6` <- 0
names(number_per_cluster)[43] <- 'rad-54.B (gt3308)'
names(number_per_cluster)[52] <- 'rad-54.B (gk340656)'
names(number_per_cluster)[44] <- 'bub-3 (gt2000)'
names(number_per_cluster)[45] <- 'bub-3 (ok3437)'
names(number_per_cluster)[23] <- 'tdpo-1'

correct_order <- c(1,5:7,13,53,31:32,23:24,
                   36,11,37:38,14,26,35,
                   34,33,19,3,8,30,22,42,
                   43,52,20,21,2,15,12,28,
                   54,18,25,44:46,9:10)

cluster_sizes$`rad-54B` <- 0
cluster_sizes$`ndx-4` <- 0
cluster_sizes$`him-6` <- 0
names(cluster_sizes)[43] <- 'rad-54.B (gt3308)'
names(cluster_sizes)[52] <- 'rad-54.B (gk340656)'
names(cluster_sizes)[44] <- 'bub-3 (gt2000)'
names(cluster_sizes)[45] <- 'bub-3 (ok3437)'
names(cluster_sizes)[23] <- 'tdpo-1'

library(beeswarm)
pdf('~/Desktop/Irradiation analysis/Mutations_per_cluster_across_genotypes.pdf',10,6)
par(mar = c(10,4,4,2))
beeswarm(number_per_cluster[correct_order], las = 2, bty = 'n', method = 'hex', ylab = 'Number of mutations / cluster',
         ylim = c(0,10), cex = 0.5, pch = 16, corral = 'gutter', main = 'Sizes of IR-induced clusters across genotypes (muts)',
         xaxt = 'n')
axis(side = 1, at = c(1:41), labels = names(number_per_cluster)[correct_order], 
     font = 3, las = 2, lty = 0)
dev.off()
pdf('~/Desktop/Irradiation analysis/Cluster_span_across_genotypes.pdf',10,6)
par(mar = c(10,4,4,2))
boxplot(cluster_sizes[correct_order], las = 2, frame = F, xaxt = 'n',
        main = 'Sizes of IR-induced clusters across genotypes (bps)',
        ylim = c(0,200), ylab = 'bp spanned by cluster', pch = 16, 
        cex = 0.5, outline = F)
axis(side = 1, at = c(1:41), labels = names(cluster_sizes)[correct_order], 
     font = 3, las = 2, lty = 0)
beeswarm(cluster_sizes[correct_order], las = 2, bty = 'n', method = 'hex', 
         ylab = 'bp spanned by cluster', col = 'gray56', add = T, 
         ylim = c(0,200), cex = 0.5, pch = 16, corral = 'gutter', 
         main = 'Sizes of IR-induced clusters across genotypes (bps)')
dev.off()


CD2Mutant.reduced <- sapply(CD2Mutant, function(x) {
  tmp <- unlist(strsplit(x, split = '[:]'))
  if (length(tmp) == 2) return(x)
  else return(paste(tmp[1:3],collapse=':'))
})

clust_summary <- data.frame(name = unique(allclust$Sample))
clust_summary$name <- as.character(clust_summary$name)
clust_summary$generation <- data$Generation[match(clust_summary$name,data$Sample)]
clust_summary$genotype <- data$Genotype.new[match(clust_summary$name,data$Sample)]
clust_summary$exposure <- data$Mutagen[match(clust_summary$name,data$Sample)]
clust_summary$dose <- data$Drug.concentration[match(clust_summary$name,data$Sample)]
clust_summary$code <- CD2Mutant.reduced[match(clust_summary$name,names(CD2Mutant.reduced))]

clust_summary$number_of_clusters <- sapply(clust_summary$name, function(x) sum(diff(allclust$start[allclust$Sample == x])>1000) + 1)
clust_summary$number_clust_muts <- sapply(clust_summary$name, function(x) sum(allclust$Sample == x))
clust_summary$total_samples_of_this_code <- sapply(clust_summary$code, function(x) sum(CD2Mutant.reduced == x))
clust_summary$total_mutations_in_sample <- sapply(clust_summary$name, function(x) sum(Y[x,1:112]))

write.csv(clust_summary, file = '~/Desktop/Irradiation analysis/Clusters_in_irradiated_samples_new.csv')
#################### now some stats ##########################

load('~/yoda2/IR/IR_adj.RData')

#experiments <- unique(CD2Mutant[which(data$Mutagen == 'Radiation' & data$Drug.concentration>0 & data$Type == 'mutagen')])

clust <- sapply(rownames(Y.IR)[which(data.IR$Mutagen == 'Radiation' & data.IR$Drug.concentration>0 )], 
                function(x) {
                  if (x %in% clust_summary$name) 
                    return(clust_summary$number_of_clusters[clust_summary$name == x])
                  else 
                    return(0)
          })
mean(clust[grep('exo-1:Rad:0',CD2Mutant[names(clust)])])
mean(clust[grep('exo-1:Rad:10',CD2Mutant[names(clust)])])
mean(clust[grep('exo-1:Rad:20',CD2Mutant[names(clust)])])
mean(clust[grep('N2:Rad:40',CD2Mutant[names(clust)])])
mean(clust[grep('N2:Rad:80',CD2Mutant[names(clust)])])

CD2Mutant.very.reduced <- sapply(CD2Mutant[names(clust)], function(x) {
  tmp <- unlist(strsplit(x, split = '[:]'))
  #if (length(tmp) == 2 | tmp[1] == 'N2') return('non')
  #else 
  return(paste(tmp[1:2],collapse=':'))
})

doses <- data.IR$Drug.concentration[match(names(clust), data.IR$Sample)]
interactions <- t(sapply(names(clust), function(x) as.numeric(unique(CD2Mutant.very.reduced) == CD2Mutant.very.reduced[match(x,names(CD2Mutant.very.reduced))])))
colnames(interactions) <- unique(CD2Mutant.very.reduced)

cd <- data.frame(interactions * doses)

cd <- cd[,-match(c('exo.1.Rad','exo.3.Rad','apn.1.Rad',
                   'parp.1.Rad','pole.4.Rad','rcq.5.Rad',
                   'ndx.4.Rad','agt.1.Rad','tdp.1.Rad'), colnames(cd))]

cd1 <- cd
cd1[cd1>0] <- 1

simplemodel <- glm( clust ~ offset(log(doses)) + 0 + .,data = cd1, family = stats::poisson())
coef(summary(simplemodel))
which(p.adjust(coef(summary(simplemodel))[,4], method = 'BH')<0.05) # everything

coeffs <- coef(summary(simplemodel))

pvclust <- NULL
for (zzz in rownames(coeffs)[-1]) {
  stat_mu = coeffs[zzz,1] - coeffs['N2.Rad',1]
  stat_sd = sqrt(coeffs[zzz,2]**2 + coeffs['N2.Rad',2]**2)
  zscore = stat_mu / stat_sd
  pvclust <- c(pvclust, 1 - pchisq(q = zscore**2, df = 1))
}
rownames(coeffs)[-1][which(p.adjust(pvclust,method='BH') < 0.05)]

rates <- exp(coeffs[,1])
rates.sd <- rates * coeffs[,2]

# Run a model on proportions (without generations)

prop.of.clust <- sapply(names(clust), function(x) {
  if (x %in% clust_summary$name) return(clust_summary$number_clust_muts[clust_summary$name == x] / sum(Y[x,c(1:112)]))
  else return(0)
})

simplemodel <- glm( prop.of.clust ~ 0 + .,data = cd1, family = gaussian())
coef(summary(simplemodel))
which(p.adjust(coef(summary(simplemodel))[,4], method = 'BH')<0.05)

coeffs <- coef(summary(simplemodel))

pvclust <- NULL
for (zzz in rownames(coeffs)[-1]) {
  stat_mu = coeffs[zzz,1] - coeffs['N2.Rad',1]
  stat_sd = sqrt(coeffs[zzz,2]**2 + coeffs['N2.Rad',2]**2)
  zscore = stat_mu / stat_sd
  pvclust <- c(pvclust, 1 - pchisq(q = zscore**2, df = 1))
}
rownames(coeffs)[-1][which(p.adjust(pvclust,method='BH') < 0.05)]

prop.rates <- coeffs[,1]
prop.rates.sd <- coeffs[,2]

# Visualize
# no. of clusters
set.seed(111)
par(mfrow = c(1,2))
boxplot(rates*80, frame = F, outline = F, ylim = c(0,max(rates*80+1.96*rates.sd*80)),
        ylab = 'No. of clusters per 80 Gy', main = 'Clusters per 80 Gy')
newx <- jitter(rep(1, length(rates)), amount = 0.1)
points(x = newx, y = rates*80, col = 'gray', pch = 16) 
#o <- which(p.adjust(pvclust,method='BH') < 0.1)
points(x = newx[1], y = rates[1]*80, col = 'darkred', pch = 16) 
#points(x = newx[o+1], y = rates[o+1], col = 'darkred', pch = 16) 
arrows(x0 = newx[1], 
       y0 = rates[1]*80 - 1.96*80*rates.sd[1], y1 = rates[1]*80 + 1.96*80*rates.sd[1],
       col = 'gray21',lwd=0.5,length=0)
abline(h = 80*rates['N2.Rad'], lty = 2)
#text(x = c(newx[o+1][1] - 0.2, newx[o+1][-1] + 0.2), y = rates[o+1], font = 3,
#     labels = names(rates)[o+1], cex = 0.7)
text(x = c(newx[which(rates > 2 * rates['N2.Rad'])] - 0.2), 
     y = 80*rates[which(rates > 2 * rates['N2.Rad'])], font = 3,
     labels = names(which(rates > 2 * rates['N2.Rad'])), cex = 0.7)
# proportion of clustered muts
boxplot(prop.rates, frame = F, outline = F, ylim = c(0,max(prop.rates+1.96*prop.rates.sd)),
        ylab = 'Prop. of clustered mut-s', main = 'Proportion of clustered mutations\n across genotypes')
newx2 <- jitter(rep(1, length(prop.rates)), amount = 0.1)
points(x = newx2, y = prop.rates, col = 'gray', pch = 16) 
#o2 <- which(p.adjust(pvprop,method='BH') < 0.1)
points(x = newx2[1], y = prop.rates[1], col = 'darkred', pch = 16) 
arrows(x0 = newx2[1],y0 = prop.rates[1] - 1.96*prop.rates.sd[1],
       y1=prop.rates[1]+1.96*prop.rates.sd[1],
       col = 'gray21',lwd=0.5,length=0)
abline(h = prop.rates['N2.Rad'], lty = 2)
#text(x = c(newx2[o2+1][1] - 0.2, newx2[o2+1][c(2:3)] + 0.22,newx2[o2+1][4] - 0.2,newx2[o2+1][5] - 0.25,newx2[o2+1][6:7] - 0.2),
#     y = prop.rates[o2+1], font = 3,
#     labels = print_names[sapply(names(prop.rates)[o2+1], function(x) grep(x,names(print_names))[1])],
#     cex = 0.7)
text(x = newx2[which(prop.rates > 2 * prop.rates['N2.Rad'])],
     y = prop.rates[which(prop.rates > 2 * prop.rates['N2.Rad'])], font = 3,
     labels = names(which(prop.rates > 2 * prop.rates['N2.Rad'])),
     cex = 0.7)
legend('topright', legend = 'significantly different\n from N2 (FDR 5%)', fill = 'darkred',bty = 'n',
       border = NA, cex = 0.7)


allclust <- do.call('rbind',huge_tot)
allclust <- allclust[allclust$clust>1,]

CD2Mutant.reduced <- sapply(CD2Mutant, function(x) {
  tmp <- unlist(strsplit(x, split = '[:]'))
  if (length(tmp) == 2) return(x)
  else return(paste(tmp[1:3],collapse=':'))
})

clust_summary <- data.frame(name = unique(allclust$Sample))
clust_summary$name <- as.character(clust_summary$name)
clust_summary$generation <- data$Generation[match(clust_summary$name,data$Sample)]
clust_summary$genotype <- data$Genotype.new[match(clust_summary$name,data$Sample)]
clust_summary$exposure <- data$Mutagen[match(clust_summary$name,data$Sample)]
clust_summary$dose <- data$Drug.concentration[match(clust_summary$name,data$Sample)]
clust_summary$code <- CD2Mutant.reduced[match(clust_summary$name,names(CD2Mutant.reduced))]

clust_summary$number_of_clusters <- sapply(clust_summary$name, function(x) sum(diff(allclust$start[allclust$Sample == x])>1000) + 1)
clust_summary$number_clust_muts <- sapply(clust_summary$name, function(x) sum(allclust$Sample == x))
clust_summary$total_samples_of_this_code <- sapply(clust_summary$code, function(x) sum(CD2Mutant.reduced == x))
clust_summary$total_mutations_in_sample <- sapply(clust_summary$name, function(x) sum(Y[x,1:112]))

clust_summary_mut <- clust_summary[!is.na(clust_summary$dose),]
clust_summary_mut <- clust_summary_mut[!(clust_summary_mut$exposure == 'Protonbeam'),]

################################################################################################################
################################### now some stats #############################################################

clust <- sapply(rownames(Y)[which(data$Drug.concentration[match(rownames(Y),data$Sample)]>0)], function(x) {
  if (x %in% clust_summary_mut$name) return(clust_summary_mut$number_of_clusters[clust_summary_mut$name == x])
  else return(0)
})

# need a design matrix cd
generation.function <- function(N) {
  if (is.na(N)) return(N)
  if (N==0) return(0)
  if (N==1) return(1)
  alpha = sum(sapply(1:N, function(i) 1/(2^(i-1)))) +
    0.25*sum(sapply(1:(N-1), function(i) sum(sapply(1:i, function(j) 1/(2^(j-1))))))
  return(alpha)
}

CD2Mutant.very.reduced <- sapply(CD2Mutant[names(clust)], function(x) {
  tmp <- unlist(strsplit(x, split = '[:]'))
  #if (length(tmp) == 2 | tmp[1] == 'N2') return('non')
  #else 
  return(paste(tmp[1:2],collapse=':'))
})

mutagens <- t(sapply(names(clust), function(x) as.numeric(unique(data$Mutagen)[-c(1,12)] == data$Mutagen[match(x,data$Sample)])))
colnames(mutagens) <- unique(data$Mutagen)[-c(1,12)]
mutagens[is.na(mutagens)] <- 0
mutagens <- mutagens[,colSums(mutagens)>0]
doses <- data$Drug.concentration[match(names(clust), data$Sample)]
#doses[is.na(doses)] <- 0
interactions <- t(sapply(names(clust), function(x) as.numeric(unique(CD2Mutant.very.reduced) == CD2Mutant.very.reduced[match(x,names(CD2Mutant.very.reduced))])))
colnames(interactions) <- unique(CD2Mutant.very.reduced)

cd <- data.frame(interactions * doses)

# median adjustment for mutagens
mutagens_in_interactions <- sapply(unique(CD2Mutant.very.reduced), function(x) unlist(strsplit(x,split='[:]'))[2])
median_per_mut <- apply(mutagens*doses,2,function(z) median(z[z>0]))
names(median_per_mut) <- unique(mutagens_in_interactions)
for (j in 1:ncol(interactions)) {
  cd[,j] <- cd[,j] / median_per_mut[mutagens_in_interactions[j]]
}

simplemodel <- glm( clust ~ ., data = cd, family = stats::poisson)
coef(summary(simplemodel))
which(p.adjust(coef(summary(simplemodel))[,4], method = 'BH')<0.05 & coef(summary(simplemodel))[,1]>0)

library(greta)
m <- ncol(cd)
# Run a model on the number of clusters vs genotype and generation
sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
mu = (cd %*% t(rates))
clust <- t(t(clust))
distribution(clust) = poisson(lambda = mu)
cl.model <- model(rates,sigma)
cl.draws <- mcmc(cl.model, warmup = 500, n_samples = 500) # do on cluster
draws_all <- do.call('rbind',cl.draws)
rates <- colMeans(draws_all[,1:m])
rates.sd <- apply(draws_all[,1:m],2,sd)
names(rates) = names(rates.sd) <- colnames(cd)

pvclust <- NULL
for (zzz in names(rates)) {
  stat_mu = rates[zzz]
  stat_sd = rates.sd[zzz]
  zscore = stat_mu / stat_sd
  pvclust <- c(pvclust, 1 - pchisq(q = zscore**2, df = 1))
}
which(p.adjust(pvclust,method='BH') < 0.05)

pvclust_ind <- NULL
for (zzz in names(rates)[-1]) {
  stat_mu = rates[zzz] - rates[1]
  stat_sd = sqrt(rates.sd[zzz]**2 + rates.sd[1]**2)
  zscore = stat_mu / stat_sd
  pvclust_ind <- c(pvclust_ind, 1 - pchisq(q = zscore**2, df = 1))
}
which(p.adjust(pvclust_ind,method='BH') < 0.05)

#################################

# Run a model on proportions (without generations)

prop.of.clust <- sapply(names(clust), function(x) {
  if (x %in% clust_summary$name) return(clust_summary$number_clust_muts[clust_summary$name == x] / sum(Y[x,c(1:112)]))
  else return(0)
})

sigma <- variable(lower = 0)
sigma2 <- variable(lower = 0)
prop.rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
cd1 <- cd / rowSums(cd)
mu = (cd1 %*% t(prop.rates))
prop.of.clust <- t(t(prop.of.clust))
distribution(prop.of.clust) = normal(mean = mu, sd = sigma2)
prop.model <- model(prop.rates,sigma,sigma2)
prop.draws <- mcmc(prop.model, warmup = 500, n_samples = 500) # do on cluster
draws_all <- do.call('rbind',prop.draws)
prop.rates <- colMeans(draws_all[,1:m])
prop.rates.sd <- apply(draws_all[,1:m],2,sd)
names(prop.rates) = names(prop.rates.sd) <- colnames(cd)


pvprop <- NULL
for (zzz in names(prop.rates)) {
  stat_mu = prop.rates[zzz]
  stat_sd = prop.rates.sd[zzz]
  zscore = stat_mu / stat_sd
  pvprop <- c(pvprop, 1 - pchisq(q = zscore**2, df = 1))
}
which(p.adjust(pvprop,method='BH') < 0.05)

pvprop_ind <- NULL
for (zzz in names(prop.rates)[-1]) {
  stat_mu = prop.rates[zzz] - prop.rates[1]
  stat_sd = sqrt(prop.rates.sd[zzz]**2 + prop.rates.sd[1]**2)
  zscore = stat_mu / stat_sd
  pvprop_ind <- c(pvprop_ind, 1 - pchisq(q = zscore**2, df = 1))
}
which(p.adjust(pvprop_ind,method='BH') < 0.05)


#################################################################

# Visualize

set.seed(123)
par(mfrow = c(1,2))
boxplot(rates, frame = F, outline = F, ylim = c(0,max(rates+1.96*rates.sd)),
        ylab = 'No. of clusters per generation', main = 'Clusters across genotypes')
newx <- jitter(rep(1, length(rates)), amount = 0.1)
points(x = newx, y = rates, col = 'gray', pch = 16) 
#o <- which(p.adjust(pvclust_ind,method='BH') < 0.1)
o <- 0
points(x = newx[o+1], y = rates[o+1], col = 'black', pch = 16) 
#arrows(x0 = newx[o+1],y0 = rates[o+1] - 1.96*rates.sd[o+1],y1=rates[o+1]+1.96*rates.sd[o+1],
#       col = 'gray21',lwd=0.5,length=0)
text(x = c(newx[o+1][1] + 0.3, newx[o+1][-1] + 0.6), y = rates[o+1] + 0.1, font = 3,
     labels = names(rates)[o+1], cex = 0.7)
abline(h = rates['N2.Rad'], lty = 2)
o <- which(p.adjust(pvclust,method='BH') < 0.1)
#arrows(x0 = newx[o],y0 = rates[o] - 1.96*rates.sd[o],y1=rates[o]+1.96*rates.sd[o],
#       col = 'gray21',lwd=0.5,length=0)

boxplot(prop.rates, frame = F, outline = F, ylim = c(0,max(prop.rates+1.96*prop.rates.sd)),
        ylab = 'Prop. of clustered mut-s', main = 'Proportion of clustered mutations\n across genotypes')
newx2 <- jitter(rep(1, length(prop.rates)), amount = 0.1)
points(x = newx2, y = prop.rates, col = 'gray', pch = 16) 
o <- 0
points(x = newx2[o+1], y = prop.rates[o+1], col = 'black', pch = 16) 
abline(h = prop.rates['N2.Rad'], lty = 2)
text(x = c(newx2[o+1][1] + 0.3, newx2[o+1][-1] + 0.6), y = prop.rates[o+1] + 0.02, font = 3,
     labels = names(prop.rates)[o+1], cex = 0.7)
o2 <- which(p.adjust(pvprop,method='BH') < 0.1)
#arrows(x0 = newx2[o2],y0 = prop.rates[o2] - 1.96*prop.rates.sd[o2],
#       y1=prop.rates[o2]+1.96*prop.rates.sd[o2],
#       col = 'gray21',lwd=0.5,length=0)
legend('topright', legend = c('IR mutation clustering in wildtype','significantly different\n from wildtype (FDR 10%)'), col = c('black','darkred'),
                              pch = 16, bty = 'n',border = NA, cex = 0.7)




# Distributions of cluster sizes and genomic span
number_per_cluster <- list()
cluster_sizes <- list()
for (gene in unique(allclust$genotype)) {
  
  distances <- diff(allclust$start[allclust$genotype == gene])
  rle(as.numeric(abs(distances) < 100)) -> dist_rle
  number_per_cluster[[gene]] <- dist_rle$lengths[dist_rle$values == 1] + 1
  cluster_sizes[[gene]] <- NULL
  cur_clust <- 0
  for (j in 1:length(distances)) {
    if (abs(distances[j]) < 100)
      cur_clust <- cur_clust + distances[j]
    else {
      if (cur_clust > 0)
        cluster_sizes[[gene]] <- c(cluster_sizes[[gene]], cur_clust)
      cur_clust <- 0
    }
  }
  if (length(distances) == 1 & cur_clust > 0) cluster_sizes[[gene]] <- c(cluster_sizes[[gene]], cur_clust)
  
}

number_per_cluster$`rad-54B` <- 0
number_per_cluster$`ndx-4` <- 0
number_per_cluster$`him-6` <- 0
names(number_per_cluster)[43] <- 'rad-54.B (gt3308)'
names(number_per_cluster)[52] <- 'rad-54.B (gk340656)'
names(number_per_cluster)[44] <- 'bub-3 (gt2000)'
names(number_per_cluster)[45] <- 'bub-3 (ok3437)'
names(number_per_cluster)[23] <- 'tdpo-1'

correct_order <- c(1,6,53,32,24,29,
                   36,11,37:38,14,26,35,
                   34,33,19,3,8,30,22,42,
                   43,52,20,21,51,2,15,12,28,
                   54,40,25,44:46,9:10)

cluster_sizes$`rad-54B` <- 0
cluster_sizes$`ndx-4` <- 0
cluster_sizes$`him-6` <- 0
names(cluster_sizes)[43] <- 'rad-54.B (gt3308)'
names(cluster_sizes)[52] <- 'rad-54.B (gk340656)'
names(cluster_sizes)[44] <- 'bub-3 (gt2000)'
names(cluster_sizes)[45] <- 'bub-3 (ok3437)'
names(cluster_sizes)[23] <- 'tdpo-1'

library(beeswarm)

# Mutations per cluster across genotypes
par(mar = c(10,4,4,2))
beeswarm(number_per_cluster[correct_order], las = 2, bty = 'n', method = 'hex', ylab = 'Number of mutations / cluster',
         ylim = c(0,10), cex = 0.5, pch = 16, corral = 'gutter', main = 'Sizes of IR-induced clusters across genotypes (muts)',
         xaxt = 'n')
axis(side = 1, at = c(1:length(correct_order)), labels = names(number_per_cluster)[correct_order], 
     font = 3, las = 2, lty = 0)

# Cluster span across genotypes
par(mar = c(10,4,4,2))
boxplot(cluster_sizes[correct_order], las = 2, frame = F, xaxt = 'n',
        main = 'Sizes of IR-induced clusters across genotypes (bps)',
        ylim = c(0,200), ylab = 'bp spanned by cluster', pch = 16, 
        cex = 0.5, outline = F)
axis(side = 1, at = c(1:length(correct_order)), labels = names(cluster_sizes)[correct_order], 
     font = 3, las = 2, lty = 0)
beeswarm(cluster_sizes[correct_order], las = 2, bty = 'n', method = 'hex', 
         ylab = 'bp spanned by cluster', col = 'gray56', add = T, 
         ylim = c(0,200), cex = 0.5, pch = 16, corral = 'gutter', 
         main = 'Sizes of IR-induced clusters across genotypes (bps)')
