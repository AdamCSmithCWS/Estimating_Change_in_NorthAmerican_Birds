#### jags function for SOCB model
library(ggplot2)
library(ggrepel)
library(ggforce)
library(rjags)
library(jagsUI)
library(stringr)
library(MCMCvis)
library(mgcv)

ind = "index"
lci = "lci"
uci = "uci"
year = "year"
base.yr = 1970
popsource = "Pop.source"

set.seed(2019)





popest = read.csv("Rosenberg et al species list.csv",
                   stringsAsFactors = F)
popest$popse = ((popest$popestuci-popest$popestlci)/(1.96*2))





indicesraw = read.csv("Rosenberg et al annual indices of abundance.csv",
                   stringsAsFactors = F)

indicesraw = indicesraw[order(indicesraw$species,indicesraw$year),]
# ###################### GAM smoothing of indices


sps = popest$species

for(ss in sps){
 wss = which(indicesraw$species == ss)
 tmp = indicesraw[wss,]
 
 fyr = unique(tmp$firstyear)
 lyr = unique(tmp$lastyear)
 
 torep = which(indicesraw$species == ss & (indicesraw$year >= fyr & indicesraw$year <= lyr))
 
 if(any(!is.na(tmp$lci.raw))){
   wprec = T
 }else{
   wprec = F
 }
 
 tmpd = tmp[which(!is.na(tmp$index.raw)),]
 

 
  nknots = min(c(11,max(floor(nrow(tmpd)/3),3)))
  if(ss %in% c("Cackling Goose","Greater White???fronted Goose")){nknots = 3} 
  
   
   
   form = as.formula(paste("lindex","~",
                           "s(year,k =",nknots,")"))
   
   


   ncounts = nrow(tmpd)
   ###### building the GAM basis function
   # gam basis functions created using function jagam from mgcv package          
   
   yminy = min(tmpd$year,na.rm = T)
   ymaxy = max(tmpd$year,na.rm = T)
   yearvec = tmpd$year-(yminy-1)
   yrs = seq(yminy,ymaxy,by = 1)
   nyears = length(yrs)
   ymin = min(yearvec)
   ymax = max(yearvec)
   tmpd$lindex = log(tmpd$index.raw)
   lindex = tmpd$lindex
   
   
   preddat = data.frame(lindex = 1,
                        year = yrs)
   
   
   
   if(wprec){
     ### if true then jags model runs 
     #### setting of the number of knots
     
     
   
   preci = 1/(((log(tmpd[,"uci.raw"])-log(tmpd[,"lci.raw"]))/(1.96*2))^2)


   
   gamprep = jagam(formula = form,
                   data = tmpd,
                   file = "tempgam.txt",
                   centred = T)
   
   gamprep.pred = jagam(formula = form,
                        data = preddat,
                        file = "tempgampred.txt",
                        centred = T)
   
   
   
   dat = list(X = gamprep$jags.data$X,
              S1 = gamprep$jags.data$S1,
              ncounts = nrow(tmpd),
              lindex = lindex,
              nknots = nknots,
              preci = preci,
              X.pred = gamprep.pred$jags.data$X,
              zero = gamprep$jags.data$zero)
   
   mg <- jags.model(data = dat,
                    file = paste0("GAM model smoothing indices jagam.txt"),
                    n.chains = 3)
   
   adaptest <- adapt(object = mg,
                     n.iter = 10)
   
   while(adaptest == F){
     adaptest <- adapt(object = mg,
                       n.iter = 1000)
     
   }
   
   nburn = 10000
   mgo = coda.samples(mg,
                      c("ind.pred"
                        #"rho"
                        #"mu",
                        #"b"
                      ),
                      n.iter = 10000,
                      n.burnin = nburn,
                      thin = 10)
   

   mgosum = summary(mgo)
   
   predg = mgosum$quantiles
   
   predgs = mgosum$statistics
   
   
   indicesraw[torep,"index"] <- exp(predg[grep(row.names(predg),pattern = "ind.pred"),"50%"])
   indicesraw[torep,"lci"] <- exp(predg[grep(row.names(predg),pattern = "ind.pred"),"2.5%"])
   indicesraw[torep,"uci"] <- exp(predg[grep(row.names(predg),pattern = "ind.pred"),"97.5%"])
   
   
    
 }else{ ### else if wprec
   
 
   m1 = gam(formula = form,
            data = tmpd)
   
   
   pred = predict(m1,newdata = preddat,
                  type = "link",
                  se.fit = T)
   

   
   
   indicesraw[torep,"index"] <- exp(pred$fit)
   indicesraw[torep,"lci"] <- exp(pred$fit-(1.96*pred$se.fit))
   indicesraw[torep,"uci"] <- exp(pred$fit+(1.96*pred$se.fit))
   
   
 } ### end if wprec
  
  
}


save(indicesraw,file = "output/post GAM indices.RDATA")
######################### end GAM smoothing of annual indices


load("output/post GAM indices.RDATA")


indices = indicesraw





yrs = sort(unique(indices$year))



indices[,"se"] <- ((indices[,uci]-indices[,lci])/(1.96*2))

indices <- indices[order(indices$firstyear,indices$species,indices$year),]

splist <- unique(indices[,c("species","firstyear","lastyear")])
splist3 = merge(splist,popest,by.x = "species",by.y = "species")

splist3$spfactor <- factor(splist3$species,
                           levels = splist3$species,
                           ordered = T)

splist3$spfact <- as.integer(splist3$spfactor) 

indices$spfactor = factor(indices$species,
                          levels = splist3$species,
                          ordered = T)






splist = splist3
base.i <- rep(NA, length = length(unique(indices$species)))
names(base.i) <- unique(indices$species)
base.se.sp <- rep(NA, length = length(unique(indices$species)))
names(base.se.sp) <- unique(indices$species)
se = "se"

for (sn in 1:nrow(splist)) {
  
  s = as.character(splist[sn,"species"])
  


  r <- which(indices[,"species"] == s)
  rmpre = which(indices[,"species"] == s &
                  is.na(indices[,ind]) &
                  indices[,"year"] < splist[sn,"firstyear"])
  rmpost = which(indices[,"species"] == s &
                   is.na(indices[,"se"]) &
                   indices[,"year"] > splist[sn,"lastyear"])
  
  
  base.s <- indices[which(indices$species == s & indices$year == base.yr),ind] #stores the base index value 
  base.se <- indices[which(indices$species == s & indices$year == base.yr),se] #stores the base se value 
  if(is.na(base.s)){
    byr <- min(indices[which(indices$species == s & !is.na(indices[,ind])),"year"],na.rm = T)
    base.s <- indices[which(indices$species == s & indices$year == byr),ind] #stores the base index value 
    base.se <- indices[which(indices$species == s & indices$year == byr),se] #stores the base se value 
    
  }
  for (y in r) {
    
    indices[y,"index.s"] <- (indices[y,ind])/base.s # standardized index WRT base year
    indices[y,"cvar.s"] <- (((indices[y,se]^2)/(indices[y,ind]^2))+((base.se^2)/(base.s^2))) 
    indices[y,"logthetahat"] <- log(indices[y,"index.s"])      
    indices[y,"prec.logthetahat"] <- 1/log(1+(indices[y,"cvar.s"]))
  }
  print(s)
}

indices = merge(indices,splist[,c("species","spfactor","spfact")],by = "species")
indices$year_i = (indices$year - base.yr)+1

indices = indices[order(indices$spfact,indices$year_i),]




splist = splist[order(splist$spfact),]
nspecies <- nrow(splist)
nyears = max(indices$year_i)

splist$g1 = as.integer(factor(splist$Winter.Biome))
splist$g2 = as.integer(factor(splist$Breeding.Biome))

grps = as.matrix(splist[,c("g1","g2")])

ngroups1 = max(grps[,1])
ngroups2 = max(grps[,2])

ident = function(x){
  return(x)
}
logthetahat <- as.matrix(tapply(indices[,"logthetahat"],indices[,c("spfact","year_i")],ident))
prec.logthetahat <- as.matrix(tapply(indices[,"prec.logthetahat"],indices[,c("spfact","year_i")],ident))

ne = splist$popest #population estimate
tau.ne = 1/splist$popse^2 #precision of population estimate








# 
# 
 yest = (splist$year_est1-base.yr)+1 # first year over which species population should be averaged
 
 yavg = 1+(splist$year_est2-splist$year_est1) ## number of years over which to average the species population estimate

 
 
 
 
 
##### species group indexing 
 
subgrps = unique(grps)
subgrps = subgrps[order(subgrps[,1],subgrps[,2]),]
nsubgrps = table(subgrps[,1])
subgrpsl = matrix(NA,ncol = ngroups1,nrow = max(nsubgrps))

nsppsubbiomesmat = table(grps[,1],grps[,2])

spsubgrpmat = array(NA,dim = c(ngroups1,ngroups2,max(nsppsubbiomesmat)))


for(i in 1:ngroups1){
  subs = subgrps[which(subgrps[,1] == i),2]
  
  subgrpsl[1:nsubgrps[i],i] = subs 
  for(j in subs){
    spsubgrpmat[i,j,1:nsppsubbiomesmat[i,j]] = which(grps[,1] == i & grps[,2] == j)
    
  }
}#i







#### breeding biome indexing

nsppbiomes = table(grps[,2])
spinbiomes = matrix(NA,nrow = max(nsppbiomes),ncol = max(ngroups2))
for(g in 1:ngroups2){
  
  spinbiomes[1:nsppbiomes[g],g] <- which(grps[,2] == g)
}




### wintering biome indexing
nsppwinters = table(grps[,1])
spinwinters = matrix(NA,nrow = max(nsppwinters),ncol = max(ngroups1))
for(g in 1:ngroups1){
  
  spinwinters[1:nsppwinters[g],g] <- which(grps[,1] == g)
}



### indexing for family summaries

splist$famfact = (factor(splist$Family))
splist$famfactn = as.integer(factor(splist$Family))

nsppfams = table(splist$famfact)
spinfams = matrix(NA,nrow = max(nsppfams),ncol = length(nsppfams))
for(f in 1:length(nsppfams)){
  fn = names(nsppfams)[f]
  spinfams[1:nsppfams[f],f] <- which(splist$Family == fn)
}
nfams = length(nsppfams)
fams = unique(splist[,c("Family","famfact","famfactn")])
fams = fams[order(fams$famfactn),]




#indexing for bird.group summaries

splist$birdgroupfact = (factor(splist$bird.group))
splist$birdgroupfactn = as.integer(factor(splist$bird.group))

nsppbirdgroups = table(splist$birdgroupfact)
spinbirdgroups = matrix(NA,nrow = max(nsppbirdgroups),ncol = length(nsppbirdgroups))
for(f in 1:length(nsppbirdgroups)){
  fn = names(nsppbirdgroups)[f]
  spinbirdgroups[1:nsppbirdgroups[f],f] <- which(splist$bird.group == fn)
}
nbirdgroups = length(nsppbirdgroups)
birdgroups = unique(splist[,c("bird.group","birdgroupfact","birdgroupfactn")])
birdgroups = birdgroups[order(birdgroups$birdgroupfactn),]



#indexing for migration summaries

splist$migratefact = (factor(splist$Migrate))
splist$migratefactn = as.integer(factor(splist$Migrate))

nsppmigrates = table(splist$migratefact)
spinmigrates = matrix(NA,nrow = max(nsppmigrates),ncol = length(nsppmigrates))
for(f in 1:length(nsppmigrates)){
  fn = names(nsppmigrates)[f]
  spinmigrates[1:nsppmigrates[f],f] <- which(splist$Migrate == fn)
}
nmigrates = length(nsppmigrates)
migrates = unique(splist[,c("Migrate","migratefact","migratefactn")])
migrates = migrates[order(migrates$migratefactn),]

#indexing for ai summaries

splist$aifact = (factor(splist$AI))
splist$aifactn = as.integer(factor(splist$AI))

nsppais = table(splist$aifact)
spinais = matrix(NA,nrow = max(nsppais),ncol = length(nsppais))
for(f in 1:length(nsppais)){
  fn = names(nsppais)[f]
  spinais[1:nsppais[f],f] <- which(splist$AI == fn)
}
nais = length(nsppais)
ais = unique(splist[,c("AI","aifact","aifactn")])
ais = ais[order(ais$aifactn),]

#indexing for native summaries

splist$nativefact = (factor(splist$native))
splist$nativefactn = as.integer(factor(splist$native))

nsppnatives = table(splist$nativefact)
spinnatives = matrix(NA,nrow = max(nsppnatives),ncol = length(nsppnatives))
for(f in 1:length(nsppnatives)){
  fn = names(nsppnatives)[f]
  spinnatives[1:nsppnatives[f],f] <- which(splist$native == fn)
}
nnatives = length(nsppnatives)
natives = unique(splist[,c("native","nativefact","nativefactn")])
natives = natives[order(natives$nativefactn),]










#imputing the missing data with assumptions of no-change from most recent year with data and gradually decreasing precision

wspecieslate = as.integer(splist[which(splist$firstyear > 1970),"spfact"])
yearswo1late = rep(1,nspecies)
yearswo2late = (splist[,"firstyear"])-1970


wspeciesearly = as.integer(splist[which(splist$lastyear < 2017),"spfact"])
yearswo1early = (splist[,"lastyear"]+1)-1969
yearswo2early = rep(nyears,nspecies)

yearsw1 = (splist[,"firstyear"])-1969
yearsw2 = (splist[,"lastyear"])-1969



prec.powdrop = 2 #exponential function for decreasing precision with years (i.e., precision decreases with square of the years since real data)

for( s in wspecieslate) { # wspecieslate = vector of species that don't have data in year-1
  for(y in yearswo1late[s]:yearswo2late[s]){

    logthetahat[s,y] <- logthetahat[s,yearsw1[s]]
    prec.logthetahat[s,y] <- prec.logthetahat[s,yearsw1[s]]/((yearsw1[s]-y)^prec.powdrop)# 
  }}


### imputing missing data for species without data at the end of the time series
## currently assumes that the precision decreases with the square of the number of years since the last data
for( s in wspeciesearly) { # wspecieslate = vector of species that don't have data in year-1
  for(y in yearswo1early[s]:yearswo2early[s]){

    logthetahat[s,y] <- logthetahat[s,yearsw2[s]]
    prec.logthetahat[s,y] <- prec.logthetahat[s,yearsw2[s]]/((y-yearsw2[s])^prec.powdrop)# 
    
    
    
  }}







  data.jags = list(nspecies = nspecies,
                   nyears = nyears,
                   logthetahat = logthetahat,
                   prec.logthetahat = prec.logthetahat,
                   ne = ne,
                   tau.ne = tau.ne,
                   yest = yest,
                   ngroups1 = ngroups1,
                   ngroups2 = ngroups2,
                   grps = grps,
                   yavg = yavg,
                   
                   spinbiomes = spinbiomes,
                   nsppbiomes = nsppbiomes,
                   
                   spinfams = spinfams,
                   nsppfams = nsppfams,
                   nfams = nfams,
                   
                   spinbirdgroups = spinbirdgroups,
                   nsppbirdgroups = nsppbirdgroups,
                   nbirdgroups = nbirdgroups,
                   
                   spinmigrates = spinmigrates,
                   nsppmigrates = nsppmigrates,
                   nmigrates = nmigrates,
                   
                   spinais = spinais,
                   nsppais = nsppais,
                   nais = nais,
                   
                   spinwinters = spinwinters,
                   nsppwinters = nsppwinters,
                   
                   spinnatives = spinnatives,
                   nsppnatives = nsppnatives,
                   nnatives = nnatives,
                   
                   subgrpsl = subgrpsl,
                   nsubgrps = nsubgrps,
                   nsppsubbiomesmat = nsppsubbiomesmat,
                   spsubgrpmat = spsubgrpmat)
  
  params = c("expmu1",
             #"expm2",
             "expmu2",
             "Psi",
             "tau",
             "NEst",
             "Nsum",
             "N",
             "Nlost",
             "Nalost",
             "Nlost.S",
             "Nlost.fam",
             "Nlost.biome",
             "Nalost.fam",
             "Nlost.migrate",
             "Nalost.migrate",
             "plost.migrate",
             "Nlost.birdgroup",
             "Nalost.birdgroup",
             "plost.birdgroup",
             "Nlost.ai",
             "Nalost.ai",
             "plost.ai",
             "plost.native",
             "Nlost.native",
             "Nalost.native",
             "plost.winter",
             "Nlost.winter",
             "Nalost.winter",
             "Nalost.biome",
             "plost.fam",
             "plost.biome",
             "plost",
             "Nsum.subgrp")
  
  mod = "Rosenberg et al model.txt"

  
  
  
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 30000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=10000           # Total number of steps to save per chain.
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps ) + burnInSteps) / nChains ) # Steps per chain.

t1 = Sys.time()


jagsMod = jags(data = data.jags,
               model.file = mod,
               n.chains = nChains,
               n.adapt = adaptSteps,
               n.burnin = burnInSteps,
               n.thin = thinSteps,
               n.iter = nIter+burnInSteps,
               parameters.to.save = params,
               parallel = T)

 #sumqalt = as.data.frame(jagsMod$summary)

q90 = function(x){
  quantile(x,probs = c(0.025,0.25,0.75,0.975))
}
sumq = MCMCsummary(jagsMod$samples,func = q90,Rhat = F,n.eff = F,func_name = c("lci","lqrt","uqrt","uci"))

sumq = data.frame(sumq)
names(sumq) <- c("mean","sd","lci95","med","uci95","lci","lqrt","uqrt","uci")
write.csv(sumq,paste0("output/population change parameters NA loss.csv"))
#write.csv(sumqalt,paste0("population change parameters NA loss w neff.csv"))


  
  faml = as.data.frame(sumq[paste0("Nlost.fam[",1:nfams,"]"),])
  famp = as.data.frame(sumq[paste0("plost.fam[",1:nfams,"]"),])
  names(famp) = paste0("plost.fam.",names(famp))
  faml = cbind(faml,famp)
  faml = cbind(faml,fams)
  

  
  
  nativel = as.data.frame(sumq[paste0("Nlost.native[",1:nnatives,"]"),])
  nativep = as.data.frame(sumq[paste0("plost.native[",1:nnatives,"]"),])
  names(nativep) = paste0("plost.native.",names(nativep))
  nativel = cbind(nativel,nativep)
  nativel = cbind(nativel,natives)
  

  migratel = as.data.frame(sumq[paste0("Nlost.migrate[",1:nmigrates,"]"),])
  migratep = as.data.frame(sumq[paste0("plost.migrate[",1:nmigrates,"]"),])
  names(migratep) = paste0("plost.migrate.",names(migratep))
  migratel = cbind(migratel,migratep)
  migratel = cbind(migratel,migrates)
  

  ail = as.data.frame(sumq[paste0("Nlost.ai[",1:nais,"]"),])
  aip = as.data.frame(sumq[paste0("plost.ai[",1:nais,"]"),])
  names(aip) = paste0("plost.ai.",names(aip))
  ail = cbind(ail,aip)
  ail = cbind(ail,ais)
  
 
  
  birdgroupl = as.data.frame(sumq[paste0("Nlost.birdgroup[",1:nbirdgroups,"]"),])
  birdgroupp = as.data.frame(sumq[paste0("plost.birdgroup[",1:nbirdgroups,"]"),])
  names(birdgroupp) = paste0("plost.birdgroup.",names(birdgroupp))
  birdgroupl = cbind(birdgroupl,birdgroupp)
  birdgroupl = cbind(birdgroupl,birdgroups)
  

  
  biomes = unique(splist[,c("g2","Breeding.Biome")])
  biomes = biomes[order(biomes$g2),]
  biomel = as.data.frame(sumq[paste0("Nlost.biome[",1:ngroups2,"]"),])
  biomep = as.data.frame(sumq[paste0("plost.biome[",1:ngroups2,"]"),])
  names(biomep) = paste0("plost.biome.",names(biomep))
  biomel = cbind(biomel,biomep)
  biomel = cbind(biomel,biomes)
  

  
  winters = unique(splist[,c("g1","Winter.Biome")])
  winters = winters[order(winters$g1),]
  winterl = as.data.frame(sumq[paste0("Nlost.winter[",1:ngroups1,"]"),])
  winterp = as.data.frame(sumq[paste0("plost.winter[",1:ngroups1,"]"),])
  names(winterp) = paste0("plost.winter.",names(winterp))
  winterl = cbind(winterl,winterp)
  winterl = cbind(winterl,winters)
  
 
  
  
  alll = data.frame(sumq[c("Nlost","plost"),])
  allp = data.frame(sumq[c("plost","Nlost"),])
  names(allp) = paste0("plost",names(allp))
  alll = cbind(alll,allp)
  alll = alll[1,]
  alll[,c(19:21)] <- NA
  biomel[,21] <- NA
  winterl[,21] <- NA
  names(alll)[20] = "Group"
  alll$nspecies = nspecies
  nativel$nspecies = nsppnatives
  migratel$nspecies = nsppmigrates
  birdgroupl$nspecies = nsppbirdgroups
  ail$nspecies = nsppais
  biomel$nspecies = nsppbiomes
  winterl$nspecies = nsppwinters
  faml$nspecies = nsppfams
  nms = names(alll)
  
  names(nativel) = nms
  names(migratel) = nms
  names(birdgroupl) = nms
  names(ail) = nms
  names(biomel) = nms
  names(winterl) = nms
  names(faml) = nms
  
  allsums = rbind(alll,
                  nativel,
                  migratel,
                  birdgroupl,
                  ail,
                  biomel,
                  winterl,
                  faml, stringsAsFactors = FALSE)
  
  allsout = allsums[,c("Group",
                       "nspecies",
                       "med",
                       "lci",
                       "uci",
                       "plostmed",
                       "plostlci",
                       "plostuci")]
  for(j in c("med","lci","uci")){
    allsout[,j] <- signif(allsout[,j]/1e6,5)
  }
  
  
  
  
  #### biome summaries
  
  biomepop = expand.grid(biome = 1:ngroups2,
                         yrs = 1:nyears)
  biomepop$param = paste0("Nalost.biome[",biomepop$biome,",",biomepop$yrs,"]")
  
  biomepopt = as.data.frame(sumq[biomepop$param,])
  biomepop = cbind(biomepop,biomepopt)
  biomepop = merge(biomepop,biomes,by.x = "biome",by.y = "g2")
  biomepop$year = biomepop$yrs + (base.yr-1)
  biomepop = biomepop[order(biomepop$Breeding.Biome,biomepop$year),]
  biomlab = biomepop[which(biomepop$year == 2017),]
  
  pmain = ggplot(data = biomepop,aes(x = year,y = med))+
    geom_ribbon(aes(x = year,ymin = lci,ymax = uci,group = Breeding.Biome,fill = Breeding.Biome),alpha = 0.2)+
    geom_line(aes(colour = Breeding.Biome))+
    geom_label_repel(data = biomlab,aes(label = Breeding.Biome,colour = Breeding.Biome),xlim = c(2017,2040))+
    labs(x = "",y = "Change in number of birds in North America")+
    xlim(1970,2035)+
    theme_minimal()+
    theme(legend.position = "none")
  
  
  pdf(paste0("Biome population population change trajectory  loss.pdf"))
  print(pmain)
  dev.off()
  
  pdf(paste0("Biome population population change trajectory facet loss.pdf"))
  
  for(jj in 1:ceiling(ngroups2/4)){
    pmain = ggplot(data = biomepop,aes(x = year,y = med))+
      
      geom_ribbon(aes(x = year,ymin = lci,ymax = uci),fill = grey(0.5),alpha = 0.2)+
      geom_line(aes(colour = "Breeding.Biome"))+
      facet_wrap_paginate(~Breeding.Biome,ncol = 2,nrow = 2,scales = "free",page = jj)
    print(pmain)
  }
  dev.off()
  write.csv(biomepop[,c("Breeding.Biome",
                        "year",
                        "med",
                        "lci",
                        "uci")],"Breeding Biome population change trajectories.csv")
  
  
  
  
  #### wintering summaries
  
  winterpop = expand.grid(winter = 1:ngroups2,
                          yrs = 1:nyears)
  winterpop$param = paste0("Nalost.winter[",winterpop$winter,",",winterpop$yrs,"]")
  
  winterpopt = as.data.frame(sumq[winterpop$param,])
  winterpop = cbind(winterpop,winterpopt)
  winterpop = merge(winterpop,winters,by.x = "winter",by.y = "g1")
  winterpop$year = winterpop$yrs + (base.yr-1)
  winterpop = winterpop[order(winterpop$Winter.Biome,winterpop$year),]
  biomlab = winterpop[which(winterpop$year == 2017),]
  
  pmain = ggplot(data = winterpop,aes(x = year,y = med))+
    geom_ribbon(aes(x = year,ymin = lci,ymax = uci,group = Winter.Biome,fill = Winter.Biome),alpha = 0.2)+
    geom_line(aes(colour = Winter.Biome))+
    geom_label_repel(data = biomlab,aes(label = Winter.Biome,colour = Winter.Biome),xlim = c(2017,2040))+
    labs(x = "",y = "Change in number of birds in North America")+
    xlim(1970,2035)+
    theme_minimal()+
    theme(legend.position = "none")
  
  
  pdf(paste0("winter population population change trajectory loss.pdf"))
  print(pmain)
  dev.off()
  
  pdf(paste0("winter population population change trajectory facet loss.pdf"))
  
  for(jj in 1:ceiling(ngroups2/4)){
    pmain = ggplot(data = winterpop,aes(x = year,y = med))+
      
      geom_ribbon(aes(x = year,ymin = lci,ymax = uci),fill = grey(0.5),alpha = 0.2)+
      geom_line(aes(colour = "Winter.Biome"))+
      facet_wrap_paginate(~Winter.Biome,ncol = 2,nrow = 2,scales = "free",page = jj)
    print(pmain)
  }
  dev.off()
  
  write.csv(winterpop[,c("Winter.Biome",
                         "year",
                         "med",
                         "lci",
                         "uci")],"wintering group population change trajectories.csv")
  
  
  #### fam summaries
  
  fampop = expand.grid(fam = 1:nfams,
                       yrs = 1:nyears)
  fampop$param = paste0("Nalost.fam[",fampop$fam,",",fampop$yrs,"]")
  
  fampopt = as.data.frame(sumq[fampop$param,])
  fampop = cbind(fampop,fampopt)
  fampop = merge(fampop,fams,by.x = "fam",by.y = "famfactn")
  fampop$year = fampop$yrs + (base.yr-1)
  biomlab = fampop[which(fampop$year == 2017),]
  
  
  
  
  
  
  
  pmain = ggplot(data = fampop,aes(x = year,y = med))+
    
    geom_ribbon(aes(x = year,ymin = lci,ymax = uci,fill = Family),alpha = 0.2)+
    geom_line(aes(colour = Family))+
    labs(x = "",y = "Change in number of birds in North America")+
  geom_label_repel(data = biomlab,aes(label = Family,colour = Family),xlim = c(2017,2050))+
    xlim(1970,2045)+
    theme_minimal()+
    theme(legend.position = "none")
  
  pdf(paste0("Family population population change trajectory loss.pdf"),
      width = 14,
      height = 10)
  print(pmain)
  dev.off()
  
  pdf(paste0("Family population population change trajectory facet loss.pdf"))
  
  for(jj in 1:ceiling(ngroups2/4)){
    pmain = ggplot(data = fampop,aes(x = year,y = med))+
      
      geom_ribbon(aes(x = year,ymin = lci,ymax = uci,fill = Family),alpha = 0.2)+
      geom_line(aes(colour = Family))+
      theme_minimal()+
      theme(legend.position = "none")+
      facet_wrap_paginate(~Family,ncol = 2,nrow = 2,scales = "fixed",page = jj)
    print(pmain)
  }
  dev.off()
  

#### overall summaries
  


totp = as.data.frame(sumq[paste0("Nsum[",1:nyears,"]"),c("med","lci","uci")])
names(totp) = paste0("N_",names(totp))
totp$year_i = 1:nyears
totp$year = totp$year_i + (base.yr-1)



lossa = sumq[paste0("Nalost[",1:nyears,"]"),c("med","lci","uci")]
names(lossa) = paste0("Loss_",names(lossa))
lossa = cbind(totp,lossa)

write.csv(lossa,row.names = F,"overall avifauna trajectories N and loss.csv")



lost = as.data.frame(sumq["Nlost",])


lost.s = as.data.frame(sumq[paste0("Nlost.S[",1:nspecies,"]"),])
names(lost.s) = paste0("Loss_",names(lost.s))
lost.s = cbind(splist,lost.s)
Psi = as.data.frame(sumq[paste0("Psi[",1:nspecies,"]"),])
names(Psi) = paste0("Psi_",names(Psi))
lost.s = cbind(lost.s,Psi)
NEst = as.data.frame(sumq[paste0("NEst[",1:nspecies,"]"),])
names(NEst) = paste0("NEst_",names(NEst))
lost.s = cbind(lost.s,NEst)




write.csv(lost.s,paste0("population change by species NA loss.csv"))


#summarize th eproportion of declining species for tables 1 and S2
lost.s$decline = F
lost.s[which(lost.s$Loss_med > 0),"decline"] = T
allsout$node = row.names(allsout)


for(j in 1:nrow(allsout)){
  if(j == 1){
    S = lost.s[,"decline"]
    tf = table(S)/length(S)
    
  }else{
    nn = gsub(gsub(allsout[j,"node"],pattern = "Nlost.",fixed = T,replacement = ""),pattern = "\\[.*",replacement = "")
    nn = paste0(nn,"fact")
    
    if(grepl(nn,pattern = "biome")){
      nn = "Breeding.Biome" 
    }
    if(grepl(nn,pattern = "winter")){
      nn = "Winter.Biome" 
    }
    S = lost.s[which(lost.s[,nn] == allsout[j,"Group"]),"decline"]
    tf = table(S)/length(S)
  }
  allsout[j,"proportion.species.decline"] = as.numeric(tf["TRUE"])
  allsout[j,"n.species.decline"] = sum(S)
}
write.csv(allsout,"table 1 and table S2.csv")



lostpy = 0.4*max(totp$N_med)
if(lost$med > 0){
  lostp = paste0(signif(lost$med/1e9,4)," Billion birds lost [",signif(lost$lci/1e9,4),"-",signif(lost$uci/1e9,4),"]")
}else{
  lostp = paste0(signif(abs(lost$med/1e9),4)," Billion birds gained [",signif(abs(lost$lci)/1e9,4),"-",signif(abs(lost$uci)/1e9,4),"]")
  
}
sppop = sumq[paste0("N[",rep(1:nspecies,each = nyears),",",rep(1:nyears,times = nspecies),"]"),]
sppop = as.data.frame(sppop)
sppop$spfact = rep(1:nspecies,each = nyears)
sppop$year_i = rep(1:nyears,times = nspecies)
sppop$year =  sppop$year_i + (base.yr-1)
sppop = merge(sppop,splist,
              by= "spfact")





  spstartpop = sppop[which(sppop$year == base.yr),]
  spstartpop$meanpopstart = spstartpop$mean
  indices2 = merge(indices,spstartpop[,c("species","meanpopstart")],by = "species")
  
  for(s in unique(indices2$species)){
    wsp = which(indices2$species == s)
    
    by = unique(indices2[wsp,"firstyear"])
    
    wspb = which(indices2$species == s & indices2$year == by)
    
    indices2[wsp,"index.s.raw"] <- indices2[wsp,"index.raw"]/indices2[wspb,"index.raw"]
    
  }
  
  indices2$rescindex = indices2$index.s*indices2$meanpopstart 
  indices2$rescindex.raw = indices2$index.s.raw*indices2$meanpopstart 
  sppop2 = merge(sppop,indices2[,c("species","rescindex","rescindex.raw","year")],
                 by = c("species","year")) 
  
  

for(j in c("N_med","N_lci","N_uci")){
  totp[,j] = totp[,j]/1e9
}

splabs = sppop2[which(sppop2$year == max(sppop2$year)),]
pmain = ggplot(data = totp,aes(x = year,y = N_med))+
  
  geom_ribbon(aes(x = year,ymin = N_lci,ymax = N_uci),fill = grey(0.5),alpha = 0.2)+
  geom_line()+
  labs(x = "",y = "Number of birds in North America (Billions)")+
  scale_y_continuous(limits = c(0,11),expand = expand_scale(mult = c(0, 0)))+
  annotate(geom = "text", x = 1990,y = lostpy/1e9,label = lostp)+
  theme_minimal()+
  theme(legend.position = "none")

pdf(paste0("total population population change trajectory.pdf"))
print(pmain)
dev.off()


for(j in c("Loss_med","Loss_lci","Loss_uci")){
  lossa[,j] = (lossa[,j]/1e9)
}

pmain = ggplot(data = lossa,aes(x = year,y = Loss_med))+
  
  geom_ribbon(aes(x = year,ymin = Loss_lci,ymax = Loss_uci),fill = grey(0.5),alpha = 0.2)+
  geom_line()+
  labs(x = "",y = "Change in number of birds in North America (Billions)")+
  annotate(geom = "text", x = 2000,y = -1,label = lostp)+
  theme_minimal()+
  theme(legend.position = "none")

pdf(paste0("overall loss trajectory.pdf"))
print(pmain)
dev.off()


lost.st = lost.s[rev(order(lost.s$Loss_med)),]
spord = unique(lost.st$species)
sppop2$spsort = factor(sppop2$species,levels = spord,ordered = T)
# 

sppop2 = sppop2[order(sppop2$spsort,sppop2$year),]

for(cl in c("lci","uci","med","lqrt","uqrt","rescindex","rescindex.raw")){
  sppop2[,cl] = sppop2[,cl]/1e6
}

spinlabs = sppop2[which(sppop2$year == 1970),]
spinlabs = merge(spinlabs,lost.s[,c("species","Loss_med","Loss_lci","Loss_uci","decline")],by = "species")
spinlabs = spinlabs[order(spinlabs$spsort),]

decs = which(spinlabs$decline)
gns = which(spinlabs$decline == F)
spinlabs[decs,"labs"] = paste0(signif(-1*(spinlabs[decs,"Loss_med"])/1e6,2),"M "," [",signif(-1*(spinlabs[decs,"Loss_uci"])/1e6,2),":",signif(-1*(spinlabs[decs,"Loss_lci"])/1e6,2),"]")

spinlabs[gns,"labs"] = paste0("+",signif(-1*(spinlabs[gns,"Loss_med"])/1e6,2),"M "," [",signif(-1*(spinlabs[gns,"Loss_uci"])/1e6,2),":",signif(-1*(spinlabs[gns,"Loss_lci"])/1e6,2),"]")

rwsnodat = which(is.na(sppop2$rescindex))

write.csv(sppop2,"output/species modeled trajectories w projections.csv")

sppop2[rwsnodat,c("lci","uci","lqrt","uqrt","med")] = NA 

write.csv(sppop2,"output/species modeled trajectories.csv")

pdf(paste0("individual populations manypage.pdf"))

for(jj in 1:ceiling(nspecies/9)){
  pmain = ggplot(data = sppop2,aes(x = year,y = med))+
    geom_ribbon(data = sppop2,aes(x = year,ymin = lci,ymax = uci),alpha = 0.2)+
    geom_ribbon(data = sppop2,aes(x = year,ymin = lqrt,ymax = uqrt),alpha = 0.2)+
    geom_line(data = sppop2,aes(x = year,y = med))+
    geom_line(data = sppop2,aes(x = year,y = rescindex,colour = Breeding.Biome))+
    geom_point(data = sppop2,aes(x = year,y = rescindex.raw,colour = Breeding.Biome),size = 0.9)+
    geom_text(data = spinlabs,aes(x = year,y = uci,label = labs),nudge_x = 20,size = 2)+
    geom_text(data = spinlabs,aes(x = year,y = lci,label = Breeding.Biome,colour = Breeding.Biome),nudge_x = 20,size = 2)+
    labs(x = "",y = "Millions of birds in North American population")+
    theme_minimal()+
    theme(legend.position = "none")+
    facet_wrap_paginate(~spsort,ncol = 3,nrow = 3,scales = "free_y",page = jj)
  print(pmain)
}


dev.off()

pdf(paste0("predicted change in species populations.pdf"))

for(jj in 1:ceiling(nspecies/9)){
  pmain = ggplot(data = sppop2,aes(x = year,y = med))+
    geom_ribbon(data = sppop2,aes(x = year,ymin = lci,ymax = uci),alpha = 0.2)+
    geom_ribbon(data = sppop2,aes(x = year,ymin = lqrt,ymax = uqrt),alpha = 0.2)+
    geom_line(data = sppop2,aes(x = year,y = med))+
    #geom_line(data = sppop2,aes(x = year,y = rescindex),colour = "red")+
    labs(x = "",y = "Millions of birds in North American population")+
    theme_minimal()+
    geom_text(data = spinlabs,aes(x = year,y = uci,label = labs),nudge_x = 20,size = 2)+
    theme(legend.position = "none")+
    facet_wrap_paginate(~spsort,ncol = 3,nrow = 3,scales = "free_y",page = jj)
  print(pmain)
}


dev.off()




## plot of the expmu values showing the group-level trajectories over time

groups1 = unique(splist[,c("g1","Winter.Biome")])
nspg1 = table(splist$g1)
groups1 = groups1[order(groups1$g1),]
groups1 = cbind(groups1,nspg1)
names(groups1)[4] = "nspecies"

expmu1gps = expand.grid(g1 = 1:ngroups1,
                        yr = 1:nyears)
expmu1gps$year = expmu1gps$yr + (base.yr-1)
expmu1gps$param = paste0("expmu1[",expmu1gps$g1,",",expmu1gps$yr,"]")
expmu1 = sumq[expmu1gps$param,]
expmu1 = cbind(expmu1,expmu1gps)
expmu1 = merge(expmu1,groups1,
               by = c("g1"))
expmu1 = expmu1[order(expmu1$g1,expmu1$year),]

pdf("Winter.Biome group trajectories.pdf")
pe1 = ggplot(data = expmu1,
             aes(x = year,
                 y = med))+
  geom_line()+
  facet_wrap(~Winter.Biome,scales = "free_y")+
  geom_ribbon(aes(ymin = lci,
                  ymax = uci),
              alpha = 0.2)

print(pe1)

dev.off()


groups12 = unique(splist[,c("g1","g2","Winter.Biome","Breeding.Biome")])
nspg12 = data.frame(table(splist[,c("g1","g2")]))
groups12 = groups12[order(groups12$g1,groups12$g2),]

groups12 = merge(groups12,nspg12,by = c("g1","g2"))
names(groups12)[5] = "nspecies"



expmu2gps = expand.grid(g1 = 1:ngroups1,
                        g2 = 1:ngroups2,
                        yr = 1:nyears)
expmu2gps$param = paste0("expmu2[",expmu2gps$g1,",",expmu2gps$g2,",",expmu2gps$yr,"]")

expmu2gps$year = expmu2gps$yr + (base.yr-1)

expmu2 = sumq[expmu2gps$param,]

expmu2 = cbind(expmu2,expmu2gps)

expmu2 = merge(expmu2,groups12,
               by = c("g1","g2"))

expmu2 = expmu2[order(expmu2$g1,expmu2$g2,expmu2$year),]
expmu2$group = paste0(expmu2$Winter.Biome,"-",expmu2$Breeding.Biome)



pdf("Combined biome group trajectories.pdf")

for(jj in 1:ceiling(nrow(groups12)/9)){
  pe12 = ggplot(data = expmu2,
                aes(x = year,
                    y = med))+
    geom_line()+
    facet_wrap_paginate(~group,scales = "free_y",ncol = 3,nrow = 3,page = jj)+
    geom_ribbon(aes(ymin = lci,
                    ymax = uci,fill = Winter.Biome),
                alpha = 0.2)+
    theme(legend.position = "none")
  
  
  print(pe12)
}

dev.off()




### plot the abundance trajectories for the subgroups

########subgroup abundance trajectories
wsubgrp = grep(row.names(sumq),pattern = "Nsum.subgrp",fixed = T)
subgpout = as.data.frame(sumq[wsubgrp,])
subgpout$param = row.names(subgpout)
newcol = c("g1","g2","yr")
for(i in 1:nrow(subgpout)){
  wrb = (str_locate_all(subgpout[i,"param"],pattern = "\\[|\\,|\\]") )
  
  for(j in 1:3){
    cl = newcol[j]
    subgpout[i,cl] <- str_sub(subgpout[i,"param"],start = wrb[[1]][j,1]+1,
                              end = wrb[[1]][j+1,2]-1)
  }
  
  
  
}

subgpout = merge(subgpout,groups12,by = c("g1","g2"),all = T)

subgpout$year = as.integer(subgpout$yr) + 1969 
subgpout = subgpout[order(subgpout$g1,subgpout$g2,subgpout$year),]  
subgpout$group = paste0(subgpout$Winter.Biome," - ",subgpout$Breeding.Biome)
subgpout$splab = paste(subgpout$nspecies,"species")
gplabs = subgpout[which(subgpout$year == 1970),]

for(i in 1:nrow(gplabs)){
  g = gplabs[i,"group"]
  delta = subgpout[which(subgpout$year == 2017 & subgpout$group == g),"med"]-subgpout[which(subgpout$year == 1970 & subgpout$group == g),"med"]
  if(delta > 0){
    lb = paste0(signif(delta/1e6,2),"M gained")
  }else{
    lb = paste0(signif((-1*delta)/1e6,2),"M lost")
  }
  
  gplabs[i,"total_change"] = paste(lb,gplabs[i,"splab"])
}


for(cl in c("med","lci","uci")){
  subgpout[,cl] = subgpout[,cl]/1e6
  gplabs[,cl] = gplabs[,cl]/1e6
}

pdf("Combined biome group abundance trajectories.pdf")

for(jj in 1:ceiling(nrow(groups12)/9)){
  pe12 = ggplot(data = subgpout,
                aes(x = year,
                    y = med))+
    geom_line(aes(colour = Breeding.Biome))+
    facet_wrap_paginate(~group,scales = "free_y",ncol = 3,nrow = 3,page = jj)+
    geom_ribbon(aes(ymin = lci,
                    ymax = uci,fill = Breeding.Biome),
                alpha = 0.2)+
    geom_label(data = gplabs,aes(label = total_change,x = year,y = uci),nudge_x = 25, size = 2.5, colour = grey(0.4))+
    labs(x = "",y = "Millions of birds in North American population")+
    theme_minimal()+
    theme(legend.position = "none",
          strip.text = element_text(size = 6))
  
  
  print(pe12)
}

dev.off()

write.csv(indices,"original data.csv")

#save.image("NA avifanual change.RData")










