
# script to compute name distribution
# http://www.census.gov/topics/population/genealogy/data/2000_surnames.html
rm(list=ls(all=T))
setwd('~/Dropbox/Projects/Microclustering/Code/')

fn <- read.table('Data/ssa-names/yob1945.txt',sep=',')
names(fn) <- c('first','sex','freq')
for (yr in 1946:2015) {
  print(yr)
  fn0 <- read.table(paste('Data/ssa-names/yob',yr,'.txt',sep=''),sep=',')
  names(fn0) <- c('first','sex','freq')
  fn <- rbind(fn,fn0)
}
fnu <- tapply(fn$freq,fn$first,sum)
fnun <- rownames(fnu)
fnu <- data.frame(fnu)
names(fnu) <- 'freq.f'
fnu$first <- fnun

ln <- read.table('Data/names/app_c.csv',sep=',',header=T)
lnu <- tapply(ln$count,ln$name,sum)
lnun <- rownames(lnu)
lnu <- data.frame(lnu)
names(lnu) <- 'freq.l'
lnu$last <- lnun

write.csv(fnu,'Data/names/first.csv')
write.csv(lnu,'Data/names/last.csv')









