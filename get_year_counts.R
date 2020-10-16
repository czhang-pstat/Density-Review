# Generate point processes for density estimation of Mortality data
fls = list.files(paste(wd, "/Data/HMD/bltper_1x1", sep = "")) # get file names
pos1 = regexpr('\\.', fls)
pos2 = regexpr('NP', fls); pos2[pos2 == -1] = Inf
ctry = substr(fls, 1, min(pos1, pos2)-1) # get contry names

x = seq(0.5, 109.5, 1)
tms = list()
yrs = list()

for(i in 1:length(ctry)){
  tmp = read.table(fls[i], skip = 2, header = T)
  yrs[[i]] = unique(tmp$Year)
  tms[[i]] = sapply(1:length(yrs[[i]]), function(j) as.numeric(tmp$dx[tmp$Year == yrs[[i]][j]]))[-111,] # getting rid of counts of those who died above 110
}
names(tms) = names(yrs) = ctry

writeMat('./age_counts.mat', ctry = ctry, x = x, tms = tms, yrs = yrs)
