# run the following in console to update dataset
# source ~/.bashrc
# conda activate
# python ~/code/casereveal/analysis/process_data.py

# whichtask <- "exculpatory_preponderous"
fpath <- paste0("~/projects/learnedjuror/task/respdata/data", whichtask)

catchdat <- fread(paste0(fpath, "_catches.csv"))
clickdat <- fread(paste0(fpath, "_answers.csv"))
scendat <- fread(paste0(fpath, "_realized_scen.csv"))
subjdat <- fread(paste0(fpath, "_subjectinfo.csv"))[,-1,with=F]
subjdat[,start:=as.POSIXct(start,format="%Y-%m-%d %H:%M:%S")]
# issue: those who completed no catch trials will not show up
setkey(subjdat,"start")
