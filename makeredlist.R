alldbs <- c("exculpatory_pilot",
              "exculpatory_burdenofproof",
              "exculpatory_preponderous",
            "exculpatory_conditional",
            "exculpatory_rateless")

uids <- vector()
for (i in 1:length(alldbs)) {
  print(i)
  whichtask <- alldbs[i]
  source("process_data.R")
  uids <- c(uids, subjdat$uid)
}

uids <- uids[!str_detect(uids,"test")]
write(uids, "../task/redlist.txt")
