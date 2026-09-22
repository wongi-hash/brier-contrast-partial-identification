# Export the public Rotterdam and GBSG datasets from R survival.
# Verified locally with survival 3.8-3. Later package versions are accepted only
# if the normalized exported-data hashes pass stage 01's locked checks.
suppressPackageStartupMessages(library(survival))
rotterdam <- survival::rotterdam
gbsg <- survival::gbsg
rotterdam <- data.frame(rownames=rownames(rotterdam), rotterdam, check.names=FALSE)
gbsg <- data.frame(rownames=rownames(gbsg), gbsg, check.names=FALSE)
write.csv(rotterdam, "rotterdam.csv", row.names=FALSE, quote=FALSE, na="")
write.csv(gbsg, "gbsg.csv", row.names=FALSE, quote=FALSE, na="")
writeLines(capture.output(sessionInfo()), "R_SESSION_INFO.txt")
writeLines(capture.output(citation("survival")), "SURVIVAL_CITATION.txt")
cat("survival", as.character(packageVersion("survival")), "\n")
cat("rotterdam", nrow(rotterdam), "gbsg", nrow(gbsg), "\n")
