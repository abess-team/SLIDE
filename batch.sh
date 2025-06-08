# Rscript simu_ws.R --type=8 --method="nodewise_logistic_gic2"
Rscript simu_ws.R --type=8 --method="RPLE_thres"
Rscript simu_ws.R --type=8 --method="RISE_thres"
Rscript simu_ws.R --type=8 --method="logRISE_thres"
Rscript simu_ws.R --type=8 --method="ELASSO_thres"

# Rscript simu_ws.R --type=10 --method="nodewise_logistic_gic2"
Rscript simu_ws.R --type=10 --method="RPLE_thres"
Rscript simu_ws.R --type=10 --method="RISE_thres"
Rscript simu_ws.R --type=10 --method="logRISE_thres"
Rscript simu_ws.R --type=10 --method="ELASSO_thres"

# Rscript simu_p.R --type=8 --method="nodewise_logistic_gic2"
Rscript simu_p.R --type=8 --method="RPLE_thres"
Rscript simu_p.R --type=8 --method="logRISE_thres"
Rscript simu_p.R --type=8 --method="RISE_thres"

# Rscript simu_p.R --type=10 --method="nodewise_logistic_gic2"
Rscript simu_p.R --type=10 --method="RPLE_thres"
Rscript simu_p.R --type=10 --method="logRISE_thres"
Rscript simu_p.R --type=10 --method="RISE_thres"

# Rscript simu_p.R --type=8 --method="ELASSO_thres"
# Rscript simu_p.R --type=10 --method="ELASSO_thres"

/usr/bin/shutdown 