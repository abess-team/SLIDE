Rscript simu_degree.R --method="RPLE_thres" --case=0
Rscript simu_degree.R --method="ELASSO_thres" --case=0
Rscript simu_degree.R --method="SLIDE" --case=0
Rscript simu_degree.R --method="RPLE_thres" --case=1
Rscript simu_degree.R --method="ELASSO_thres" --case=1
Rscript simu_degree.R --method="SLIDE" --case=1

Rscript simu_beta.R --method="RPLE_thres" --type=3
Rscript simu_beta.R --method="RISE_thres" --type=3
Rscript simu_beta.R --method="logRISE_thres" --type=3
Rscript simu_beta.R --method="ELASSO_thres" --type=3
Rscript simu_beta.R --method="SLIDE" --type=3
Rscript simu_beta.R --method="RPLE_thres" --type=1
Rscript simu_beta.R --method="RISE_thres" --type=1
Rscript simu_beta.R --method="logRISE_thres" --type=1
Rscript simu_beta.R --method="ELASSO_thres" --type=1
Rscript simu_beta.R --method="SLIDE" --type=1
Rscript simu_beta.R --method="RPLE_thres" --type=4
Rscript simu_beta.R --method="RISE_thres" --type=4
Rscript simu_beta.R --method="logRISE_thres" --type=4
Rscript simu_beta.R --method="ELASSO_thres" --type=4
Rscript simu_beta.R --method="SLIDE" --type=4
Rscript simu_beta.R --method="RPLE_thres" --type=2
Rscript simu_beta.R --method="RISE_thres" --type=2
Rscript simu_beta.R --method="logRISE_thres" --type=2
Rscript simu_beta.R --method="ELASSO_thres" --type=2
Rscript simu_beta.R --method="SLIDE" --type=2
Rscript simu_beta.R --method="RPLE_thres" --type=5
Rscript simu_beta.R --method="RISE_thres" --type=5
Rscript simu_beta.R --method="logRISE_thres" --type=5
Rscript simu_beta.R --method="ELASSO_thres" --type=5
Rscript simu_beta.R --method="SLIDE" --type=5

Rscript simu_ws.R --method="RPLE_thres" --type=1
Rscript simu_ws.R --method="RISE_thres" --type=1
Rscript simu_ws.R --method="logRISE_thres" --type=1
Rscript simu_ws.R --method="ELASSO_thres" --type=1
Rscript simu_ws.R --method="SLIDE" --type=1
Rscript simu_ws.R --method="RPLE_thres" --type=3
Rscript simu_ws.R --method="RISE_thres" --type=3
Rscript simu_ws.R --method="logRISE_thres" --type=3
Rscript simu_ws.R --method="ELASSO_thres" --type=3
Rscript simu_ws.R --method="SLIDE" --type=3

Rscript simu_p.R --method="RPLE_thres" --type=1
Rscript simu_p.R --method="logRISE_thres" --type=1
Rscript simu_p.R --method="RISE_thres" --type=1
Rscript simu_p.R --method="ELASSO_thres" --type=1
Rscript simu_p.R --method="SLIDE" --type=1
Rscript simu_p.R --method="RPLE_thres" --type=3
Rscript simu_p.R --method="logRISE_thres" --type=3
Rscript simu_p.R --method="RISE_thres" --type=3
Rscript simu_p.R --method="ELASSO_thres" --type=3
Rscript simu_p.R --method="SLIDE" --type=3

Rscript simu_high.R --method="RPLE_cv_thres"
Rscript simu_high.R --method="RISE_cv_thres"
Rscript simu_high.R --method="logRISE_cv_thres"
Rscript simu_high.R --method="LogRelax"
Rscript simu_high.R --method="ELASSO_thres"
Rscript simu_high.R --method="SLIDE"

/usr/bin/shutdown 