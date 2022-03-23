
jid1=$(sbatch Script10_mgs_bp_hb/Correlation_gmm_metabolites.txt|cut -c21-)

sbatch --dependency=afterok:$jid1 Script10_mgs_bp_hb/Enrichment_gmm_subclass.txt