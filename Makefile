all: Output/final_ranked_targets.csv

Data/Processed/expression_matrix.csv: Data/Input/raw_data.csv Scripts/01_extract_expression_data.R
	Rscript Scripts/01_extract_expression_data.R

Data/Processed/cutoffs.csv: Data/Processed/expression_matrix.csv Scripts/02_cal_cutoff_percentage_above_cutoff.R
	Rscript Scripts/02_cal_cutoff.R

Output/final_ranked_targets.csv: Data/Processed/cutoffs.csv Scripts/03_risk_benefit_pipeline.py
	python Scripts/03_risk_benefit_pipeline.py
