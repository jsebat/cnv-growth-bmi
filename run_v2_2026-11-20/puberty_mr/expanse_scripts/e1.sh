cat /expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v9e_chop_genotype_bands/run_v9e.sbatch
M=/expanse/projects/sebat1/msacks/UKBB/comparitiveBodySize
for f in field1787_maternal_smoking_around_birth.tsv field1747_hair_colour.tsv field1697_comparative_height_size_age10.tsv field20022_1777_birthweight_multiplebirth.tsv; do head -1 $M/$f; done
R=/expanse/projects/sebat1/a1sriniv/pgs_sumstats_organized/raw_sumstats
zcat $R/GCST90002409_childhood_BMI.tsv.gz | wc -l
squeue -u jsebat | head -20
