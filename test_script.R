#Script for playing around

options(digits.secs = 3)

parquet_file <- '~/EAS_shared/cross_sleep/working/Data/Acc/spidermonkey/spidermonkey_1_Albus.parquet'

df <- read_parquet(parquet_file)

df_stand <- standardize_acc_to_uniform_sampling(df)

df_vedba <- acc_to_vedba(df_stand, rolling_mean_width = 20)
