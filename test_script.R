#Script for playing around

options(digits.secs = 3)

parquet_file <- '~/EAS_shared/cross_sleep/working/Data/Acc/spidermonkey/spidermonkey_1_Albus.parquet'

df <- read_parquet(parquet_file)

df_stand <- standardize_acc_to_uniform_sampling(df)

df_vedba <- acc_to_vedba(df_stand, rolling_mean_width = 20)

df_vedba_sub <- df_vedba[seq(1,nrow(df_vedba),20),]

#ßwrite_parquet(x = df_stand, sink = '~/EAS_shared/cross_sleep/working/Data/Acc_standard/spidermonkey_1_Albus_vedba_standard.parquet')
write_parquet(x = df_vedba, sink = '~/EAS_shared/cross_sleep/working/Data/VeDBA_cont/spidermonkey_1_Albus_vedba_standard.parquet')
