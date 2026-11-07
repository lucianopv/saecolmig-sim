# Per-fit seeds (R/fit_seed.R).

# Seeds of both configurations for every method, outcome and three replicates,
# each one the formula 1e7 * config + 1e6 * method + 1e4 * outcome + rep with the
# indexes written out: geo_bld carries the index of geo (2) and geo_blk_bld that
# of geo_blk (9); methods fh .. megb_sub are 1 .. 8; outcomes are 1 .. 6.
pinned_fit_seeds <- function() {
  utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
  config      method         outcome                    rep seed
  geo_bld     fh             migrated_last_12m            1 21010001
  geo_bld     fh             migrated_last_12m           57 21010057
  geo_bld     fh             migrated_last_12m          100 21010100
  geo_bld     fh             migrated_internal_12m        1 21020001
  geo_bld     fh             migrated_internal_12m       57 21020057
  geo_bld     fh             migrated_internal_12m      100 21020100
  geo_bld     fh             migrated_intl_12m            1 21030001
  geo_bld     fh             migrated_intl_12m           57 21030057
  geo_bld     fh             migrated_intl_12m          100 21030100
  geo_bld     fh             pers_migrated_12m            1 21040001
  geo_bld     fh             pers_migrated_12m           57 21040057
  geo_bld     fh             pers_migrated_12m          100 21040100
  geo_bld     fh             pers_migrated_internal_12m   1 21050001
  geo_bld     fh             pers_migrated_internal_12m  57 21050057
  geo_bld     fh             pers_migrated_internal_12m 100 21050100
  geo_bld     fh             pers_migrated_intl_12m       1 21060001
  geo_bld     fh             pers_migrated_intl_12m      57 21060057
  geo_bld     fh             pers_migrated_intl_12m     100 21060100
  geo_bld     fh_arcsin      migrated_last_12m            1 22010001
  geo_bld     fh_arcsin      migrated_last_12m           57 22010057
  geo_bld     fh_arcsin      migrated_last_12m          100 22010100
  geo_bld     fh_arcsin      migrated_internal_12m        1 22020001
  geo_bld     fh_arcsin      migrated_internal_12m       57 22020057
  geo_bld     fh_arcsin      migrated_internal_12m      100 22020100
  geo_bld     fh_arcsin      migrated_intl_12m            1 22030001
  geo_bld     fh_arcsin      migrated_intl_12m           57 22030057
  geo_bld     fh_arcsin      migrated_intl_12m          100 22030100
  geo_bld     fh_arcsin      pers_migrated_12m            1 22040001
  geo_bld     fh_arcsin      pers_migrated_12m           57 22040057
  geo_bld     fh_arcsin      pers_migrated_12m          100 22040100
  geo_bld     fh_arcsin      pers_migrated_internal_12m   1 22050001
  geo_bld     fh_arcsin      pers_migrated_internal_12m  57 22050057
  geo_bld     fh_arcsin      pers_migrated_internal_12m 100 22050100
  geo_bld     fh_arcsin      pers_migrated_intl_12m       1 22060001
  geo_bld     fh_arcsin      pers_migrated_intl_12m      57 22060057
  geo_bld     fh_arcsin      pers_migrated_intl_12m     100 22060100
  geo_bld     fh_spatial     migrated_last_12m            1 23010001
  geo_bld     fh_spatial     migrated_last_12m           57 23010057
  geo_bld     fh_spatial     migrated_last_12m          100 23010100
  geo_bld     fh_spatial     migrated_internal_12m        1 23020001
  geo_bld     fh_spatial     migrated_internal_12m       57 23020057
  geo_bld     fh_spatial     migrated_internal_12m      100 23020100
  geo_bld     fh_spatial     migrated_intl_12m            1 23030001
  geo_bld     fh_spatial     migrated_intl_12m           57 23030057
  geo_bld     fh_spatial     migrated_intl_12m          100 23030100
  geo_bld     fh_spatial     pers_migrated_12m            1 23040001
  geo_bld     fh_spatial     pers_migrated_12m           57 23040057
  geo_bld     fh_spatial     pers_migrated_12m          100 23040100
  geo_bld     fh_spatial     pers_migrated_internal_12m   1 23050001
  geo_bld     fh_spatial     pers_migrated_internal_12m  57 23050057
  geo_bld     fh_spatial     pers_migrated_internal_12m 100 23050100
  geo_bld     fh_spatial     pers_migrated_intl_12m       1 23060001
  geo_bld     fh_spatial     pers_migrated_intl_12m      57 23060057
  geo_bld     fh_spatial     pers_migrated_intl_12m     100 23060100
  geo_bld     fh_robust      migrated_last_12m            1 24010001
  geo_bld     fh_robust      migrated_last_12m           57 24010057
  geo_bld     fh_robust      migrated_last_12m          100 24010100
  geo_bld     fh_robust      migrated_internal_12m        1 24020001
  geo_bld     fh_robust      migrated_internal_12m       57 24020057
  geo_bld     fh_robust      migrated_internal_12m      100 24020100
  geo_bld     fh_robust      migrated_intl_12m            1 24030001
  geo_bld     fh_robust      migrated_intl_12m           57 24030057
  geo_bld     fh_robust      migrated_intl_12m          100 24030100
  geo_bld     fh_robust      pers_migrated_12m            1 24040001
  geo_bld     fh_robust      pers_migrated_12m           57 24040057
  geo_bld     fh_robust      pers_migrated_12m          100 24040100
  geo_bld     fh_robust      pers_migrated_internal_12m   1 24050001
  geo_bld     fh_robust      pers_migrated_internal_12m  57 24050057
  geo_bld     fh_robust      pers_migrated_internal_12m 100 24050100
  geo_bld     fh_robust      pers_migrated_intl_12m       1 24060001
  geo_bld     fh_robust      pers_migrated_intl_12m      57 24060057
  geo_bld     fh_robust      pers_migrated_intl_12m     100 24060100
  geo_bld     bayes_fh_beta  migrated_last_12m            1 25010001
  geo_bld     bayes_fh_beta  migrated_last_12m           57 25010057
  geo_bld     bayes_fh_beta  migrated_last_12m          100 25010100
  geo_bld     bayes_fh_beta  migrated_internal_12m        1 25020001
  geo_bld     bayes_fh_beta  migrated_internal_12m       57 25020057
  geo_bld     bayes_fh_beta  migrated_internal_12m      100 25020100
  geo_bld     bayes_fh_beta  migrated_intl_12m            1 25030001
  geo_bld     bayes_fh_beta  migrated_intl_12m           57 25030057
  geo_bld     bayes_fh_beta  migrated_intl_12m          100 25030100
  geo_bld     bayes_fh_beta  pers_migrated_12m            1 25040001
  geo_bld     bayes_fh_beta  pers_migrated_12m           57 25040057
  geo_bld     bayes_fh_beta  pers_migrated_12m          100 25040100
  geo_bld     bayes_fh_beta  pers_migrated_internal_12m   1 25050001
  geo_bld     bayes_fh_beta  pers_migrated_internal_12m  57 25050057
  geo_bld     bayes_fh_beta  pers_migrated_internal_12m 100 25050100
  geo_bld     bayes_fh_beta  pers_migrated_intl_12m       1 25060001
  geo_bld     bayes_fh_beta  pers_migrated_intl_12m      57 25060057
  geo_bld     bayes_fh_beta  pers_migrated_intl_12m     100 25060100
  geo_bld     bayes_fh_heavy migrated_last_12m            1 26010001
  geo_bld     bayes_fh_heavy migrated_last_12m           57 26010057
  geo_bld     bayes_fh_heavy migrated_last_12m          100 26010100
  geo_bld     bayes_fh_heavy migrated_internal_12m        1 26020001
  geo_bld     bayes_fh_heavy migrated_internal_12m       57 26020057
  geo_bld     bayes_fh_heavy migrated_internal_12m      100 26020100
  geo_bld     bayes_fh_heavy migrated_intl_12m            1 26030001
  geo_bld     bayes_fh_heavy migrated_intl_12m           57 26030057
  geo_bld     bayes_fh_heavy migrated_intl_12m          100 26030100
  geo_bld     bayes_fh_heavy pers_migrated_12m            1 26040001
  geo_bld     bayes_fh_heavy pers_migrated_12m           57 26040057
  geo_bld     bayes_fh_heavy pers_migrated_12m          100 26040100
  geo_bld     bayes_fh_heavy pers_migrated_internal_12m   1 26050001
  geo_bld     bayes_fh_heavy pers_migrated_internal_12m  57 26050057
  geo_bld     bayes_fh_heavy pers_migrated_internal_12m 100 26050100
  geo_bld     bayes_fh_heavy pers_migrated_intl_12m       1 26060001
  geo_bld     bayes_fh_heavy pers_migrated_intl_12m      57 26060057
  geo_bld     bayes_fh_heavy pers_migrated_intl_12m     100 26060100
  geo_bld     xgb_sub        migrated_last_12m            1 27010001
  geo_bld     xgb_sub        migrated_last_12m           57 27010057
  geo_bld     xgb_sub        migrated_last_12m          100 27010100
  geo_bld     xgb_sub        migrated_internal_12m        1 27020001
  geo_bld     xgb_sub        migrated_internal_12m       57 27020057
  geo_bld     xgb_sub        migrated_internal_12m      100 27020100
  geo_bld     xgb_sub        migrated_intl_12m            1 27030001
  geo_bld     xgb_sub        migrated_intl_12m           57 27030057
  geo_bld     xgb_sub        migrated_intl_12m          100 27030100
  geo_bld     xgb_sub        pers_migrated_12m            1 27040001
  geo_bld     xgb_sub        pers_migrated_12m           57 27040057
  geo_bld     xgb_sub        pers_migrated_12m          100 27040100
  geo_bld     xgb_sub        pers_migrated_internal_12m   1 27050001
  geo_bld     xgb_sub        pers_migrated_internal_12m  57 27050057
  geo_bld     xgb_sub        pers_migrated_internal_12m 100 27050100
  geo_bld     xgb_sub        pers_migrated_intl_12m       1 27060001
  geo_bld     xgb_sub        pers_migrated_intl_12m      57 27060057
  geo_bld     xgb_sub        pers_migrated_intl_12m     100 27060100
  geo_bld     megb_sub       migrated_last_12m            1 28010001
  geo_bld     megb_sub       migrated_last_12m           57 28010057
  geo_bld     megb_sub       migrated_last_12m          100 28010100
  geo_bld     megb_sub       migrated_internal_12m        1 28020001
  geo_bld     megb_sub       migrated_internal_12m       57 28020057
  geo_bld     megb_sub       migrated_internal_12m      100 28020100
  geo_bld     megb_sub       migrated_intl_12m            1 28030001
  geo_bld     megb_sub       migrated_intl_12m           57 28030057
  geo_bld     megb_sub       migrated_intl_12m          100 28030100
  geo_bld     megb_sub       pers_migrated_12m            1 28040001
  geo_bld     megb_sub       pers_migrated_12m           57 28040057
  geo_bld     megb_sub       pers_migrated_12m          100 28040100
  geo_bld     megb_sub       pers_migrated_internal_12m   1 28050001
  geo_bld     megb_sub       pers_migrated_internal_12m  57 28050057
  geo_bld     megb_sub       pers_migrated_internal_12m 100 28050100
  geo_bld     megb_sub       pers_migrated_intl_12m       1 28060001
  geo_bld     megb_sub       pers_migrated_intl_12m      57 28060057
  geo_bld     megb_sub       pers_migrated_intl_12m     100 28060100
  geo_blk_bld fh             migrated_last_12m            1 91010001
  geo_blk_bld fh             migrated_last_12m           57 91010057
  geo_blk_bld fh             migrated_last_12m          100 91010100
  geo_blk_bld fh             migrated_internal_12m        1 91020001
  geo_blk_bld fh             migrated_internal_12m       57 91020057
  geo_blk_bld fh             migrated_internal_12m      100 91020100
  geo_blk_bld fh             migrated_intl_12m            1 91030001
  geo_blk_bld fh             migrated_intl_12m           57 91030057
  geo_blk_bld fh             migrated_intl_12m          100 91030100
  geo_blk_bld fh             pers_migrated_12m            1 91040001
  geo_blk_bld fh             pers_migrated_12m           57 91040057
  geo_blk_bld fh             pers_migrated_12m          100 91040100
  geo_blk_bld fh             pers_migrated_internal_12m   1 91050001
  geo_blk_bld fh             pers_migrated_internal_12m  57 91050057
  geo_blk_bld fh             pers_migrated_internal_12m 100 91050100
  geo_blk_bld fh             pers_migrated_intl_12m       1 91060001
  geo_blk_bld fh             pers_migrated_intl_12m      57 91060057
  geo_blk_bld fh             pers_migrated_intl_12m     100 91060100
  geo_blk_bld fh_arcsin      migrated_last_12m            1 92010001
  geo_blk_bld fh_arcsin      migrated_last_12m           57 92010057
  geo_blk_bld fh_arcsin      migrated_last_12m          100 92010100
  geo_blk_bld fh_arcsin      migrated_internal_12m        1 92020001
  geo_blk_bld fh_arcsin      migrated_internal_12m       57 92020057
  geo_blk_bld fh_arcsin      migrated_internal_12m      100 92020100
  geo_blk_bld fh_arcsin      migrated_intl_12m            1 92030001
  geo_blk_bld fh_arcsin      migrated_intl_12m           57 92030057
  geo_blk_bld fh_arcsin      migrated_intl_12m          100 92030100
  geo_blk_bld fh_arcsin      pers_migrated_12m            1 92040001
  geo_blk_bld fh_arcsin      pers_migrated_12m           57 92040057
  geo_blk_bld fh_arcsin      pers_migrated_12m          100 92040100
  geo_blk_bld fh_arcsin      pers_migrated_internal_12m   1 92050001
  geo_blk_bld fh_arcsin      pers_migrated_internal_12m  57 92050057
  geo_blk_bld fh_arcsin      pers_migrated_internal_12m 100 92050100
  geo_blk_bld fh_arcsin      pers_migrated_intl_12m       1 92060001
  geo_blk_bld fh_arcsin      pers_migrated_intl_12m      57 92060057
  geo_blk_bld fh_arcsin      pers_migrated_intl_12m     100 92060100
  geo_blk_bld fh_spatial     migrated_last_12m            1 93010001
  geo_blk_bld fh_spatial     migrated_last_12m           57 93010057
  geo_blk_bld fh_spatial     migrated_last_12m          100 93010100
  geo_blk_bld fh_spatial     migrated_internal_12m        1 93020001
  geo_blk_bld fh_spatial     migrated_internal_12m       57 93020057
  geo_blk_bld fh_spatial     migrated_internal_12m      100 93020100
  geo_blk_bld fh_spatial     migrated_intl_12m            1 93030001
  geo_blk_bld fh_spatial     migrated_intl_12m           57 93030057
  geo_blk_bld fh_spatial     migrated_intl_12m          100 93030100
  geo_blk_bld fh_spatial     pers_migrated_12m            1 93040001
  geo_blk_bld fh_spatial     pers_migrated_12m           57 93040057
  geo_blk_bld fh_spatial     pers_migrated_12m          100 93040100
  geo_blk_bld fh_spatial     pers_migrated_internal_12m   1 93050001
  geo_blk_bld fh_spatial     pers_migrated_internal_12m  57 93050057
  geo_blk_bld fh_spatial     pers_migrated_internal_12m 100 93050100
  geo_blk_bld fh_spatial     pers_migrated_intl_12m       1 93060001
  geo_blk_bld fh_spatial     pers_migrated_intl_12m      57 93060057
  geo_blk_bld fh_spatial     pers_migrated_intl_12m     100 93060100
  geo_blk_bld fh_robust      migrated_last_12m            1 94010001
  geo_blk_bld fh_robust      migrated_last_12m           57 94010057
  geo_blk_bld fh_robust      migrated_last_12m          100 94010100
  geo_blk_bld fh_robust      migrated_internal_12m        1 94020001
  geo_blk_bld fh_robust      migrated_internal_12m       57 94020057
  geo_blk_bld fh_robust      migrated_internal_12m      100 94020100
  geo_blk_bld fh_robust      migrated_intl_12m            1 94030001
  geo_blk_bld fh_robust      migrated_intl_12m           57 94030057
  geo_blk_bld fh_robust      migrated_intl_12m          100 94030100
  geo_blk_bld fh_robust      pers_migrated_12m            1 94040001
  geo_blk_bld fh_robust      pers_migrated_12m           57 94040057
  geo_blk_bld fh_robust      pers_migrated_12m          100 94040100
  geo_blk_bld fh_robust      pers_migrated_internal_12m   1 94050001
  geo_blk_bld fh_robust      pers_migrated_internal_12m  57 94050057
  geo_blk_bld fh_robust      pers_migrated_internal_12m 100 94050100
  geo_blk_bld fh_robust      pers_migrated_intl_12m       1 94060001
  geo_blk_bld fh_robust      pers_migrated_intl_12m      57 94060057
  geo_blk_bld fh_robust      pers_migrated_intl_12m     100 94060100
  geo_blk_bld bayes_fh_beta  migrated_last_12m            1 95010001
  geo_blk_bld bayes_fh_beta  migrated_last_12m           57 95010057
  geo_blk_bld bayes_fh_beta  migrated_last_12m          100 95010100
  geo_blk_bld bayes_fh_beta  migrated_internal_12m        1 95020001
  geo_blk_bld bayes_fh_beta  migrated_internal_12m       57 95020057
  geo_blk_bld bayes_fh_beta  migrated_internal_12m      100 95020100
  geo_blk_bld bayes_fh_beta  migrated_intl_12m            1 95030001
  geo_blk_bld bayes_fh_beta  migrated_intl_12m           57 95030057
  geo_blk_bld bayes_fh_beta  migrated_intl_12m          100 95030100
  geo_blk_bld bayes_fh_beta  pers_migrated_12m            1 95040001
  geo_blk_bld bayes_fh_beta  pers_migrated_12m           57 95040057
  geo_blk_bld bayes_fh_beta  pers_migrated_12m          100 95040100
  geo_blk_bld bayes_fh_beta  pers_migrated_internal_12m   1 95050001
  geo_blk_bld bayes_fh_beta  pers_migrated_internal_12m  57 95050057
  geo_blk_bld bayes_fh_beta  pers_migrated_internal_12m 100 95050100
  geo_blk_bld bayes_fh_beta  pers_migrated_intl_12m       1 95060001
  geo_blk_bld bayes_fh_beta  pers_migrated_intl_12m      57 95060057
  geo_blk_bld bayes_fh_beta  pers_migrated_intl_12m     100 95060100
  geo_blk_bld bayes_fh_heavy migrated_last_12m            1 96010001
  geo_blk_bld bayes_fh_heavy migrated_last_12m           57 96010057
  geo_blk_bld bayes_fh_heavy migrated_last_12m          100 96010100
  geo_blk_bld bayes_fh_heavy migrated_internal_12m        1 96020001
  geo_blk_bld bayes_fh_heavy migrated_internal_12m       57 96020057
  geo_blk_bld bayes_fh_heavy migrated_internal_12m      100 96020100
  geo_blk_bld bayes_fh_heavy migrated_intl_12m            1 96030001
  geo_blk_bld bayes_fh_heavy migrated_intl_12m           57 96030057
  geo_blk_bld bayes_fh_heavy migrated_intl_12m          100 96030100
  geo_blk_bld bayes_fh_heavy pers_migrated_12m            1 96040001
  geo_blk_bld bayes_fh_heavy pers_migrated_12m           57 96040057
  geo_blk_bld bayes_fh_heavy pers_migrated_12m          100 96040100
  geo_blk_bld bayes_fh_heavy pers_migrated_internal_12m   1 96050001
  geo_blk_bld bayes_fh_heavy pers_migrated_internal_12m  57 96050057
  geo_blk_bld bayes_fh_heavy pers_migrated_internal_12m 100 96050100
  geo_blk_bld bayes_fh_heavy pers_migrated_intl_12m       1 96060001
  geo_blk_bld bayes_fh_heavy pers_migrated_intl_12m      57 96060057
  geo_blk_bld bayes_fh_heavy pers_migrated_intl_12m     100 96060100
  geo_blk_bld xgb_sub        migrated_last_12m            1 97010001
  geo_blk_bld xgb_sub        migrated_last_12m           57 97010057
  geo_blk_bld xgb_sub        migrated_last_12m          100 97010100
  geo_blk_bld xgb_sub        migrated_internal_12m        1 97020001
  geo_blk_bld xgb_sub        migrated_internal_12m       57 97020057
  geo_blk_bld xgb_sub        migrated_internal_12m      100 97020100
  geo_blk_bld xgb_sub        migrated_intl_12m            1 97030001
  geo_blk_bld xgb_sub        migrated_intl_12m           57 97030057
  geo_blk_bld xgb_sub        migrated_intl_12m          100 97030100
  geo_blk_bld xgb_sub        pers_migrated_12m            1 97040001
  geo_blk_bld xgb_sub        pers_migrated_12m           57 97040057
  geo_blk_bld xgb_sub        pers_migrated_12m          100 97040100
  geo_blk_bld xgb_sub        pers_migrated_internal_12m   1 97050001
  geo_blk_bld xgb_sub        pers_migrated_internal_12m  57 97050057
  geo_blk_bld xgb_sub        pers_migrated_internal_12m 100 97050100
  geo_blk_bld xgb_sub        pers_migrated_intl_12m       1 97060001
  geo_blk_bld xgb_sub        pers_migrated_intl_12m      57 97060057
  geo_blk_bld xgb_sub        pers_migrated_intl_12m     100 97060100
  geo_blk_bld megb_sub       migrated_last_12m            1 98010001
  geo_blk_bld megb_sub       migrated_last_12m           57 98010057
  geo_blk_bld megb_sub       migrated_last_12m          100 98010100
  geo_blk_bld megb_sub       migrated_internal_12m        1 98020001
  geo_blk_bld megb_sub       migrated_internal_12m       57 98020057
  geo_blk_bld megb_sub       migrated_internal_12m      100 98020100
  geo_blk_bld megb_sub       migrated_intl_12m            1 98030001
  geo_blk_bld megb_sub       migrated_intl_12m           57 98030057
  geo_blk_bld megb_sub       migrated_intl_12m          100 98030100
  geo_blk_bld megb_sub       pers_migrated_12m            1 98040001
  geo_blk_bld megb_sub       pers_migrated_12m           57 98040057
  geo_blk_bld megb_sub       pers_migrated_12m          100 98040100
  geo_blk_bld megb_sub       pers_migrated_internal_12m   1 98050001
  geo_blk_bld megb_sub       pers_migrated_internal_12m  57 98050057
  geo_blk_bld megb_sub       pers_migrated_internal_12m 100 98050100
  geo_blk_bld megb_sub       pers_migrated_intl_12m       1 98060001
  geo_blk_bld megb_sub       pers_migrated_intl_12m      57 98060057
  geo_blk_bld megb_sub       pers_migrated_intl_12m     100 98060100
")
}

test_that("sae_fit_seed reproduces the pinned seeds for both configurations", {
  pin <- pinned_fit_seeds()
  expect_identical(nrow(pin), 288L)
  expect_identical(anyDuplicated(pin[c("config", "method", "outcome", "rep")]), 0L)
  expect_identical(
    sae_fit_seed(pin$config, pin$method, pin$outcome, pin$rep),
    pin$seed)
})

test_that("a seed is the digit-partitioned sum of the index of each part", {
  # 1e7 * 2 + 1e6 * 7 + 1e4 * 3 + 57
  expect_identical(sae_fit_seed("geo_bld", "xgb_sub", "migrated_intl_12m", 57L),
                   27030057L)
  # 1e7 * 9 + 1e6 * 8 + 1e4 * 6 + 100
  expect_identical(
    sae_fit_seed("geo_blk_bld", "megb_sub", "pers_migrated_intl_12m", 100L),
    98060100L)
  expect_identical(sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 1L), 21010001L)
})

test_that("the two configurations differ in the seed by their config index only", {
  # geo_bld has config index 2 and geo_blk_bld index 9: 7e7 apart, all else equal
  for (m in names(sae_seed_method_index)) {
    for (o in names(sae_seed_outcome_index)) {
      expect_identical(
        sae_fit_seed("geo_blk_bld", m, o, 7L) - sae_fit_seed("geo_bld", m, o, 7L),
        70000000L, info = paste(m, o))
    }
  }
  expect_identical(unname(sae_seed_config_index), c(2L, 9L))
})

test_that("sae_fit_seed is vectorised and recycles its arguments", {
  expect_identical(
    sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 1:3),
    sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 1L) + 0:2)
  expect_identical(
    sae_fit_seed("geo_bld", c("fh", "fh_arcsin"), "migrated_last_12m", 1L),
    c(21010001L, 22010001L))
})

test_that("an unknown config, method or outcome is refused by name", {
  expect_error(sae_fit_seed("geo_x", "fh", "migrated_last_12m", 1L),
               "Unknown config: geo_x", fixed = TRUE)
  # the parent configurations are not seeded here, only the with-buildings ones
  expect_error(sae_fit_seed("geo", "fh", "migrated_last_12m", 1L),
               "Unknown config: geo", fixed = TRUE)
  expect_error(sae_fit_seed("geo_blk", "fh", "migrated_last_12m", 1L),
               "Unknown config: geo_blk", fixed = TRUE)
  expect_error(sae_fit_seed("geo_bld", "xgb", "migrated_last_12m", 1L),
               "Unknown method: xgb", fixed = TRUE)
  expect_error(sae_fit_seed("geo_bld", "fh", "not_an_outcome", 1L),
               "Unknown outcome: not_an_outcome", fixed = TRUE)
  expect_error(sae_fit_seed(c("geo_bld", "geo_x"), "fh", "migrated_last_12m", 1L),
               "geo_x", fixed = TRUE)
})

test_that("a replicate outside 1-9999 is refused", {
  expect_error(sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 0L))
  expect_error(sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 10000L))
  expect_error(sae_fit_seed("geo_bld", "fh", "migrated_last_12m", NA_integer_))
})

test_that("the outcome index follows the outcome registry", {
  expect_identical(names(sae_seed_outcome_index), sae_outcomes$outcome)
  expect_identical(unname(sae_seed_outcome_index), seq_along(sae_outcomes$outcome))
})

test_that("the method index names every registered sub-area method", {
  expect_true(all(sae_methods_subarea %in% names(sae_seed_method_index)))
  expect_identical(anyDuplicated(sae_seed_method_index), 0L)
})

test_that("method_args forwards the seed to the Bayesian methods", {
  d  <- data.frame(area_id = c("05001", "05002"), outcome = "migrated_last_12m",
                   theta_direct = c(0.1, 0.2))
  X  <- data.frame(municipality = c("05001", "05002"), x1 = c(1, 2))
  a  <- method_args("bayes_fh_beta", direct = d, aux_in = X, aux_oos = X[0, ],
                    seed = 123L)
  expect_identical(a$seed, 123L)
  h  <- method_args("bayes_fh_heavy", direct = d, aux_in = X, aux_oos = X[0, ],
                    seed = 123L)
  expect_identical(h$seed, 123L)
  # no seed given: the heavy sampler keeps its own NULL default
  h0 <- method_args("bayes_fh_heavy", direct = d, aux_in = X, aux_oos = X[0, ])
  expect_null(h0$seed)
})
