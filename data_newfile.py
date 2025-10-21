from data_getvals import MFG_out


''' Note: first iterate over masses, not run sets'''

# for run in ('A','F'):
#     for iso in [162]:
for iso in [164]:
    for run in ('F'):
        dim = MFG_out('dim',**{'iso':iso,'run':run}).kv_init
