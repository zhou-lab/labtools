
def pd_print_n(x, n):
    import pandas as pd
    pd.set_option('display.max_rows', n)
    print(x)
    pd.reset_option('display.max_rows')
    
def pd_print_full(x):
    import pandas as pd
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')
