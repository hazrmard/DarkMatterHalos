# Configuration options for halos

import numpy as np

COORDS = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])    # coordinates data type
COORDSwID = np.dtype([('id', np.int64), ('x', np.float32), ('y', np.float32), ('z', np.float32)])

_default_ascii_settings = {'dtype': COORDS,
                           'comments':'#',
                           'delimiter':None,
                           'skip_header': 1,
                           'skip_footer':0,   #   cannot be used together with 'max_rows' argument
                           'converters':None,
                           'missing_values':None,
                           'filling_values':None,
                           'usecols': (1, 2, 3),
                           'names':None,
                           'excludelist':None,
                           'deletechars':None,
                           'replace_space':'_',
                           'autostrip':False,
                           'case_sensitive':True,
                           'defaultfmt':'f%i',
                           'unpack':None,
                           'usemask':False,
                           'loose':True,
                           'invalid_raise':True,
                           # 'max_rows':200
                           }
