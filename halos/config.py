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

# fields in header that have to be same for 2 bgc2 files to merge
MANDATORY_HEADER_FIELDS = ('magic',\
                           'version',\
                           'num_files',\
                           'snapshot',\
                           'format_group_data',\
                           'format_part_data',\
                           'group_type',\
                           'ngroups_total',\
                           'npart_total',\
                           'max_npart_total',\
                           'min_group_part',\
                           'valid_part_ids',\
                           'linkinglength',\
                           'overdensity',\
                           'time',\
                           'redshift',\
                           'box_size',\
                           'box_min_x',\
                           'box_min_y',\
                           'box_min_z',\
                           'part_mass',\
                           'Omega0',\
                           'OmegaLambda',\
                           'Hubble0',\
                           'GravConst',)
