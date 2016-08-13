# Configuration options for halos

import numpy as np

DT_COORDS = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])    # coordinates data type
DT_COORDSwID = np.dtype([('id', np.int64), ('x', np.float32), ('y', np.float32), ('z', np.float32)])

DT_HEADER = np.dtype([('magic', np.uint64), \
                      ('version', np.int64), \
                      ('num_files', np.int64), \
                      ('file_id', np.int64), \
                      ('snapshot', np.int64), \
                      ('format_group_data', np.int64), \
                      ('format_part_data', np.int64), \
                      ('group_type', np.int64), \
                      ('ngroups', np.int64), \
                      ('ngroups_total', np.int64), \
                      ('npart', np.int64), \
                      ('npart_total', np.int64), \
                      ('npart_orig', np.int64), \
                      ('max_npart', np.int64), \
                      ('max_npart_total', np.int64), \
                      ('min_group_part', np.int64), \
                      ('valid_part_ids', np.int64), \
                      ('linkinglength', np.float64), \
                      ('overdensity', np.float64), \
                      ('time', np.float64), \
                      ('redshift', np.float64), \
                      ('box_size', np.float64), \
                      ('box_min_x', np.float64), \
                      ('box_min_y', np.float64), \
                      ('box_min_z', np.float64), \
                      ('bounds_xmin', np.float64), \
                      ('bounds_xmax', np.float64), \
                      ('bounds_ymin', np.float64), \
                      ('bounds_ymax', np.float64), \
                      ('bounds_zmin', np.float64), \
                      ('bounds_zmax', np.float64), \
                      ('part_mass', np.float64), \
                      ('Omega0', np.float64), \
                      ('OmegaLambda', np.float64), \
                      ('Hubble0', np.float64), \
                      ('GravConst', np.float64)])

DT_GROUPS = np.dtype([('id', np.int64), \
                      ('parent_id', np.int64), \
                      ('npart', np.uint64), \
                      ('npart_self', np.uint64), \
                      ('radius', np.float32), \
                      ('mass', np.float32), \
                      ('x', np.float32), \
                      ('y', np.float32), \
                      ('z', np.float32), \
                      ('vx', np.float32), \
                      ('vy', np.float32), \
                      ('vz', np.float32), \
                      ('vmax', np.float32), \
                      ('rvmax', np.float32)])

DT_PARTICLES = np.dtype([('id', np.int64), \
                         ('x', np.float32), \
                         ('y', np.float32), \
                         ('z', np.float32), \
                         ('vx', np.float32), \
                         ('vy', np.float32), \
                         ('vz', np.float32)])

_default_ascii_settings = {'dtype': DT_COORDS,
                           'comments':'#',
                           'delimiter':None,
                           'skip_header': 1,
                           'skip_footer':0,   #   cannot be used together with 'max_rows' argument
                           'converters':None,
                           'missing_values':None,
                           'filling_values':None,
                           # 'usecols': (1, 2, 3),
                           'names': True,
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

# fields in header that must be different
EXCLUSIVE_HEADER_FIELDS = ('file_id',\
                           )
