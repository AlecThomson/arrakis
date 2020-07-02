#!/group/askap/askapcli/current/bin/python

from __future__ import division, print_function
import os
import subprocess
from subprocess import Popen, PIPE
import shlex
from glob import glob
import multiprocessing as mp
from functools import partial

def worker(src, info_dict, scriptdir):
    if src == 'src%d':
        return
    else:
        field = info_dict[src]['field_name']
        pol_lst = info_dict[src]['pol_axis']
        pol_ax_rot = float(pol_lst[1])
        foot_rot = float(info_dict['src%d']['footprint']['rotation'])
        rotation = pol_ax_rot+foot_rot
        direction = info_dict[src]['field_direction']
        cmd = "footprint calculate -d " + \
            str(direction[0]) + "," + str(direction[1]) + \
                " -r " + str(rotation) + " " + \
                    str(info_dict['src%d']['footprint']['name']) + \
                        " -p " + str(info_dict['src%d']['footprint']['pitch'])
        p = Popen(shlex.split(cmd), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        outfile = scriptdir + "/../askap_surveys/RACS/admin/epoch_0/" + "beam_sep-" + field + '.csv'
        csv = output.replace("  ", ",").replace(" -", ",-").replace(" ", "").replace("(", "").replace(")", "")
        csv = csv.splitlines()
        new_csv = []
        for line in csv:
            line +=  "," + str(direction[0]) + "," + str(direction[1])
            new_csv.append(line)
        
        header = ["BEAM,DELTA_RA,DELTA_DEC,BEAM_RA,BEAM_DEC,FOOTPRINT_RA,FOOTPRINT_DEC"]
        new_csv = header + new_csv
        with open(outfile, 'w') as f:
            for line in new_csv:
                f.write("%s\n" % line)
        print('Written to', outfile)

def main():
    # Use ASKAPcli to get beam separations for PB correction
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    schedfile = scriptdir + "/../askap_surveys/RACS/db-inputs/RACS_SkyScan_10.parset"


    # Get the SBID info
    f = open(schedfile, 'r')
    info = f.readlines() 
    keys, vals = [], []
    for line in info: 
        temp = line.strip() 
        key = temp.split()[0].split('.')

        keys.append(key)

        idx = temp.find('=') 
        val = temp[idx+2:] 
        if temp[-1] == '=': 
            val = None
        elif val[0] == '[': 
            val = val.strip('][').split(', ')
        vals.append(val)    

    srcs = []
    for key, val in zip(keys, vals):
        if key[0] == 'craft':
            continue
        elif key[1] == 'target':
            srcs.append(key[2])

    # Store data in dict stucture
    srcs = sorted(set(srcs))
    info_dict = {}
    for src in srcs:
        info_dict.update({src: {}})

    info_dict['src%d']['footprint'] = {}

    for key, val in zip(keys, vals):
        if key[0] == 'craft':
            continue
        elif key[1] == 'target':
            src = key[2]
            if len(key[3:]) == 1:
                info_dict[src].update(
                    {
                        key[3]: val
                    }
                )
            elif key[3] == 'footprint':
                info_dict[src]['footprint'].update(
                    {
                        key[4]: val
                    }
                )
    count = 0
    worker_part = partial(worker, info_dict=info_dict, scriptdir=scriptdir)
    #map(worker_part, info_dict.keys())
    pool = mp.Pool(mp.cpu_count()//2)
    pool.map(worker_part, info_dict.keys())
    pool.close()


if __name__ == "__main__":
    main()