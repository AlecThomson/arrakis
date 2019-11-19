"""Generate parsets and sbatch scripts to run selavy on RACS single beam images.
Assumes pipeline naming i.e. image.i.SBXXXX.cont.RACS_test4_1.05_XXXXXXXX.beamXX.taylor.X.restored.fits.
Output images are not written to conserve disk space.
Each beam in the given field will be run in the same job serially.
"""
import re
from pathlib import Path
import argparse
import textwrap
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("parset_template", help="Selavy parset template file. Must contain format placeholders {{image_base}} and {{sbid}}.")
parser.add_argument("sbatch_template", help="Selavy sbatch template file. Must contain format placeholders {{clusters}}, {{account}}, {{time}}, {{job_name}} and {{parset}}.")
parser.add_argument("field_dir", type=Path, help="Directory containing input FITS images for the same RACS field.")
parser.add_argument("--submit", action="store_true", help="Submit the job to the cluster.")
parser.add_argument("--clusters", default="magnus", help="Cluster to run the job. Default: magnus.")
parser.add_argument("--account", default="ja3", help="Account to run the job under. Default: ja3.")
parser.add_argument("--time", default="08:00:00", help="Time limit for the job. Default 08:00:00.")

args = parser.parse_args()

# read the templates
with open(args.parset_template) as parset_template_fd:
    print("reading parset template file: {}".format(args.parset_template))
    parset_template = parset_template_fd.read()

images = args.field_dir.glob("image.i.*.taylor.0.restored.fits")
field = args.field_dir.name

# make field dir
output_field_dir = Path("fields") / Path(field)
output_field_dir.mkdir(parents=True, exist_ok=True)

# format the sbatch script header
script_header = """
    #!/bin/bash -l
    #SBATCH --partition=workq
    #SBATCH --clusters={clusters}
    #SBATCH --account={account}
    #SBATCH --export=NONE
    #SBATCH --ntasks=19
    #SBATCH --ntasks-per-node=20
    #SBATCH --time={time}
    #SBATCH --job-name={job_name}

    module load slurm
    module unload python
    module load python
    module use /group/askap/modulefiles
    module load askapdata
    module unload askapsoft
    module load askapsoft/0.24.4
    module unload askappipeline
    module load askappipeline/0.24.4
""".format(
    clusters=args.clusters,
    account=args.account,
    time=args.time,
    job_name="selavy_{field}".format(field=field)
)
script = textwrap.dedent(script_header)[1:]  # remove intentation and first blank line

for image in images:
    sbid, beam = re.match(r"image.i.SB(\d{4}).cont.RACS_test4_1.05_.+?.beam(\d{2}).taylor.\d.restored", image.stem).groups()
    image_base = image.stem

    # write parset
    parset = parset_template.format(image=image.resolve(strict=True), image_base=image_base, sbid=sbid)
    parset_file = output_field_dir / Path("selavy_{}.in".format(image_base))
    with parset_file.open('w') as parset_file_fd:
        parset_file_fd.write(parset)

    srun_cmd = "srun --export=ALL --ntasks=19 --ntasks-per-node=20 selavy -c {parset} >> selavy_${{SLURM_JOB_ID}}.log".format(parset=parset_file.name)
    script += srun_cmd + "\n"

script_file = output_field_dir / Path("selavy_{}.sbatch".format(field))
with script_file.open('w') as script_file_fd:
    script_file_fd.write(script)

# submit job
if args.submit:
    subprocess.Popen("sbatch {}".format(script_file.name), shell=True, cwd=output_field_dir)

